#!/usr/bin/env perl

=pod
###### 2021-01-04 13:47 %Yong Li%

This script is to correct frameshift based on SPALN alignment and vcf file.

#1. Align protein sequences to genome using SPALN.
#2. Identify frameshift regions from SPALN results.
#3. Search frameshift regions and its variants from VCF file (pre-generated), then correct frameshift.
=cut

use strict;
use warnings;
use DB_File;
use DB_File::Lock;
use Fcntl qw(:flock O_RDWR O_CREAT O_RDONLY O_WRONLY);
use Time::Interval;
use Parallel::ForkManager;

my $thisScript = $0;
$thisScript = $1 if ( $0 =~ /\\([^\\]+)$/ || $0 =~ /\/([^\/]+)$/ );

if ( @ARGV < 5 || $ARGV[0] =~ /^-{1,2}h/ ) {
    print "\n$0:\n";
    print "Usage: ./$thisScript genome.fa genome_spaln_db protein.fa variant.vcf out_prefix [num_threads].\n";
    print "This script is to correct frameshift based on SPALN alignment and vcf file.\n\n";
    exit(1);
}

print "\nRunning $thisScript" . ( defined( $ARGV[0] ) ? " for $ARGV[0] " : " " ) . "now ...\n";
my $start_time = localtime(time);

my $num_threads = 60;
$num_threads = $ARGV[5] if defined( $ARGV[5] );

my $chroms_fa   = shift @ARGV;
my $chroms_db   = shift @ARGV;
my $protein_fa  = shift @ARGV;
my $variant_vcf = shift @ARGV;
my $out_prefix  = shift @ARGV;

my $tie_file      = "genome_DBFile_tmp";
my $out_spaln_gff = "$out_prefix.spaln.gff3";
my $out_spaln_log = "$out_prefix.spaln.log";
my $out_corr_fa   = "$out_prefix.frameshift_fix.fa";
my $log_file      = "$out_prefix.fix.log";
my $log           = 1;
my $verbose       = 0;

### find frameshifts
# spaln -W -KP -t100 genome_test.fa && mv genome_test.??? /path/to/spaln/seqdb/
# spaln -Q6 -S3 -H20 -M -O0 -t100 -o test.gff3 -dgenome_test test.pep.fa 2>&1 | tee test.spaln.log
# transcript_1 Chr01 FrameShift: 39258367 -1
# transcript_1 Chr01 FrameShift: 39258317 2
print "\nRunning spaln now ...\n";
if ( !-e $out_spaln_log ) {
    system( "spaln -Q6 -S3 -H20 -M -O0 -t100 -o $out_spaln_gff -d$chroms_db $protein_fa" . ' 2>&1 | tee ' . $out_spaln_log ) >> 8
        and die "Error running spaln: exit code $?";
}
else {
    print "Spaln result file $out_spaln_log exists, skipping running spaln.\n";
}
print "    ... spaln running finished.\n";

print "\nReading in spaln result now ...\n";
my (%frameshifts);
open( SPALN_OUTPUT, $out_spaln_log ) or die "Cannot open spaln output file $out_spaln_log for reading: $!\n";
while ( my $line = <SPALN_OUTPUT> ) {
    chomp $line;
    next unless ( $line =~ /FrameShift/ );
    my @line = split /\s/, $line;
    $frameshifts{ $line[1] }->{ $line[3] } = $line[4];
    print "Frameshift: Seq:$line[1]\tPos:$line[3]\tType:$line[4]\n" if $verbose;
}
close(SPALN_OUTPUT);
print "    ... spaln result reading finished.\n";

### read in genome sequences
print "\nReading in genome now ...\n";
my %chroms;
unlink $tie_file;
tie( %chroms, 'DB_File::Lock', $tie_file, O_CREAT | O_RDWR, 0666, $DB_BTREE, "write" ) || die "Cannot open $tie_file for tie: $!\n";
open( CHROMS_FASTA, $chroms_fa ) or die "Cannot open $chroms_fa for reading: $!\n";
$/ = "\n>";    # record separator for FASTA
while (<CHROMS_FASTA>) {
    my ( $desc, $seq ) = />?(.*?)\n(.*)/s;
    $desc =~ /^([^\s]+)/;
    my $chrom_id = $1;
    $seq =~ tr/ \n//d;
    $seq =~ s/>$//;
    $chroms{$chrom_id} = $seq;    # for editing later
    print "Read in genome: " . $chrom_id . "\n" . $seq . "\n" if $verbose;
}
close(CHROMS_FASTA);
$/ = "\n";
untie %chroms;
print "    ... genome reading finished.\n";

### read in variants
# Chr00	77	.	G	A
# Chr00	86	.	CT	C
print "\nReading in VCF now ...\n";
open( VCF, $variant_vcf ) or die "Cannot open $variant_vcf for reading: $!\n";
my %var;
while (<VCF>) {
    chomp;
    my $variant = $_;
    next if ( $variant =~ /^#/ );                                                   # skip comment lines
    my @variant = split /\t/, $variant;
    next if ( $variant[4] eq '.' );                                                 # skip non-variant position
    $var{ $variant[0] }->{ $variant[1] }[0] = $variant[3];
    $var{ $variant[0] }->{ $variant[1] }[1] = $variant[4];
    $var{ $variant[0] }->{ $variant[1] }[1] = $1 if $variant[4] =~ /^([A-Z]+),/;    # multiple allele, use the 1st variant
    print "$variant[0]\t$variant[1]\t$variant[3]\t$variant[4]\n" if $verbose;
}
close(VCF);
print "    ... VCF reading finished.\n";

### fix frameshift
print "\nFixing frameshifts now ...\n";
open( LOG, ">$log_file" ) || die "Cannot open $log_file: $!\n" if $log;
my $pm = new Parallel::ForkManager($num_threads);
for my $chrom_id ( keys %frameshifts ) {
    my $pid = $pm->start and next;                                                  # Forks and returns the pid for the child
    tie( %chroms, "DB_File::Lock", $tie_file, O_RDONLY, 0666, $DB_BTREE, "read" ) || die "Cannot open $tie_file for tie: $!\n";
    my $seq = $chroms{$chrom_id};
    untie %chroms;
    print "Fix: " . $chrom_id . "\n" . $seq . "\n" if $verbose;

    # Editing only works with the original coordinates if you process from 3' to 5' on the positive ref strand
    for my $pos ( sort { $b <=> $a } keys %{ $frameshifts{$chrom_id} } ) {
        my $shift_type = $frameshifts{$chrom_id}->{$pos};

        # search 5 positions (up-stream 1 nt + 3 nt (a codon) + down-stream 1 nt), 1-based for seq, 0-based for perl string.
        # +1 frameshift in VCF: GA G, -1 frameshift in VCF: T TC, so it may show at adjacent nt.
        # correct frameshift from 3' to 5'.
        for ( my $i = $pos + 3; $i >= $pos - 1; $i-- ) {
            next unless exists( $var{$chrom_id}->{$i}[0] );
            print "  var: $chrom_id: $i " . $var{$chrom_id}->{$i}[0] . " -> " . $var{$chrom_id}->{$i}[1] . "\n" if $verbose;
            if ( $shift_type == 0 and length( $var{$chrom_id}->{$i}[0] ) == length( $var{$chrom_id}->{$i}[1] ) ) {    # stop codon gain
                print LOG "$chrom_id\t$i\t" . substr( $seq, $i - 1, length( $var{$chrom_id}->{$i}[0] ) ) . "\t" . $var{$chrom_id}->{$i}[1] . "\n" if $log;
                substr( $seq, $i - 1, length( $var{$chrom_id}->{$i}[0] ) ) = $var{$chrom_id}->{$i}[1];
            }
            elsif ( $shift_type != 0 and length( $var{$chrom_id}->{$i}[0] ) != length( $var{$chrom_id}->{$i}[1] ) ) {    # frameshift
                print LOG "$chrom_id\t$i\t" . substr( $seq, $i - 1, length( $var{$chrom_id}->{$i}[0] ) ) . "\t" . $var{$chrom_id}->{$i}[1] . "\n" if $log;
                substr( $seq, $i - 1, length( $var{$chrom_id}->{$i}[0] ) ) = $var{$chrom_id}->{$i}[1];
                last;                                                                                                    # quit if changed any position.
            }
        }
    }

    # save corrected seq
    tie( %chroms, 'DB_File::Lock', $tie_file, O_CREAT | O_WRONLY, 0666, $DB_BTREE, "write" ) || die "Cannot open $tie_file for tie: $!\n";
    $chroms{$chrom_id} = $seq;
    untie %chroms;

    # Terminates the child process
    $pm->finish;
}
$pm->wait_all_children;
close(LOG) if $log;
print "    ... frameshifts fixing finished.\n";

### write out corrected fasta sequences
print "\nWriting out result now ...\n";
tie( %chroms, "DB_File::Lock", $tie_file, O_RDONLY, 0666, $DB_BTREE, "read" ) || die "Cannot open $tie_file for tie: $!\n";
open( FASTA_OUT, ">$out_corr_fa" ) or die "Cannot open $out_corr_fa for writing: $!\n";
for my $contig_id ( sort keys %chroms ) {
    my $seq = $chroms{$contig_id};
    $seq =~ s/(.{1,100})/$1\n/g;
    print FASTA_OUT ">", $contig_id, "\n", $seq;
}
close(FASTA_OUT);
untie %chroms;
unlink $tie_file;
unlink $tie_file . ".lock";
print "    ... result writing finished.\n";

my $end_time     = localtime(time);
my $elapsed_time = getInterval( $start_time, $end_time, 'String' => 1 );
print "    ... $thisScript running finished, used time: $elapsed_time.\n\n";
