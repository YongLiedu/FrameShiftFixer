#!/usr/bin/env perl

=pod
###### 2021-01-04 13:47 %Yong Li%

This script is to correct frameshift based on Diamond BLASTX search and vcf file.

#1. Align protein sequences to genome using Diamond (BLASTX search).
#2. Identify frameshift regions from Diamond results.
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

if ( @ARGV < 4 || $ARGV[0] =~ /^-{1,2}h/ ) {
    print "\n$0:\n";
    print "Usage: ./$thisScript genome.fa protein.dmnd variant.vcf out_prefix [num_threads].\n";
    print "This script is to correct frameshift based on Diamond BLASTX search and vcf file.\n\n";
    exit(1);
}

print "\nRunning $thisScript" . ( defined( $ARGV[0] ) ? " for $ARGV[0] " : " " ) . "now ...\n";
my $start_time = localtime(time);

my $matrix         = "PAM30";    # better suited for closely related, short matches
my $evalue_diamond = 1e-5;       # evalue threshold for diamond output
my $evalue_corr    = 1e-6;       # evalue threshold for frameshift search
my $num_threads    = 60;
$num_threads = $ARGV[4] if defined( $ARGV[4] );

my $chroms_fa    = shift @ARGV;
my $protein_dmnd = shift @ARGV;
my $variant_vcf  = shift @ARGV;
my $out_prefix   = shift @ARGV;

my $tie_file    = "genome_DBFile_tmp";
my $out_diamond = "$out_prefix.diamond.txt";
my $out_corr_fa = "$out_prefix.frameshift_fix.fa";
my $log_file    = "$out_prefix.fix.log";
my $log         = 1;
my $verbose     = 0;

### find frameshifts
# diamond makedb --in test_aa.fa -d test_aa
# diamond blastx -d $protein_dmnd -q $chroms_fa -o $out_diamond -e $evalue_diamond -p $num_threads --matrix $matrix -F 10 --gapopen 8 --gapextend 1 --outfmt 6 qseqid sseqid qstart qend sstart send qframe btop evalue
# ori_luc_nt	LUC_aa	36	1	1	12	-1	12         # complete match
# LUC_nt	LUC_aa	36	1	1	12	-1	8*K3           # stop codon gain
# rep_nt1	rep_aa	1	108	1	36	1	2/-PT13\-20    # -1 nt frameshift and +1 nt frameshift
# rep_nt2	rep_aa	1	106	1	36	1	17-H\-18       # -2 nt frameshift
# rep_nt3	rep_aa	1	110	1	36	1	17H-/-19       # +2 nt frameshift
print "\nRunning DIAMOND now ...\n";
if ( !-e $out_diamond ) {
    system(   "diamond blastx -d $protein_dmnd -q $chroms_fa -o $out_diamond "
            . "-e $evalue_diamond -p $num_threads --matrix $matrix -F 10 --gapopen 8 --gapextend 1 --masking 0 "
            . "-b12 -c1 --ultra-sensitive --max-hsps 100 -k 100000 --quiet "
            . "--outfmt 6 qseqid sseqid qstart qend sstart send qframe btop evalue" ) >> 8
        and die "Error running diamond blastx: exit code $?";
}
else {
    print "Diamond result file $out_diamond exists, skipping running Diamond.\n";
}
print "    ... DIAMOND running finished.\n";

print "\nReading in DIAMOND result now ...\n";
my (%frameshifts);
open( DIAMOND_OUTPUT, $out_diamond ) or die "Cannot open Diamond output file $out_diamond for reading: $!\n";
while (<DIAMOND_OUTPUT>) {
    chomp;
    my ( $qseqid, $sseqid, $qstart, $qend, $sstart, $send, $qframe, $btop, $evalue ) = split /\t/, $_;
    next if $evalue > $evalue_corr;

    # The BTOP (BLAST traceback operation) will have a / if there is a -1 frameshift,
    # and a \ if there is a +1 frameshift, and a * if there is a gained stop codon.
    next unless ( $btop =~ /\/|\\|\*/ );
    my @btop_parts = split /(\/\-[A-Z][A-Z]|\-[A-Z]\\\-|[A-Z]\-\/\-|\*[A-Z]|\\\-)/, $btop;
    my $offset     = 0;
    for ( my $i = 0; $i < @btop_parts - 1; $i += 2 ) {
        while ( $btop_parts[$i] =~ /(\d+|..)/g ) {    # match number or mismatch pair
            my $move = $1;
            if ( $move =~ /^\d+$/ ) {                 # match number of aa
                $offset += $move;
            }
            elsif ( $move =~ /^[A-Z]/ ) {             # mismatch and query/DNA uses one aa
                $offset++;
            }
        }

        # frameshift's location (1-based) in query DNA
        my $frameshift_location = ( $qframe > 0 ) ? ( 3 * $offset + $qstart ) : ( $qstart - ( 3 * $offset ) - 2 );
        my $frameshift_type     = "";
        if ( $btop_parts[ $i + 1 ] =~ /^\// ) {
            $frameshift_type = -1;
        }
        elsif ( $btop_parts[ $i + 1 ] =~ /^\\/ ) {
            $frameshift_type = 1;
        }
        elsif ( $btop_parts[ $i + 1 ] =~ /^\*/ ) {
            $frameshift_type = 0;
        }
        elsif ( $btop_parts[ $i + 1 ] =~ /\\/ ) {
            $frameshift_type = -2;
        }
        elsif ( $btop_parts[ $i + 1 ] =~ /\// ) {
            $frameshift_type = 2;
        }
        $offset += $frameshift_type;    # correct offset for next variant.
        $frameshifts{$qseqid}->{$frameshift_location} = $frameshift_type;
        print "Frameshift: Seq:$qseqid\tPos:$frameshift_location\tType:$frameshift_type\tPreBTOP:" . $btop_parts[$i] . "\tShiftBTOP:" . $btop_parts[ $i + 1 ] . "\n" if $verbose;
    }
}
close(DIAMOND_OUTPUT);
print "    ... DIAMOND result reading finished.\n";

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
