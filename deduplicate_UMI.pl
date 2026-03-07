#!/usr/bin/perl -w

# Deduplicate FASTQ sequences using UMIs
#
# DedupUMI removes exact duplicate FASTQ read pairs using UMI sequences without requiring alignment to a reference.
#
# Duplicates are identified using concatenation of R1, R2 and UMI sequences.
# For duplicate molecules the read pair with the highest total base quality score is retained.
#
# ## Why deduplicate?
#
# Exact duplicates (incl UMI) are typically the result of over-amplification / too many amplification cycles in sequence library prepping 'PCR'.  
#
# ## Duplicate definition
#
# Reads are considered duplicates when the following combination is identical:
# R1 sequence + R2 sequence + UMI
# *For the older R3 system this becomes: R1 sequence + R2 sequence + R3 sequence*
#
# ## Memory usage
#
# Deduplication stores unique read keys in memory.
# Memory usage scales with the number of unique molecules in the dataset.
# In practice MAXIMUM memory usage is approximately:
#  - 2 × size of the R1 sequence file
#  - + optional R3 UMI file
#
# ## Input and output
#
# NEW system R1 R2 UMI in header:
#  - Input: R1 and R2 FASTQ readset. UMI should be present in the header of at least R1!
#  - Output: R1 and R2 filtered for EXACT duplicates based on concatenated seq of R1 R2 keeping highest TOTAL qual score.
#
# Older system R1 R2 R3 UMI in R2:
#  - Input: R1 R2 and R3 FASTQ readset. For the best total quality to work the UMI should be present in R2!
#  - Output: R1 R2 R3 filtered for EXACT duplicates based on concatenated seq of R1 R2 R3 keeping highest TOTAL qual score.
#
# ## Requirements
#
#  - Unix system with `zcat` to allow reading/writing of gzipped files
#  - FASTQ files in fixed 4-line format having read-pairs in separate files (do not use interleaved files)
#
# Optional but recommended:
#  - pigz (parallel gzip) for faster compression
#
#
# ## Important notes
#
# - It writes out the last sequences found of a duplicate set having the highest TOTAL qualityscore.
# - Input FASTQ files must follow the standard 4-line FASTQ format, starting at the first record.
# - Sequences in R1 R2 (and R3) should be in same order and NOT INTERLEAVED!!
#
# ## Author
# a.bossers@uu.nl // alex.bossers@wur.nl
# @author: a.bossers@uu.nl // alex.bossers@wur.nl
#
# Latest version ad documentation: https://github.com/alxelerator/dedupUMI
#
#
# Disclaimer:  
#         Script is provided AS IS under GPL-3.
#         We did our best to verify that the results are legitimate. However, the output should be considered erroneous, so you should check your results.
#         The authors nor their institutions/employers are in any way direct or indirect responsible for the direct or indirect damages caused by using this script.
#         Use at your own responsibility.
#

my $versionno     = "1.4";
my $versiondate   = "2026-03-07";

# Version history
#       1.4      07-03-2026  Integrated R2 and R3 system to remove all duplicated code.
#       1.3      05-03-2026  Improved file out of sync detection and removed duplicated code (file2 or file3 system)
#       1.2b     27-02-2026  Fixed incoming reads counter.
#       1.2      25-02-2026  Added an option to write out the counts incoming and written to file (APPEND!)
#       1.1b     25-02-2026  Verbosity printing to console bit organised
#       1.1      29-08-2023  Add option to use UMI in fastq headers instead of R3-system.
#       0.9=1.0  12-10-2022  Added quality check (Wouter + streamline Alexc)
#       0.8      15-09-2022  first working concept 
#       0.1      xx-xx-2010  template

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
#use IO::Compress::Gzip qw(gzip $GzipError) ;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use POSIX qw/strftime/;
use Data::Dumper;

#handle cmd line options
my ( $input_fastq1, $input_fastq2, $input_fastq3, $output_fastq1, $output_fastq2, $output_fastq3, $output_counts, $noUMI, $verbose, $help, $version );

GetOptions ('input-fastq1=s'   => \$input_fastq1,
            'input-fastq2=s'   => \$input_fastq2,
            'input-fastq3=s'   => \$input_fastq3,
            'output-fastq1=s'  => \$output_fastq1,
            'output-fastq2=s'  => \$output_fastq2,
            'output-fastq3=s'  => \$output_fastq3,
            'output-counts=s'  => \$output_counts,
            'noUMI'            => \$noUMI,
            'verbose'          => \$verbose,
            'help'             => \$help,
            'version'          => \$version);

if ( defined($help) ) {
    print "DedupUMI version $versionno ($versiondate)\n";
    print "https://github.com/alxelerator/dedupUMI\n";
    print "By a.bossers\@uu.nl / alex.bossers\@wur.nl\n";

    print "\nUsage:\n";
    print "  deduplicate_UMI.pl \\\n";
    print "      --input-fastq1 <R1.fastq[.gz]> \\\n";
    print "      --input-fastq2 <R2.fastq[.gz]> \\\n";
    print "      --output-fastq1 <R1_out.fastq[.gz]> \\\n";
    print "      --output-fastq2 <R2_out.fastq[.gz]> \\\n";
    print "      [--input-fastq3 <UMI.fastq[.gz]>] \\\n";
    print "      [--output-fastq3 <UMI_out.fastq[.gz]>] \\\n";
    print "      [--output_counts <counts.tab>] \\\n";
    print "      [--noUMI] \\\n";
    print "      [--verbose] \n";

    print "\nRequired arguments:\n";
    print "  --input-fastq1 <file>     FASTQ read1 input file (plain or gz)\n";
    print "  --input-fastq2 <file>     FASTQ read2 input file (plain or gz)\n";
    print "  --output-fastq1 <file>    FASTQ read1 output file (plain or gz)\n";
    print "  --output-fastq2 <file>    FASTQ read2 output file (plain or gz)\n";

    print "\nOptional arguments:\n";
    print "  --input-fastq3 <file>     UMI FASTQ file (R3 system). If omitted, UMI is\n";
    print "                            extracted from the read header.\n";
    print "  --output-fastq3 <file>    Output FASTQ for R3 system UMIs\n";
    print "  --output_counts <file>    Append input/output read counts to tab file\n";
    print "  --noUMI                   Two-file mode without UMI. Deduplicate on R1+R2 only\n";
    print "  --verbose                 Print sequences to STDOUT (debug output)\n";
    print "  --help                    Show this help message\n";

    print "\nNotes:\n";
    print "  - UMI location (header or separate FASTQ file) is detected automatically unless --noUMI is used\n";
    print "  - If you have problems reading/writing gzipped files, ensure the 'zcat'\n";
    print "    command is available on your system.\n";
    print "  - On Unix systems gzip performance can often be improved by replacing\n";
    print "    gzip with 'pigz' (parallel gzip).\n";
    print "  - If using --noUMI we just filter exact duplicates. This most likely filters to harsh since\n";
    print "    exact duplicates can occurr in natural good quality datasets. Especially at high sequencing depths.\n";

    print "\nExamples:\n";
    print "  UMI in header:\n";
    print "    deduplicate_UMI.pl --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz \\\n";
    print "                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz\n";
    print "  R3 UMI system:\n";
    print "    deduplicate_UMI.pl --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --input-fastq3 s1_R3.fastq.gz \\\n";
    print "                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz --output-fastq3 s1_R3.dedup.fastq.gz\n";
    print "  No UMI:\n";
    print "    deduplicate_UMI.pl --input-fastq1 s1_R1.fastq.gz --input-fastq2 s1_R2.fastq.gz --noUMI \\\n";
    print "                       --output-fastq1 s1_R1.dedup.fastq.gz --output-fastq2 s1_R2.dedup.fastq.gz\n\n";
    exit 0;
}
    
if( defined($version) ) {
   print "DedupUMI version $versionno ($versiondate)\nhttps://github.com/alxelerator/dedupUMI\nBy a.bossers\@uu.nl / alex.bossers\@wur.nl\n";
   exit 0;
}

#error state
if( ! defined($input_fastq1) || ! defined($input_fastq2) || ! defined($output_fastq1 )|| ! defined($output_fastq2) ) {
    print STDERR "ERROR: One of the required arguments --input-fastq1, --input-fastq2, --input-fastq3 or --output-fastq1, --output-fastq2, --output-fastq3 is missing!\n" if !defined($help);
    print "DedupUMI version $versionno ($versiondate)\n";
    print "Use: deduplicate_UMI.pl --help for options\n";
    exit 1;
}

# UMI compatibility chaeck of options
if ( defined($noUMI) && defined($input_fastq3) ) {
    print STDERR "ERROR: --noUMI cannot be combined with --input-fastq3\n";
    exit 1;
}
if ( defined($noUMI) && defined($output_fastq3) ) {
    print STDERR "ERROR: --noUMI cannot be combined with --output-fastq3\n";
    exit 1;
}

# create tabular output file if not exists otherwise append the data to it
my $COUNTS;
if( defined($output_counts) ) {

    if (-e $output_counts) {
        print "dedupUMI: Appending read count data to existing file!\n";
        open ( $COUNTS, ">>$output_counts" ) || die "Output file error: $output_counts\n$!\n" ;
    } else {
        open ( $COUNTS, ">$output_counts" ) || die "Output file error: $output_counts\n$!\n" ;
        print $COUNTS "Sample\tReads_in\tReads_out\n";
    }
}


########################################################################################################
#
# Lets go!
#
########################################################################################################

# Check if we have UMI in the headers (new), in a separate file (R3-system), or no UMI at all
my $umiheader = 0;

if ( defined($input_fastq3) ) {
    # R3-system processing (R1 R2 R3) where UMI is in R3 file
    if ( ! defined($output_fastq3) ) {
        print STDERR "DeDupe ERROR: R3-system (R1 R2 R3) detected, but output file R3 is not defined! Define by using --output-fastq3\n" if !defined($help);
        exit 1;
    }
    print "DeDupe: UMI sequences processing R3-system (R1 R2 R3)\n";

} elsif ( defined($noUMI) ) {
    # noUMI processing (R1 R2 only)
    print "DeDupe: Processing noUMI mode (R1 R2 exact duplicate read-pairs)\n";

} else {
    # UMI in header
    print "DeDupe: UMI sequences processing UMI in the R1 R2 headers!\n";
    $umiheader = 1;
}



# Unified processing for R3-system and UMI-in-header system

my @infiles  = ($input_fastq1, $input_fastq2);
my @outfiles = ($output_fastq1, $output_fastq2);

if ( defined($input_fastq3) ) {
    push @infiles,  $input_fastq3;
    push @outfiles, $output_fastq3;
}

my $nfiles = scalar(@infiles);
my (@INFQ, @OUT);

# open INPUT files
for my $file (@infiles) {
    my $fh;
    if ( substr($file,-3) eq ".gz" || substr($file,-5) eq ".gzip" ) {
        open($fh, "zcat $file |") || die "Input file error: $file\n$!\n";
    } else {
        open($fh, $file) || die "Input file error: $file\n$!\n";
    }
    push @INFQ, $fh;
}

# open OUTPUT files
for my $file (@outfiles) {
    my $fh;
    if ( substr($file,-3) eq ".gz" || substr($file,-5) eq ".gzip" ) {
        open($fh, "|-", "gzip >$file") || die "Output file error: $file\n$!\n";
    } else {
        open($fh, ">$file") || die "Output file error: $file\n$!\n";
    }
    push @OUT, $fh;
}

# read lines and store best duplicates
my $reads_in = 0;
my $stored   = 0;
my %uniqseq  = ();

while ( my $fq1line1 = readline($INFQ[0]) ) {
    $reads_in++;

    # read header lines
    my @line1 = ($fq1line1);
    for my $i (1 .. $nfiles-1) {
        my $line = readline($INFQ[$i]);
        if (!defined $line) {
            die "ERROR: Paired FASTQ ended early / files out of sync (file ".($i+1)." missing header while R1 still has data)\n";
        }
        push @line1, $line;
    }

    if (!defined $line1[0] || $line1[0] eq '') {
        warn "WARNING: Empty header line in R1 (extra newline at EOF?)\n";
        last;
    }

    chomp @line1;

    # read sequence lines
    my @line2;
    for my $i (0 .. $nfiles-1) {
        my $line = readline($INFQ[$i]);
        die "ERROR: FASTQ ended early / missing sequence line in file ".($i+1)."\n" if !defined $line;
        chomp $line;
        push @line2, $line;
    }

    # get umi from header if needed
    my $umi;
    if ($umiheader) {
        if ( $line1[0] =~ /^@.+:[0-9]+:[0-9]+:([ATGCNatgcn]+) [0-9]+:.+/ ) {
            $umi = $1;
        } else {
            print STDERR "  ERROR: UMI not found in fastq header or fastq header format has changed!?\n    Header: '$line1[0]'\n";
            print STDERR "  Consider using the option --noUMI if your sequence library did not contain UMI sequences!\n";
            exit 1;
        }
    }

    # skip plus lines
    for my $i (0 .. $nfiles-1) {
        my $line = readline($INFQ[$i]);
        die "ERROR: FASTQ ended early / missing '+' line in file ".($i+1)."\n" if !defined $line;
    }

    # read quality lines
    my @line4;
    for my $i (0 .. $nfiles-1) {
        my $line = readline($INFQ[$i]);
        die "ERROR: FASTQ ended early / missing quality line in file ".($i+1)."\n" if !defined $line;
        chomp $line;
        push @line4, $line;
    }

    # store seq header and quality
    my %contents = ();
    $contents{"head"} = join("#alx#", @line1);
    $contents{"qual"} = join(" ", @line4);

    my $key;
    if ($umiheader) {
        $key = "$line2[0] $line2[1] $umi";
    } elsif ( defined($input_fastq3) ) {
        $key = "$line2[0] $line2[1] $line2[2]";
    } else {
        $key = "$line2[0] $line2[1]";
    }

    # Check if current sequence is already stored. If new seq quality is better replace existing.
    if ( exists( $uniqseq{$key} ) ) {

        my $total = 0;
        my $orgtotal = 0;
        my @orgquality = split(" ", $uniqseq{$key}{"qual"});

        if ( $umiheader || !defined($input_fastq3) ) {
            my @fq1line4num = unpack("C*", $line4[0]);
            my @fq2line4num = unpack("C*", $line4[1]);
            $total += $_ for(@fq1line4num);
            $total += $_ for(@fq2line4num);

            my @orgfq1line4num = unpack("C*", $orgquality[0]);
            my @orgfq2line4num = unpack("C*", $orgquality[1]);
            $orgtotal += $_ for(@orgfq1line4num);
            $orgtotal += $_ for(@orgfq2line4num);
        } else {
            my @fq1line4num = unpack("C*", $line4[0]);
            my @fq3line4num = unpack("C*", $line4[2]);
            $total += $_ for(@fq1line4num);
            $total += $_ for(@fq3line4num);

            my @orgfq1line4num = unpack("C*", $orgquality[0]);
            my @orgfq3line4num = unpack("C*", $orgquality[2]);
            $orgtotal += $_ for(@orgfq1line4num);
            $orgtotal += $_ for(@orgfq3line4num);
        }

        if ($total <= $orgtotal) {
            next;
        }
    }

    # store the new sequence in the hash
    $uniqseq{$key} = \%contents;
    $stored++;
}

print "DeDupe: Read    : $reads_in sequences\n";

##################################################
# write all unique sequences out                 #
##################################################

my $seqcount_w = 0;
foreach my $seqkey (keys %uniqseq) {

    print "  $seqkey\n" if $verbose;

    my @seq     = split(" ", $seqkey);
    my @headers = split("#alx#", $uniqseq{$seqkey}{"head"});
    my @quality = split(" ", $uniqseq{$seqkey}{"qual"});

    for my $i (0 .. $nfiles-1) {
        print {$OUT[$i]} $headers[$i]."\n";
        print {$OUT[$i]} $seq[$i]."\n";
        print {$OUT[$i]} "+\n";
        print {$OUT[$i]} $quality[$i]."\n";
    }

    $seqcount_w++;
}

print "Dedupe: Written : $seqcount_w sequences\n";
print "DeDupe: Finished ".(strftime "%m/%d/%Y %H:%M:%S", localtime)."\n\n";

if( defined($output_counts) ) {
    print $COUNTS "$input_fastq1\t$reads_in\t$seqcount_w\n";
    close ($COUNTS);
}

# end and close
for my $fh (@INFQ) {
    close($fh);
}
for my $fh (@OUT) {
    close($fh);
}


#end script
