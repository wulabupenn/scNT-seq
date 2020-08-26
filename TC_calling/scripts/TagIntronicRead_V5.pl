#!/usr/bin/env perl

##############################################################
#  Author     : Peng Hu
#  Email      ：penghu41@gmail.com
#  Last Edited: 08/26/2020
##############################################################

use warnings;
use strict;
$|++;
use Getopt::Long;
#use Cwd;
use Cwd 'abs_path';
use Carp;
use FindBin qw($Bin);
use lib "$Bin/../lib";

my %opt;


GetOptions (\%opt,"read:s","bam:s","help");


my $help=<<USAGE;
This script could annotate TC convertion by changing "GE:Z:gene" to GE:Z:gene--T(unconverted)  or GE:Z:gene--C(converted) in the end of the line.

Usage: perl $0 -read read.tsv -bam file.bam

You could check the results by:
samtools view file.bam | grep "GE:Z" | more

USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help\n";
	exit();
}

#my $samtools=`which samtools`;
#print $samtools,"\n";

my $gtf = $opt{read};
my $file = $opt{bam};
#my $outfile = $opt{out};

my $bn=`basename $opt{bam}`;

my $file_bn=$1 if ($bn =~/(.*)\.bam$/);

my (%readlist);

open (GTF, $gtf) or die "Unable to open file $gtf: $! \n";
while (<GTF>) {
   chomp;
   my @array = split "\t";
   if (! $readlist{$array[0]}) {
       $readlist{$array[0]} = "exist";
   }
}
close GTF;


my $dir1 = abs_path($file);
my @paths = split "/",$dir1;
pop @paths;
$dir1 = join "/",@paths;
$dir1 .= "/";
#print "my dir1 is $dir1 \n";


open (DETERMINE,"samtools view -h $file |") or die "Unable to read from BAM file $file: $!\n";
#open (DETERMINE,"$file") or die "Unable to read from file $file: $!\n";
#open (OUT,">${file_bn}.sam") or die $!;
open (OUT,"| samtools view -bSh 2>/dev/null - >${file_bn}.TagTC.corrected.bam") or die $!;

my $line_count=0;
my $line;

while ($line=<DETERMINE>){
	if ($line !~ /^@/) {
		++$line_count;
	    warn "Add TC convertion to gene! Processed lines: $line_count\n" if ($line_count%1000000==0);
	}
    if ($line=~/(GE:Z:\S+)/) {
		my $geneinfo = $1;
		chomp $line;
        my $read = (split("\t",$line))[0];
        if ($readlist{$read}) {
		   $line =~ s/$geneinfo/$geneinfo--C/;
		} else {
		   $line =~ s/$geneinfo/$geneinfo--T/;
		}
		print OUT "$line\n";
    } else {
		print OUT "$line";
	}
}
close DETERMINE;
close OUT;

