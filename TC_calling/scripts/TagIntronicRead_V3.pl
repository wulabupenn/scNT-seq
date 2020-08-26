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


GetOptions (\%opt,"gtf:s","bam:s","help");


my $help=<<USAGE;
This script could annotate the intronic mapped reads by adding "GE:Z:gene" and "GS:Z:strand" in the end of the line, because the information of GE and GS is required by the Dropseq pipeline.

Usage: perl $0 -gtf gtf.file -bam file.bam

You could check the results by:
samtools view file.TagIntronic.bam | grep "XF:Z:INTRONIC" | grep "GE:Z" | more

USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help\n";
	exit();
}

#my $samtools=`which samtools`;
#print $samtools,"\n";

my $gtf = $opt{gtf};
my $file = $opt{bam};
#my $outfile = $opt{out};

my $bn=`basename $opt{bam}`;

my $file_bn=$1 if ($bn =~/(.*)\.bam$/);

my (%GeneList,%locS,%locE);

open (GTF, $gtf) or die "Unable to open file $gtf: $! \n";
while (<GTF>) {
	if (/gene/) {
		my @array = split "\t";  ## chr: $array[0], $strand: $array[6]
		if ($array[2] eq 'gene') {
			my $gene;
			if (/gene_name\s\"(\S+?)\"/) {
				$gene = $1;
			}
            $locS{$gene} = $array[3];
            $locE{$gene} = $array[4];

			my $unit=$array[0]."_".$array[6];
            if (! $GeneList{$unit}) {
				$GeneList{$unit} = $gene;
            } else {
                $GeneList{$unit} = $GeneList{$unit}.";".$gene;
			}
		}
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
open (OUT,"| samtools view -bSh 2>/dev/null - >${dir1}${file_bn}.TagIntronic.bam") or die $!;
open (OUT2,"| samtools view -bSh 2>/dev/null - >${dir1}${file_bn}.TagIntronicOnly.bam") or die $!;

my $line_count=0;
my $intronic_line_count=0;
my $line;
my $uncount_intronic=0;
while ($line=<DETERMINE>){
	if ($line !~ /^@/) {
		++$line_count;
	    warn "Add Gene information into the intronic Reads! Processed lines: $line_count\n" if ($line_count%1000000==0);
	}
    if ($line=~/XF:Z:INTRONIC/) {
		chomp $line;
        my ($id,$strand,$chrom,$start,$cigar) = (split("\t",$line))[0,1,2,3,5];
        my $strand2;
		## determine the strand information
        if ($strand & 16) {
		    $strand2 = "-";
        } else {
		    $strand2 = "+";
	    }
      
	    ## detemine the start and end location 
		
		### The part of below script about how to extract information from CIGAR was inspired by Felix Krueger's Bismark pipeline (bismark_methylation_extractor)
		### https://github.com/FelixKrueger/Bismark/blob/ee6952a0aaf078fda9dffa9b926716529cc0eaf6/bismark_methylation_extractor
		
        my @len = split (/\D+/,$cigar); # storing the length per operation
	    my @ops = split (/\d+/,$cigar); # storing the operation
	    shift @ops;		# remove the empty first element
	    die "CIGAR string contained a non-matching number of lengths and operations\n" unless (scalar @len == scalar @ops);

        ## obtain the end location, considering the M, D, and N in the cigar
         my $adjust_count=0;
	     foreach my $index (0..$#ops) {
	       if ($ops[$index] =~ /[MDN]/) {
		      $adjust_count += $len[$index];
	       }
	     }

         #print "$adjust_count\n";
         my $end = $start + $adjust_count - 1;
         
		 my $unit = $chrom."_".$strand2;

         my @genelis = split (";",$GeneList{$unit});
        
		my $index;

        #print sort(@genelis),"\n";

        ## reterive the infromation from gtf file
		#print $unit,"\n";
        #print $locS{"Rgs20"},"\n";
		#print $locE{"Rgs20"},"\n";
		#print $start,"\n";
		#print $end,"\n";

		foreach my $gene(@genelis) {
			if ( $locS{$gene} <= $start && $end <= $locE{$gene} ) {
				print OUT "$line\tGE:Z:$gene\tGS:Z:$strand2\n";
                                print OUT2 "$line\tGE:Z:$gene\tGS:Z:$strand2\n";
				++$intronic_line_count;
				$index=0;
				last;
			} else {
				$index=1;
			}
		}

		if ($index) {
			++$uncount_intronic;
			#print "$id\t$start\t$adjust_count\n";
			print OUT "$line\n";
		}
  
    } else {
		print OUT "$line";
                print OUT2 "$line" if ($line=~/^@/);
	}
}
close DETERMINE;
close OUT;
close OUT2;

print "################# SUMMARY #################### \n";

print "$uncount_intronic intronic reads cannot be determined due to out the boundaries of a single-gene or mapped to the another strand \n";
print "$intronic_line_count intronic reads were tagged with the gene name and strand information \n";
print "In total, $line_count reads were processed \n";

print "############################################## \n";

