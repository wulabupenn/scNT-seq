#!/usr/bin/perl

##############################################################
#  Author     : Peng Hu
#  Email      ：penghu41@gmail.com
#  Last Edited: 08/26/2020
##############################################################

use strict;
use warnings;

my $in = shift @ARGV;
my %hash;

my $bn=`basename $in`;

my $file_bn=$1 if ($bn =~/(.*)\.bam$/);

open (IN,"samtools view -h $in |") or die "Unable to read from BAM file $in: $!\n";
while (<IN>) {
   chomp;
   if (/GE:Z:(\S+)/) {
      my $gene = $1;
      my @array = split "\t";
	  if ($array[4] > 10) {
	  if (/XC:Z:(\w+)/) {
	      my $cell_barcode = $1;
		  if (/XM:Z:(\w+)/) {
	          my $umi = $1;
			  my $index = join "\t",($gene,$cell_barcode,$umi);
			  if (! $hash{$index}) {
			      $hash{$index} = 1;
			     } else {
				  $hash{$index} += 1;
				 }
			  }
		}
      }
   }
}
close IN;

#my %hash2;
open (OUT1,">${file_bn}_gene_cell_UMI_read.txt") or die $!;

foreach my $index(keys %hash) {
    print OUT1 "$index\t$hash{$index}\n";
#	my @array = split "\t",$index;
#	pop @array; # remove the last element of UMI
#    my $index2 = join "\t",@array;
#	if (! $hash2{$index2}) {
#	    $hash2{$index2} = 1;
#	} else {
#	   $hash2{$index2} += 1;
#	}
}
close OUT1;

#open (OUT2,">${file_bn}_gene_cell_UMI.txt") or die $!;
#foreach my $index2(keys %hash2) {
#   print OUT2 "$index2\t$hash2{$index2}\n";
#}
