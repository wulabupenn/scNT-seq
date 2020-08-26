#!/usr/bin/perl

##############################################################
#  Author     : Peng Hu
#  Email      ：penghu41@gmail.com
#  Last Edited: 08/26/2020
##############################################################

use strict;
use warnings;
use Getopt::Long;
my %opt;


GetOptions (\%opt,"bg:s","in:s","help");

my $help=<<USAGE;
This script could correct TFEA mRNA by background mRNA

Usage: perl $0 -bg bg.tsv -in in.tsv

USAGE

if ($opt{help} or keys %opt < 1) {
	print "$help\n";
	exit();
}

my $in1 = $opt{bg};
my $in2 = $opt{in};

my %hash;

open (IN1,$in1) or die $!;
#open (IN2,$in2) or die $!;
#open (OUT,">${in}_corrected.tsv") or die $!;
#my $line= <IN>;
#print $line;
while (my $line=<IN1>) {
	chomp $line;
	my @array = split "\t",$line;
        my $loc = join "_",($array[1],$array[2],$array[6]);
        $hash{$loc} = "exist";
}
close IN1;


open (IN2,$in2) or die $!;
open (OUT,">${in2}_corrected.tsv") or die $!;
open (OUT2,">${in2}_discard.tsv") or die $!;
while (my $line=<IN2>) {
       chomp $line;
       my @array = split "\t",$line;
       my $loc = join "_",($array[1],$array[2],$array[6]);
       if (! $hash{$loc}) {
          print OUT $line,"\n";
       } else {
          print OUT2 $line,"\n";
       }
}
close OUT;
close OUT2;
close IN2;




