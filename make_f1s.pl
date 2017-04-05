#!/bin/perl
use strict;
use warnings;

#This script creates a number of F1 individuals equal to the population size
#Requires chromosome file with list
#Chr_number\tlength_of_chr_in_cm
my $n_size = $ARGV[0]; #Need to put in the number of individuals.
my $chrfile = "/home/owens/ref/HanXRQr1.0-20151230.bp_to_cM.280x801.chrsizes.txt";
my %chrlengths;
my @chrs;
open CHR, $chrfile;
while(<CHR>){
    chomp;
    my @a = split(/\t/,$_);
    if ($a[0] eq "Total"){next;}
    $chrlengths{$a[0]} = $a[1];
    push(@chrs,$a[0]);
}
#For simulated 1 chromosome
#$chrlengths{"TMP"} = 100;
#push(@chrs, "TMP");
#Sample format
#sample_ID\tstrand\tchr\tstart\tend\ttype
print "#\tsamplesize\t$n_size\n";
print "sample_ID\tstrand\tchr\tinfo";
foreach my $i (0..($n_size-1)){
    my $n = $i+1;
    foreach my $chr (@chrs){
        foreach my $j (0..1){
            my $type;
            if ($j == 0){
                $type = "P1";
            }else{
                $type = "P2";
            }
            print "\n$n\t$j\t$chr\t0;$type;$chrlengths{$chr}";
        }
    }
}

