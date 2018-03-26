#!/bin/perl
use warnings;
use strict;
use Math::Round;
#This script makes F1s
#This version of the hybrid speciation keeps track of each block instead of junctions.
my $n_size = $ARGV[0];

my $smallest_unit = $ARGV[1];
my %c;
my @chrs;
my $chrfile = "/home/gowens/bin/ref/HanXRQr1.0-20151230.bp_to_cM.280x801.chrsizes.txt";
open CHR, $chrfile;
while(<CHR>){
  chomp;
  my @a = split(/\t/,$_);
  if ($a[0] eq "Total"){next;}
  my $rounded = nearest($smallest_unit, $a[1]);
  $c{$a[0]} = $rounded;
  push(@chrs,$a[0]);
}
close CHR;
print "#\tsamplesize\t$n_size";
foreach my $i (0..($n_size-1)){
    my $n = $i+1;
    foreach my $chr (@chrs){
	my $rand = int(rand(1));
        foreach my $j (0..1){
            my $type;
	    if ($rand == 1){
              if ($j == 0){
                  $type = "0";
              }else{
                  $type = "1";
              }
	    }else{
              if ($j == 0){
                  $type = "1";
              }else{
                  $type = "0";
              }
	    }
            print "\n$n\t$j\t$chr\t";
            for (my $i=0; $i <= $c{$chr}; $i+=$smallest_unit){
              print "$type";
            }
        }
    }
}


