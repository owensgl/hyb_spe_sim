#!/bin/perl
use warnings;
use strict;

#This script takes in two chromosomes and measures interspecific heterozygosity.

sub heterozygosity{
  my ($data1, $data2) = @_;
  my @data1 = split(/\t/,$data1);
  my @data2 = split(/\t/,$data2);
  my @info1 = split(//,$data1[3]);
  my @info2 = split(//,$data2[3]);
  my $chr = $data1[2];
  my $heterozygosity = 0;
  foreach my $i (0..$#info1){
    if ($info1[$i] != $info2[$i]){
      $heterozygosity++;
    }
  }
  my $percent_heterozygosity = $heterozygosity / ($#info1+1);
  return($percent_heterozygosity)
}
1;
