#!/bin/perl
use strict;
use warnings;
#This script takes the simulation of hybrid species with blocks and measures the average ancestry over windows for the purpose of measuring correlation coefficients between origins.
my $window_size = 5;
my @stored_strand;
print "sample\tchr\twindow\tancestry";
while(<STDIN>){
  chomp;
  if ($. == 1){next;
  }
  my @a = split(/\t/,$_);
  my $sample = $a[0];
  my $strand = $a[1];
  my $chr = $a[2];
  my $data = $a[3];
  my $value = 0;
  my $sites = 0;
  if ($sample > 1){exit;} #Only do one sample;
  if ($strand == 0){
    @stored_strand = split(//,$a[3]);
  }else{
    my @current_strand = split(//,$a[3]);
    foreach my $i (0..$#current_strand){
      if (($i % 5 == 0) and ($i > 0)){
        my $window = int($i / 5);
        my $ancestry = $value/ $sites;
        print "\n$sample\t$chr\t$window\t$ancestry";
        $value = 0;
        $sites = 0;
      }
      $sites+=2;
      $value+=$stored_strand[$i];
      $value+=$current_strand[$i];
    }
    my $i = $#current_strand;
    my $window = int($i / 5);
    my $ancestry = $value/ $sites;
    print "\n$sample\t$chr\t$window\t$ancestry";
  }
}
