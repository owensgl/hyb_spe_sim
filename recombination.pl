#!/bin/perl
use warnings;
use strict;
use Math::Random qw(random_poisson);

sub recombination {
  #Get lengths of chromosomes
  my ($data1, $data2) = @_;
  my @data1 = split(/\t/,$data1);
  my @data2 = split(/\t/,$data2);
  my $length = length($data1[3]);
  my $n_recomb = Math::Random::random_poisson(1, 1);
  #    print STDERR "\nThe number of recombinations is $n_recomb";
  if ($n_recomb == 0){ 
    #If there are no recombination events, just return one random chromosome.
    my $strand = int(rand(2));
    my $result;
    if ($strand == 1){
      $result = $data1[3];
    }else{
      $result = $data2[3];
    }
    return($result);
  }
  my %reco; #The location of all recombination events;
  for (1..$n_recomb){
    my $pos_in_series = int(rand($length+1));
    $reco{$pos_in_series}++;
  }
  
  my @info1 = split(//,$data1[3]);
  my @info2 = split(//,$data2[3]);
  my $rand_start = int(rand(2));
  my $print_string;
  foreach my $i (0..$#info1){
    if ($rand_start){
      $print_string .= $info1[$i];
    }else{
      $print_string .= $info2[$i];
    }
    if ($reco{$i}){
      if ($rand_start){
        $rand_start = 0;
      }else{
        $rand_start = 1;
      }
    }
  }
  return($print_string);
}
1;
