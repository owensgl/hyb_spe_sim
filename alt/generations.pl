#!/bin/perl
use warnings;
use strict;
use Math::Random qw(random_poisson);
use Math::NumberCruncher;
use Math::Round;
use lib '/home/owens/bin/hyb_spe_sim/alt/';
require "recombination.pl";
require "heterozygosity.pl";
#ARGUMENTS
my $generations = $ARGV[0];
my $smallest_unit = $ARGV[1];
my $loci_selection = $ARGV[2];
my $selection_strength = $ARGV[3];
my $n_selected_loci = $ARGV[4];
my $selection_type = $ARGV[5];
my $sample_size = $ARGV[6];

my $chrfile = "/home/owens/ref/HanXRQr1.0-20151230.bp_to_cM.280x801.chrsizes.txt";
my %chrlengths;
my @chrs;
my $print_every = 10;
open CHR, $chrfile;
while(<CHR>){
    chomp;
    my @a = split(/\t/,$_);
    if ($a[0] eq "Total"){next;}
    $chrlengths{$a[0]} = $a[1];
    push(@chrs,$a[0]);
}

my %previous_generation;
my $n_chromosomes = scalar @chrs;
#Pick the locations of the selected loci.
#print "#s:$selection_strength";
#print "\n#n_s:$n_selected_loci";
my $total_s = 1 -(1-$selection_strength)**($n_selected_loci * $n_chromosomes);
print STDERR "\n#total_s = $total_s";
my $direction = 0;
my %loci;
my %loci_direction;
#This picks loci evenly distributed across the chromosome.
foreach my $chr (@chrs){
  foreach my $i (1..$n_selected_loci){
    my $n = (($chrlengths{$chr}/($n_selected_loci+1)) * $i);
    my $rounded_n = nearest($smallest_unit, $n);
    $loci{$chr}{$rounded_n}++;
    $loci_direction{$chr}{$rounded_n} = $direction;
    if($direction == 1){
      $direction = 0;
    }else{
      $direction = 1;
    }
  }
}

while (<STDIN>){
  chomp;
  my $line = $_;
  if ($line =~ /^#/){
    next;
  }
  my @a = split(/\t/,$line);
  my $sample_ID = $a[0];
  my $strand = $a[1];
  my $chr = $a[2];
  my $info = $a[3];
  $previous_generation{$sample_ID}{$chr}{$strand} = $info;
}
foreach my $gen (1..$generations){
  my $selected_loci_het;
  my $selected_loci_homo;
  my $current_sample_size = 0;
  my %current_generation;
  my @current_heterozygosity;
  until ($current_sample_size == $sample_size){
    my $n_offspring = Math::Random::random_poisson(1, 2); #Randomly chooses the number of offspring.
    if ($n_offspring == 0){ next;} #Skips cases where they produce no offspring.
    my $parent1 = int(rand($sample_size))+1;
    my $parent2 = int(rand($sample_size))+1;
    if ($parent1 eq $parent2){next;}
    foreach my $i (1..$n_offspring){
      my $fitness = 1;
      my %tmp_gamete_storage;
      my $current_sample_ID = $current_sample_size+1;
      #For each chromosome
      foreach my $chr (@chrs){

        #Generate both gametes
        my $input_chr_1_0 = "X\t0\t$chr\t$previous_generation{$parent1}{$chr}{0}";
        my $input_chr_1_1 = "X\t1\t$chr\t$previous_generation{$parent1}{$chr}{1}";
        my $input_chr_2_0 = "X\t0\t$chr\t$previous_generation{$parent2}{$chr}{0}";
        my $input_chr_2_1 = "X\t1\t$chr\t$previous_generation{$parent2}{$chr}{1}";
        my $gamete1 = recombination($input_chr_1_0,$input_chr_1_1);
        my $gamete2 = recombination($input_chr_2_0,$input_chr_2_1);
        #measure heterozygosity
        if ($gen % $print_every == 0){
          my $gamete1_full = "$current_sample_ID\t0\t$chr\t$gamete1";
          my $gamete2_full = "$current_sample_ID\t1\t$chr\t$gamete2";
          my $het = heterozygosity($gamete1_full,$gamete2_full);
          #                       print STDERR "My heterozygosity measured is $het\n";
          push(@current_heterozygosity,$het);
        }
        if ($loci_selection == 1){
	  my @current_selected_loci = keys %{$loci{$chr}};
          #SELECTION FUNCTION
	  my @bases1 = split(//,$gamete1);
          my @bases2 = split(//,$gamete2);
	  foreach my $selected_locus (@current_selected_loci){
            my $state;
            if ($bases1[$selected_locus] ne $bases2[$selected_locus]){
	      $state = "H";
	      $selected_loci_het++;
            }elsif ($bases1[$selected_locus] eq "0"){
	      $state = "0";
	      $selected_loci_homo++;
	    }else{
	      $state = "1";
	      $selected_loci_homo++;
	    }
	    if ($selection_type eq "U"){
	      if ($state eq "H"){
		$fitness = $fitness * (1 - $selection_strength);
	      }
	    }elsif ($selection_type eq "D"){
	      if ($state eq "H"){
		$fitness = $fitness * (1 - ($selection_strength/2));
	      }elsif ($state eq $loci_direction{$chr}{$selected_locus}){
		#No reduction in fitness
	      }else{
		$fitness = $fitness * (1 - $selection_strength);
	      }
	    }else{die "Selection type must be U or D"};
          }
	}
        $tmp_gamete_storage{$chr}{0} = $gamete1;
        $tmp_gamete_storage{$chr}{1} = $gamete2;
      }
      my $luck = rand(1);
      if ($luck > $fitness){
        goto NEXTOFFSPRING;
      }
      foreach my $chr (@chrs){
        #If it does survive, load the chromosomes to the next generation.
        $current_generation{$current_sample_ID}{$chr}{0} = $tmp_gamete_storage{$chr}{0};
        $current_generation{$current_sample_ID}{$chr}{1} = $tmp_gamete_storage{$chr}{1};
      }
      #Add to sample size tracker and check to make sure it doesn't go over size (because of multiple children in one pair)
      $current_sample_size++;
      if ($current_sample_size == $sample_size){
        goto FULLSIZE;
      }
      NEXTOFFSPRING:
    }

  }
  FULLSIZE:
  if ($gen % $print_every == 0){
    #calculate stats for heterozygosity.
    my $mean_het = Math::NumberCruncher::Mean(\@current_heterozygosity);
    my $std_het = Math::NumberCruncher::StandardDeviation(\@current_heterozygosity);
    print STDERR "Generation $gen: $mean_het +- $std_het";
    if ($loci_selection == 1){
      unless($selected_loci_het){$selected_loci_het = 0;}
      my $selected_loci_heterozygosity = $selected_loci_het / ($selected_loci_homo + $selected_loci_het);
      print STDERR "\tSelected het = $selected_loci_heterozygosity";
    }
    print STDERR "\n";
  }
  #refresh the variables for the next generation
  %previous_generation = %current_generation;
}
foreach my $n (1..$sample_size){
  foreach my $chr (@chrs){
    print "\n$n\t0\t$chr\t$previous_generation{$n}{$chr}{0}";
    print "\n$n\t1\t$chr\t$previous_generation{$n}{$chr}{1}";
  }
}

