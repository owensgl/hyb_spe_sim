#!/bin/perl
use Getopt::Long;
#This script takes all the variables for the hybrid speciation simulation and runs the scripts

my $n_generations   = 100;
my $n_samples = 100;
my $selection_type = "U"; #U for underdominant or D for directional
my $selection = "0";
my $s = "0";
my $n_selected_loci = "2";
my $block_size = "0.01";
GetOptions (
	"n_generations=i" => \$n_generations,    
	"n_samples=i"   => \$n_samples,
	"selection_type=s"  => \$selection_type,
	"selection" => \$selection,
	"s=f" => \$s,
	"n_selected_loci=i" => \$n_selected_loci,
	"block_size=f" => \$block_size
		);
print STDERR "n_generations = $n_generations\n";
print STDERR "n_samples = $n_samples\n";
print STDERR "selection_type = $selection_type\n";
print STDERR "selection = $selection\n";
print STDERR "selection_strength = $s\n";
print STDERR "n_selected_loci = $n_selected_loci\n";
print STDERR "block_size = $block_size\n";

my $command = "perl ./make_f1s.pl $n_samples $block_size | perl ./generations.pl $n_generations $block_size $selection $s $n_selected_loci $selection_type $n_samples";
system($command); 
