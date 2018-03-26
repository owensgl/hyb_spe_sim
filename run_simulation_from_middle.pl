#!/bin/perl
use Getopt::Long;
#This script takes all the variables for the hybrid speciation simulation and runs the scripts

my $n_generations   = 2; #Number of generations to run the simulation
my $n_samples = 5000; #Number of individuals in population
my $selection_type = "D"; #U for underdominant or D for directional
my $selection = "T"; #T for true, F for false. 
my $s = "0.1"; #This is the selection_strength.
my $n_selected_loci = "2"; #Number of selected loci
my $block_size = "0.1"; #The size of genetic boxes in cm size
my $epistasis_selection = "T"; #Whether there is epistatic selection
my $input_file; #The simulated individuals to start the simulation with
GetOptions (
	"n_generations=i" => \$n_generations,    
	"n_samples=i"   => \$n_samples,
	"selection_type=s"  => \$selection_type,
	"selection=s" => \$selection,
	"s=f" => \$s,
	"n_selected_loci=i" => \$n_selected_loci,
	"block_size=f" => \$block_size,
	"epistasis_selection=s" => \$epistasis_selection,
	"input=s" => \$input_file
		);
print STDERR "n_generations = $n_generations\n";
print STDERR "n_samples = $n_samples\n";
print STDERR "selection_type = $selection_type\n";
print STDERR "selection = $selection\n";
print STDERR "selection_strength = $s\n";
print STDERR "n_selected_loci = $n_selected_loci\n";
print STDERR "block_size = $block_size\n";
print STDERR "epistasis_selection = $epistasis_selection\n";
print STDERR "input = $input_file\n"; #Should be gzipped
unless ($input_file){
	kill "Need an input file";
}
my $command = "zcat $input_file | perl ./generations.pl $n_generations $block_size $selection $s $n_selected_loci $selection_type $n_samples $epistasis_selection";
system($command); 
