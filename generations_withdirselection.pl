#!/bin/perl
use warnings;
use strict;
use Math::Random qw(random_poisson);
use Math::NumberCruncher;
use lib '/home/owens/bin/hyb_spe_sim/';
require "recombination.pl";
require "heterozygosity.pl";

#This script takes a set of samples and then randomly mates them together to produce a new set of samples (including recombination). Later, I will add selection on specific genomic locations
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
#$chrlengths{"TMP"} = 100;
#push(@chrs, "TMP");
my $generations = $ARGV[0];
my $n_selected_loci = $ARGV[1];
my $selection_strength = $ARGV[2];
my $sample_size;
my %previous_generation;
my $print_every = 10; #print every x generations
my $loci_selection = "TRUE"; #Turns on selection at specific loci
#my $n_selected_loci = "1"; #The number of selected loci per chromosome
#my $selection_strength = "0.1"; #The probability of death with heterozygote at selected loci. Multaplicative with multiple loci.
my $n_chromosomes = 17;
my @loci;
#Pick the locations of the selected loci. Set up for fake 100cm chromosomes
print "#s:$selection_strength";
print "\n#n_s:$n_selected_loci";
my $total_s = (1-$selection_strength)**($n_selected_loci * $n_chromosomes);
print "\n#total_s = $total_s";
#This picks loci every 10 cm.
foreach my $i (1..$n_selected_loci){
	my $n = ($i+1) * 20;
	push(@loci, $n);
	print "\n#l:$n";
}

while (<STDIN>){
    chomp;
    my $line = $_;
    if ($line =~ /^#/){
        my @a = split(/\t/,$line);
        $sample_size = $a[2];
        next;
    }
    if ($line =~ /^sample_ID/){
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
#    print "Working on generation $gen\n";
    my $selected_loci_homo = 0;
    my $selected_loci_het = 0;
    my %current_generation;
    my $current_sample_size = 0;
    my @current_heterozygosity;
    until ($current_sample_size == $sample_size){
        my $n_offspring = Math::Random::random_poisson(1, 2); #Randomly chooses the number of offspring.
        if ($n_offspring == 0){ next;} #Skips cases where they produce no offspring.
        my $parent1 = int(rand($sample_size))+1;
        my $parent2 = int(rand($sample_size))+1;
        if ($parent1 eq $parent2){next;}
#	print STDERR "The parents chosen are $parent1 and $parent2\n";
        #Insert selection function here if needed.

        foreach my $i (1..$n_offspring){
	    my $fitness = 1;
	    my %tmp_gamete_storage;
	    my $current_sample_ID = $current_sample_size+1;
            foreach my $chr (@chrs){ #For each chromosome
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
#			print STDERR "My heterozygosity measured is $het\n";
                	push(@current_heterozygosity,$het);
		}
		#Selection function against specific region
		if ($loci_selection eq "TRUE"){
			foreach my $loc (@loci){
				my %type;
				my @blocks1 = split(/,/,$gamete1);
				foreach my $block (@blocks1){
					my @ends = split(/;/,$block);
					if (($ends[0] <= $loc) and ($ends[2] >= $loc)){
						$type{1} = $ends[1];
					}
				}
                                my @blocks2 = split(/,/,$gamete2);
                                foreach my $block (@blocks2){
                                        my @ends = split(/;/,$block);
                                        if (($ends[0] <= $loc) and ($ends[2] >= $loc)){
                                                $type{2} = $ends[1];
                                        }
                                }
				if ($type{1} ne $type{2}){
					$fitness = $fitness * (1 - $selection_strength); #lowers fitness if heterozygote
					$selected_loci_het++;
				}else{
					$selected_loci_homo++;
				}
			}
		}
		$tmp_gamete_storage{$chr}{0} = $gamete1;
		$tmp_gamete_storage{$chr}{1} = $gamete2;
            }
	    my $luck = rand(1);
	    if ($luck > $fitness){
		goto NEXTOFFSPRING;
	    }
	    foreach my $chr (@chrs){ #If it does survive, load the chromosomes to the next generation.
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
	if ($loci_selection eq "TRUE"){
	  my $selected_loci_heterozygosity = $selected_loci_het / ($selected_loci_homo + $selected_loci_het);
	  print STDERR "\tSelected het = $selected_loci_heterozygosity";
	}
	print STDERR "\n";
    }
    #refresh the variables for the next generation
    %previous_generation = %current_generation;
}
#print "sample_ID\tstrand\tchr\tinfo";
foreach my $n (1..$sample_size){
	foreach my $chr (@chrs){
		print "\n$n\t0\t$chr\t$previous_generation{$n}{$chr}{0}";
		print "\n$n\t1\t$chr\t$previous_generation{$n}{$chr}{1}";
	}
}
