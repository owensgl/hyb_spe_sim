#!/bin/perl
use warnings;
use strict;
use lib '/home/owens/bin/hyb_spe_sim/';

require "recombination.pl";
require "heterozygosity.pl";
my $previous_line;
while (<STDIN>){
    chomp;
    my $line = $_;
    if ($line =~ /^#/){
        next;
    }
    if ($line =~ /^sample_ID/){
        next;
    }
    unless($previous_line){
        $previous_line = $line;
    }else{
  #      my $newchrom_1 = recombination($line, $previous_line);
#	my $newchrom_2 = recombination($line, $previous_line);
#	my $newline_1 = "1\t0\tHa1\t$newchrom_1";
 #       my $newline_2 = "1\t1\tHa1\t$newchrom_2";
#	my $het = heterozygosity($newline_1, $newline_2);
#	print STDERR "Chrom1 = $newchrom_1\n";
#	print STDERR "Chrom2 = $newchrom_2\n";
 #       print STDERR "\nThe heterozygosity is $het\n";
	my $chrom_1 = "X\t0\tHa2\t0;P1;83.614";
	my $chrom_2 = "X\t1\tHa2\t0;P1;51.1554961304622,51.1554961304622;P2;83.614";
	my $new_chrom = recombination($chrom_1, $chrom_2);
	print STDERR "$new_chrom\n";
	exit;
    }
}
