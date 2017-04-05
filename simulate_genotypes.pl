#!/bin/perl
use warnings;
use strict;
my $genofile = $ARGV[0]; #SNP table in ctab format
my $parentfile = $ARGV[1];

my %pop;
open PARENT, $parentfile;
while(<PARENT>){
    chomp;
    my @a = split(/\t/,$_);
    $pop{$a[0]} = $a[1];
}
close PARENT;

my %sampleinfo;
my $samplecounter = 1;
my $HA17_counter = 0;
while(<STDIN>){
    chomp;
    next if /^\s*$/;
    my @a = split(/\t/,$_);
    my $sample_ID = $samplecounter;
    my $strand = $a[1];
    my $chr = $a[2];
    $chr =~ s/chr//;
    $chr = sprintf("%02d",$chr);
    $chr =~ s/^/HanXRQChr/; 
    my $info = $a[3];
    $sampleinfo{$sample_ID}{$chr}{$strand} = $info;
    #This is to track when you're moving from one sample into another.
    if ($chr eq "HanXRQChr17"){
        $HA17_counter++;
    }
    if ($HA17_counter == 2){
        $samplecounter++;
        $HA17_counter = 0;
    }
}
$samplecounter--;
#Start the final output file.
print "chr\tpos";
my %samplepop;
my %samplename;
my %good_number_hash;
open GENO, $genofile;
while (<GENO>){
    chomp;
	my @a = split(/\t/,$_);
	if ($. == 1){ #Load in sample names associated with column numbers, as well as population
		foreach my $i(3..$#a){
			if ($pop{$a[$i]}){
				$samplepop{$i} = $pop{$a[$i]};
				$samplename{$i} = $a[$i];
		                if (($pop{$a[$i]} eq "P1") or ($pop{$a[$i]} eq "P2")){
                		    $good_number_hash{$i}++;
                		}
			}
		}
		foreach my $i (sort keys %good_number_hash){
			print "\t$a[$i]";
		}
		foreach my $i (1..$samplecounter){
    			print "\t$i";
		}
    }else{
        my $chr = $a[0];
        my $pos = $a[1];
        my $cm = $a[2];
	if ($cm eq "NA"){ next;}
	my $P1count = 0;
	my $P2count = 0;
        foreach my $i (keys %good_number_hash){ #Load up parental alleles
            if ($a[$i] ne "NN"){
                if ($samplepop{$i} eq "P1"){
                    $P1count++;
                }elsif($samplepop{$i} eq "P2"){
                	$P2count++;
                }
            }
        }
		unless(($P1count >=5) and ($P2count >= 5)){
            next;
        }
        my $min_count;
		if ($P1count < $P2count){
                	$min_count = $P1count;
        	}else{
                	$min_count = $P2count;
        	}
		my @P1alleles;
		my @P2alleles;
		my $P1count2 = 0;
		my $P2count2 = 0;
        	print "\n$chr\t$pos";
		foreach my $i (sort keys %good_number_hash){
			print "\t$a[$i]";
		}
			
		foreach my $i(keys %good_number_hash){
			if ($samplepop{$i}){
				unless (($a[$i] eq "NN")or($a[$i] eq "XX")){
					my @bases = split(//, $a[$i]);
					if (($samplepop{$i} eq "P1") and ($P1count2 < $min_count)){
                        push(@P1alleles,$bases[0]);
                        push(@P1alleles,$bases[1]);
						$P1count2++;
					}elsif(($samplepop{$i} eq "P2") and ($P2count2 < $min_count)){
                        push(@P2alleles,$bases[0]);
                        push(@P2alleles,$bases[1]);
                        $P2count2++;
					}
				}
			}
		}
        #determine the  parentage at this location and draw a random allele
        foreach my $sample_ID (1..$samplecounter){
            my $genotype;
            foreach my $strand (0..1){
                my $info = $sampleinfo{$sample_ID}{$chr}{$strand};
                my @fields = split(/,/,$info);
                my $current_type;
                foreach my $field(@fields){
                    my @tmp = split(/;/,$field);
                    my $start = $tmp[0];
                    my $type = $tmp[1];
                    my $end = $tmp[2];
                    if (($cm >= $start) and ($cm <= $end)){
                        $current_type = $type;
                        goto NEXTSTRAND;
                    }
                }
                NEXTSTRAND:
		if ($current_type){
	                if ($current_type eq "P1"){
        	            my $tmp_base = $P1alleles[rand @P1alleles];
                	    $genotype .= $tmp_base;
	                }elsif ($current_type eq "P2"){
        	            my $tmp_base = $P2alleles[rand @P2alleles];
               	     	    $genotype .= $tmp_base;
			}
                }else{
			$genotype .= "N";
		}
            }
            print "\t$genotype";
        }
    }
}
