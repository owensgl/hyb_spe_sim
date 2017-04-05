#!/bin/perl
use warnings;
use strict;
use Math::Random qw(random_poisson);
use Math::Round;
#This script takes the genotype data an individual (both strands, one chromosome), and creates a recombined chromosome.
#Uses a poisson distribution of recombination events (mean of 1 per chr)


sub recombination {
    #Get lengths of chromosomes
    my %c;
    my $chrfile = "/home/owens/ref/HanXRQr1.0-20151230.bp_to_cM.280x801.chrsizes.txt";
    open CHR, $chrfile;
    while(<CHR>){
      chomp;
      my @a = split(/\t/,$_);
      if ($a[0] eq "Total"){next;}
      $c{$a[0]} = $a[1];
    }
    close CHR;
#    $c{"Ha1"} = "81.75";
#    $c{"Ha2"} = "83.614";
#    $c{"Ha3"} = "75.704";
#    $c{"Ha4"} = "95.992";
#    $c{"Ha5"} = "88.057";
#    $c{"Ha6"} = "59.68";
#    $c{"Ha7"} = "54.031";
#    $c{"Ha8"} = "68.464";
#    $c{"Ha9"} = "91.984";
#    $c{"Ha10"} = "87.89";
#    $c{"Ha11"} = "84.686";
#    $c{"Ha12"} = "70.222";
#    $c{"Ha13"} = "70.556";
#    $c{"Ha14"} = "76.295";
#    $c{"Ha15"} = "75.343";
#    $c{"Ha16"} = "99.127";
#    $c{"Ha17"} = "100.772";
#   $c{"TMP"} = "100";
    #Load up data passed to subroutine (2 strands of 1 chromosome)
    my ($data1, $data2) = @_;
    my @data1 = split(/\t/,$data1);
    my @data2 = split(/\t/,$data2);
    my $chr = $data1[2];
    my $n_recomb = Math::Random::random_poisson(1, 1);
#    print STDERR "\nThe number of recombinations is $n_recomb";
    if ($n_recomb == 0){ #If there are no recombination events, just return one random chromosome.
	my $strand = int(rand(2))+1;
	my $result;
	if ($strand == 1){
	    $result = $data1[3];
	}else{
	    $result = $data2[3];
	}
	return($result);
    }


    my @reco; #The location of all recombination events;
    for (1..$n_recomb){
        my $loc = rand($c{$chr});
	$loc = nearest(0.001, $loc); #This rounds the value to 0.0001
        push(@reco,$loc);
    }
    my @recombination_locations = sort {$a <=> $b } @reco; #Sorts them in order


#    print STDERR "The recombination locations are @recombination_locations\n";
    #pull out the actual chromosome information, split out chunks.
    my @info1 = split(/,/,$data1[3]);
    my %info;
    foreach my $i (0..$#info1){
        my @tmp = split(/;/,$info1[$i]);
        $info{1}{$i}{"start"} = $tmp[0];
        $info{1}{$i}{"ID"} = $tmp[1];
        $info{1}{$i}{"end"} = $tmp[2];
    }
    my @info2 = split(/,/,$data2[3]);
    foreach my $i (0..$#info2){
        my @tmp = split(/;/,$info2[$i]);
        $info{2}{$i}{"start"} = $tmp[0];
        $info{2}{$i}{"ID"} = $tmp[1];
        $info{2}{$i}{"end"} = $tmp[2];
    }
    my %n_breaks;
    $n_breaks{1} = $#info1;
    $n_breaks{2} = $#info2;
    my $strand = int(rand(2))+1;
    my $result;
    my $current_pos = 0;
    RETRY:
    until($current_pos == $c{$chr}){
        foreach my $n (0..$n_breaks{$strand}){
            if (($current_pos >= $info{$strand}{$n}{"start"}) and
            ($current_pos < $info{$strand}{$n}{"end"})){ #Checks if the current position is within this window.
		if ($recombination_locations[0]){
	                if ($recombination_locations[0] < $info{$strand}{$n}{"end"}){ #The recombination occurs within this chunk of parentage
        	            $result .= qq($current_pos;$info{$strand}{$n}{"ID"};$recombination_locations[0],);

                	    $current_pos = shift @recombination_locations; #Remove the recombination location because it's done, and put it in the current position.
                    	    if ($strand eq "1"){ #Switch to other strand
                        	$strand = 2;
                    	    }else{
                        	$strand = 1;
                    	    }
                    	    goto RETRY;
			}else{
	             	    $result .= qq($current_pos;$info{$strand}{$n}{"ID"};$info{$strand}{$n}{"end"},);
                            $current_pos = $info{$strand}{$n}{"end"};
                            goto RETRY;
			}
			
                }else{ #The recombination event is not here.
                    $result .= qq($current_pos;$info{$strand}{$n}{"ID"};$info{$strand}{$n}{"end"},);
                    $current_pos = $info{$strand}{$n}{"end"};
                    goto RETRY;
                }
            } 
		#
        }
    }
    chop($result); #Remove the hanging comma
#    print STDERR "The unmerged result is $result\n";
    #This is to merge together multiple blocks of the same ancestry
    my $merged_result;
    my @final_parts = split(/,/,$result);
    my $part_length = $#final_parts;
    my %final_hash;
    foreach my $i (0..$part_length){
	my @tmp = split(/;/,$final_parts[$i]);
	$final_hash{$i}{"start"} = $tmp[0];
	$final_hash{$i}{"ID"} = $tmp[1];
	$final_hash{$i}{"end"} = $tmp[2];
    }
    my $current = 0;
    until ($current == $part_length+1){
	no warnings 'uninitialized'; #Be careful of this suppression of warnings. It is specifically for $final_hash{$current+1}{"ID"} in the last case.
	$merged_result .= qq($final_hash{$current}{"start"};$final_hash{$current}{"ID"});
	until($final_hash{$current}{"ID"} ne $final_hash{$current+1}{"ID"}){
	    $current++;
        }
	$merged_result .= qq(;$final_hash{$current}{"end"},);
	$current++;
    }
    chop($merged_result);
    return($merged_result);
}
1; #This makes sure it returns a true value on test
