#!/bin/perl
use warnings;
use strict;
use Scalar::Util qw(looks_like_number);

#This script takes in two chromosomes and measures interspecific heterozygosity.
#It looks at each pair of blocks between the strands that overlaps, and also has the same type, it measures the total length.
sub heterozygosity{
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
#    $c{"Ha1"} = "78.515";
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
#	$c{"TMP"} = "100";

    my ($data1, $data2) = @_;
    my @data1 = split(/\t/,$data1);
    my @data2 = split(/\t/,$data2);
    my @info1 = split(/,/,$data1[3]);
    my @info2 = split(/,/,$data2[3]);
    my $chr = $data1[2];
    my %info;

    foreach my $i (0..$#info1){
        my @tmp = split(/;/,$info1[$i]);
        $info{1}{$i}{"start"} = $tmp[0];
        $info{1}{$i}{"ID"} = $tmp[1];
        $info{1}{$i}{"end"} = $tmp[2];
    }
    foreach my $i (0..$#info2){
        my @tmp = split(/;/,$info2[$i]);
        $info{2}{$i}{"start"} = $tmp[0];
        $info{2}{$i}{"ID"} = $tmp[1];
        $info{2}{$i}{"end"} = $tmp[2];
    }
    my %n_breaks;
    $n_breaks{1} = $#info1;
    $n_breaks{2} = $#info2;
    my $hom_length = 0;

    foreach my $i (0..$n_breaks{1}){ #For each strand one block
        foreach my $j (0..$n_breaks{2}){ #compared to each strand 2 block
            if ($info{1}{$i}{"ID"} eq $info{2}{$j}{"ID"}){ #If they're the same
		if (($info{1}{$i}{"start"} > $info{2}{$j}{"end"}) or ($info{2}{$j}{"start"} > $info{1}{$i}{"end"})){
			#Checks if the blocks overlap at all.
			next;
		}
		#This is the start value that is the most forward.
                my $forward_start = ($info{1}{$i}{"start"}, $info{2}{$j}{"start"})[$info{1}{$i}{"start"} < $info{2}{$j}{"start"}];
#		print STDERR qq(Between $info{1}{$i}{"start"} and $info{2}{$j}{"start"}, $forward_start is bigger\n);
                my $backward_end = ($info{1}{$i}{"end"}, $info{2}{$j}{"end"})[$info{1}{$i}{"end"} > $info{2}{$j}{"end"}];
                my $overlap = $backward_end - $forward_start;
                $hom_length+=$overlap;
            }
        }
    }
    my $percent_het = ($c{$chr} - $hom_length) / $c{$chr};
#    print STDERR "The homozygosity is $hom_length, with full length $c{$chr}\n";
    return($percent_het);

}
1;
