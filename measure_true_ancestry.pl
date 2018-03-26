#!/bin/perl
use strict;
use warnings;

#This script takes the simulated ancestry from a hybrid species simulation and calculates the true ancestry proportion for 10MB windows of the genome.
my $ref_file = "/home/gowens/ref/HanXRQr1.0-20151230.bp_to_cM.280x801.windows.txt";

my $samplecounter;
my $HA17_counter= 0;
my %sampleinfo;
while(<STDIN>){
    chomp;
    next if /^\s*$/;
    my @a = split(/\t/,$_);
    my $strand = $a[1];
    my $chr = $a[2];
    $chr =~ s/chr//;
    $chr = sprintf("%02d",$chr);
    $chr =~ s/^/HanXRQChr/;
    my $info = $a[3];
    $sampleinfo{$chr}{$strand} = $info;
    #This is to track when you're moving from one sample into another.
    if ($chr eq "HanXRQChr17"){
        $HA17_counter++;
    }
    if ($HA17_counter == 2){
        $samplecounter++;
        $HA17_counter = 0;
    }
}

print "chr\tstart\tend\ttrue_admixture";
open REF, $ref_file;
while(<REF>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos_start = $a[1];
  my $pos_end = $a[2];
  my $cm_start = $a[3];
  my $cm_end = $a[4];
  my $cm_size = $a[5];
  my $value = 0;
  foreach my $strand (0..1){
    my $info = $sampleinfo{$chr}{$strand};
    my @fields = split(/,/,$info);
    my $current_type;
    foreach my $field(@fields){
      my @tmp = split(/;/,$field);
      my $start = $tmp[0];
      my $type = $tmp[1];
      my $end = $tmp[2];
      #Skip windows that don't overlap at all;
      if (($start > $cm_end) or ($end < $cm_start)){next;}
      #Find the start and end of the overlap region;
      my $overlap_start;
      my $overlap_end;
      if ($start > $cm_start){
        $overlap_start = $start;
      }else{
        $overlap_start = $cm_start;
      }
      if ($end > $cm_end){
        $overlap_end = $cm_end;
      }else{
        $overlap_end = $end;
      }
      my $overlap_dist = $overlap_end - $overlap_start;
      if ($type eq "P2"){
        $value += $overlap_dist;
      }
    }
    
  }
  my $final_value = $value / (2 * $cm_size);
  print "\n$chr\t$pos_start\t$pos_end\t$final_value";
}
