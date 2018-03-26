#!/bin/bash
pop_size=5000
generations=10

n_loci=$5
s=$1
selection="T"
rep=$3
selection_type=$2
block_size=0.1
epistasis_selection=$4
label=$6
perl run_simulation.pl \
	--n_generations $generations \
	--block_size $block_size \
	--n_samples $pop_size \
	--selection $selection \
	--s $s \
	--selection_type $selection_type \
	--n_selected_loci $n_loci \
	--epistasis_selection $epistasis_selection \
	2>log/log_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${generations}.txt |\
gzip > simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${generations}.txt.gz

