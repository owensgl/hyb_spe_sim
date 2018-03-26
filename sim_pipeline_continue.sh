#! /bin/bash
#Run using cat simulation_stats.txt | parallel -j 10 --colsep ' ' bash ./full_sim_pipeline.sh {1} {2} {3} {4} {5} {6}
pop_size=5000
max_generation=3000
generations=10
n_loci=$5
s=$1
selection="T"
rep=$3
selection_type=$2
block_size=0.1
epistasis_selection=$4
label=$6
find simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci} -size  0 -print0 |xargs -0 rm
existing=$(ls simulations | grep "ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}")
previous_start_gen=$(echo $existing | cut -d "_" -f 8 | cut -f 1 -d '.')
new_end_gen=$((previous_start_gen + generations))
until [ $previous_start_gen == $max_generation ]; do

	if [[ $(find simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${previous_start_gen}.txt.gz -type f -size +5000c 2>/dev/null) ]]
	then
		perl run_simulation_from_middle.pl \
		--n_generations $generations \
		--block_size $block_size \
		--n_samples $pop_size \
		--selection $selection \
		--s $s \
		--selection_type $selection_type \
		--n_selected_loci $n_loci \
		--epistasis_selection $epistasis_selection \
		--input simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${previous_start_gen}.txt.gz \
		2>log/log_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${new_end_gen}.txt |\
gzip > simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${new_end_gen}.txt.gz
		if [[ $(find simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${new_end_gen}.txt.gz -type f -size +5000c 2>/dev/null) ]]
		then
			rm simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_${previous_start_gen}.txt.gz
		fi
	else
		echo "${previous_start_gen} step didn't work. Stopping"
	fi
	previous_start_gen=$new_end_gen
	new_end_gen=$((previous_start_gen + generations))
done
