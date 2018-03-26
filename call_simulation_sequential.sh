#!/bin/bash
steps=100
scenarios=180
steps_per_run=10
for i in `seq $scenarios`; do
   sample_jobs=()
   input_parameters=$(cat simulation_parameters.txt  | sed -n ${i}p)
   srun_conf="--time=01:00:00 --nodes=1 --ntasks=1 --mem=5G --account=rrg-rieseber-ac --output=log/%A.$i.log"

   n_loci=$(echo $input_parameters | cut -d " " -f 5)
   s=$(echo $input_parameters | cut -d " " -f 1)
   selection="T" 
   rep=$(echo $input_parameters | cut -d " " -f 3)
   selection_type=$(echo $input_parameters | cut -d " " -f 2)
   block_size=0.1
   epistasis_selection=$(echo $input_parameters | cut -d " " -f 4)
   label=$(echo $input_parameters | cut -d " " -f 6)

   existing=$(ls simulations | grep "ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}")
   #Run the first step of the simulation if no previous runs exist
   if [[ ! $(find simulations/ancestry_${label}_${s}_${selection_type}_${rep}_${epistasis_selection}_${n_loci}_* -type f -size +510c 2>/dev/null) ]]; then
     tmp=$(sbatch ${srun_conf} --job-name "$i.1.sim" ./sim_pipeline_start.sh $input_parameters )
     sample_jobs+=(${tmp##* })
     echo "Starting job $i.1.sim, ${tmp##* }"
     sleep 1
   fi
   starting_gen=$(echo $existing | cut -d "_" -f 8 | cut -f 1 -d '.')
   if [ -z "$starting_gen" ]; then
      starting_gen=$steps_per_run
   fi
   echo "Subsequent simulations are starting at $starting_gen"
   #Run more steps in the simulation
   for j in `seq 1 $steps`; do
       n=$((j - 1))
       start_gen=$((n * $steps_per_run + starting_gen))
       tmp=""
       if [ -z ${sample_jobs[$n]} ]; then
	  tmp=$(sbatch ${srun_conf} --job-name "$i.$start_gen.sim"  ./sim_pipeline_continue.sh $input_parameters $start_gen )
	  echo "Starting job $i.$start_gen.sim, ${tmp##* }. No dependency"
	  sample_jobs+=(${tmp##* }) #Add it twice to the list to make sure that the array has enough numbers to keep up.
       else
          tmp=$(sbatch ${srun_conf} --job-name "$i.$start_gen.sim"   --dependency=afterok:${sample_jobs[$n]}  ./sim_pipeline_continue.sh $input_parameters $start_gen )
	  echo "Starting job $i.$start_gen.sim, ${tmp##* }. After ${sample_jobs[$n]}"
       fi
       sample_jobs+=(${tmp##* })
       sleep 0.5s
   done
done

