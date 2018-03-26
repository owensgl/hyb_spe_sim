#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=01:00:00
#SBATCH --job-name=simulations
#SBATCH --array=1%140
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive=user
#SBATCH --output=/home/gowens/bin/logs/%A.%a.hybsim_log
input_contig=$(cat simulation_parameters.txt  | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $input_contig | parallel -j 1 --colsep ' ' bash ./full_sim_pipeline.sh {1} {2} {3} {4} {5} {6} && echo "Simulation complete" 
