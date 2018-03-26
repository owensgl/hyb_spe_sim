#!/bin/bash
#SBATCH --account=rrg-rieseber-ac
#SBATCH --time=24:00:00
#SBATCH --job-name=simulations
#SBATCH --array=0-5%140
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --exclusive=user
#SBATCH --output=/home/gowens/bin/logs/%A.%a.hybsim_log
cd /scratch/gowens/hyb_sim
startline=$(( SLURM_ARRAY_TASK_ID * 32 + 1 )) 
endline=$(( startline + 31 ))
input_contig=$(cat simulation_parameters.txt  | sed -n ${startline},${endline}p)

cat simulation_parameters.txt  | sed -n ${startline},${endline}p| parallel -j 32 --colsep ' ' bash ./sim_pipeline_continue.sh {1} {2} {3} {4} {5} {6}
