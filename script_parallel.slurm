#!/bin/bash
#SBATCH --job-name=plots
#SBATCH --ntasks=14
#SBATCH --time=6:00:00		#after this time the job will be killed
#SBATCH --mem-per-cpu=4000	#memory in MB per core
#SBATCH -D /home/mmarshal/simulation_codes/Paper1Plots
#SBATCH -o plots.out		#output file
#SBATCH -e plots.err		#error file, will go to outfile if omitted


srun --multi-prog program_list.conf
