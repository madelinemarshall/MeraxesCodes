#!/bin/bash
#SBATCH --job-name=history
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=12:00:00		#after this time the job will be killed
#SBATCH --mem-per-cpu=4000	#memory in MB per core
#SBATCH -D /home/mmarshal/simulation_codes
#SBATCH -o history.out		#output file
#SBATCH -e history.err		#error file, will go to outfile if omitted


srun python extract_gal_hist.py 2 10 2 9 9.1
