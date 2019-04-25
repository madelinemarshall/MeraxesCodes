#!/bin/bash
#SBATCH --job-name=low_edd
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --time=12:00:00		#after this time the job will be killed
#SBATCH --mem-per-cpu=4000	#memory in MB per core
#SBATCH -D /home/mmarshal/data_dragons/paper2_low_edd
#SBATCH -o meraxes.out		#output file
#SBATCH -e meraxes.err		#error file, will go to outfile if omitted


srun meraxes input.par
