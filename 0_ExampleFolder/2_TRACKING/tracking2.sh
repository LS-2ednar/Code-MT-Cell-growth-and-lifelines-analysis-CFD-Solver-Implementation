#!/bin/bash
#
#SBATCH --job-name=EXAMPLE_TRACKING_2
#
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --partition=tiny
#SBATCH --qos=tiny
#SBATCH --mail-type=ALL				# BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-user=YOUR@MAIL.ADDRESS		# Email adress

module load USS/2020 
module load slurm
module load gcc/7.3.0
module load openmpi/3.1.3
module load openfoam-org/9.0


particleTracks
