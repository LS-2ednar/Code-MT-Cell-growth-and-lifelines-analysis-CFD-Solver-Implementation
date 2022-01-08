#!/bin/bash
#
#SBATCH --job-name=EXAMPLE_GROWTHSIM2
#
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --partition=parallel
#SBATCH --qos=parallel
#SBATCH --mail-type=ALL				# BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-user=YOUR@MAIL.ADDRESS		# Email adress

module load USS/2020 
module load slurm
module load gcc/7.3.0
module load openmpi/3.1.3
module load openfoam-org/9.0


mpi run -np 8 particleFoam 2>&1 | tee particleFoam.log
