#!/bin/bash
#
#SBATCH --job-name=EXAMPLE_TRACKING_1
#
#SBATCH --nodes=8
#SBATCH --ntasks=256
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


blockMesh 2>&1 | tee blockMesh.log
surfaceFeatures 2>&1 | tee surfaceFeatures.log
snappyHexMesh -overwirte 2>&1 | tee snappyHexMesh.log 
decomposePar 2>&1 | tee decomposPar.log
mpirun -np 256 particleFoam -parallel 2>&1 | tee particleFoam.log
reconstructPar 2>&1 | tee reconsturctPar.log
rm -r processor*
