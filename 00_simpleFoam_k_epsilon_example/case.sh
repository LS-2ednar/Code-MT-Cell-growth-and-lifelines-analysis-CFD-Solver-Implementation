#!/bin/bash
#
#SBATCH --job-name=100RPM_8C
#
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --partition=single
#SBATCH --qos=single
#SBATCH --mail-type=ALL								# BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-user=schauluk@students.zhaw.ch		# Email adress

module load USS/2020 
module load slurm
module load gcc/7.3.0
module load openmpi/3.1.3
module load openfoam-org/9.0
 
blockMesh 2>&1 | tee blockMesh.log
surfaceFeatures | tee surfaceFeatures.log
snappyHexMesh -overwrite 2>&1  | tee snappyHexMesh.log
checkMesh 2>&1 | tee checkMesh.log
decomposePar 2>&1 | tee decomposPar.log
mpirun -np 8 simpleFoam -parallel 2>&1 | tee simpleFoam.log
reconstructPar 2>&1 | tee reconsturctPar.log
rm -r processor*
