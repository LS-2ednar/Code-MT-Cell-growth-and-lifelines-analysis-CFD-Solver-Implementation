#!/bin/bash
#
#SBATCH --job-name=EXAMPLE_TRACKING3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --partition=tiny
#SBATCH --qos=tiny
#SBATCH --mail-type=ALL				# BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-user=YOUR@MAIL.ADDRESS		# Email adress

module load USS/2020
module load gcc/7.3.0
module load python/3.6.5-pe5.26
module load py-pip/19.3-py3.6-pe5.26
module load openblas/0.3.5-haswell-ep
module load py-numpy/1.17.3-openblas-py3.6-pe5.26
module load gcc/7.3.0 python/3.6.5-pe5.26 openblas/0.3.5-haswell-ep py-pandas/0.25.1-openblas-py3.6-pe5.26 py-matplotlib/3.1.3-openblas-py3.6-pe5.26
module load slurm
module load gcc/7.3.0
module load openmpi/3.1.3
module load openfoam-org/9.0

python particlelifelines.py 2>&1 | tee particlelifelines.log
