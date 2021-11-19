#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3GB
#SBATCH --time=00:20:00
#SBATCH --partition=conti
#SBATCH --account=dconti_624
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=jagoodri@usc.edu 

module purge
module load gcc/11.2.0
module load openblas/0.3.18
module load jags
module load openmpi
module load pmix
module load r/4.1.2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun --mpi=pmix_v2 -n $SLURM_NTASKS Rscript --vanilla 0_pbdmpi.R