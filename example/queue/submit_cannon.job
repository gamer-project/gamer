#!/bin/bash
#SBATCH -n 2
#SBATCH -c 4
#SBATCH -t 0-03:10          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p itc_gpu          # Partition to submit to
#SBATCH -o gamer.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e gamer.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mem-per-cpu=4000
#SBATCH --gres=gpu:2
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=4

# load modules
module load gcc/12.2.0-fasrc01
module load mpich
module load hdf5
module load cuda

export CUDA_VISIBLE_DEVICES=0,1,2,3
export MPI_TYPE_MAX=655360
export MPI_TYPE_DEPTH=32
export MPI_MSGS_MAX=10485760
export MPI_BUFS_PER_PROC=256
export MPI_BUFS_PER_HOST=512
export MPI_USE_CUDA=1
export MPI_DSM_DISTRIBUTE=0
export KMP_AFFINITY=disabled
export HDF5_DISABLE_VERSION_CHECK=1

# run code

srun -n $SLURM_NTASKS --cpus-per-task=4 --mpi=pmix ./gamer
