#!/bin/bash

#SBATCH --job-name=YOUR_JOB_NAME          ## job name
#SBATCH --mail-type=ALL                   ## mail notification type
#SBATCH --mail-user=YOUR_EMAIL            ## mail address
#SBATCH --account=YOUR_ACCOUNT            ## account (https://www.twcc.ai/user/dashboard)
#SBATCH --partition=gp1d                  ## queue name (gtest/gp1d/gp2d/gp4d)
#SBATCH --nodes=2                         ## number of nodes
#SBATCH --gres=gpu:8                      ## number of GPUs per node
#SBATCH --ntasks-per-node=8               ## number of tasks per node
#SBATCH --cpus-per-task=4                 ## number of CPUs per task
#SBATCH --time=02:00:00                   ## time limit (hh:mm:ss)
#SBATCH --output=log-%j                   ## stdout output file

#Load module
module purge
module load compiler/gnu/4.8.5 nvidia/cuda/10.0 openmpi/3.1.4

srun ./gamer 1>>log 2>&1
