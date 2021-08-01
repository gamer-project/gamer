#!/bin/bash
#PBS -N GAMER
#PBS -M PUT_YOUR_EMAIL_HERE
#PBS -l nodes=1:ppn=16:xk
#PBS -l walltime=8:00:00
#PBS -m bea
#PBS -q normal
##PBS -q high
##PBS -q debug
##PBS -e $PBS_JOBID.err
##PBS -o $PBS_JOBID.out

cd $PBS_O_WORKDIR

. /opt/modules/default/init/bash
module load cray-hdf5/1.8.16
module load cudatoolkit/7.5.18-1.0502.10743.2.1
module load craype-hugepages2M

aprun -n 1 -N 1 -d 16 --cc=none ./gamer 1>>log 2>&1 </dev/null

