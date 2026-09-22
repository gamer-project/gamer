#!/bin/bash

#PBS -N yt_analysis
#PBS -M your_email@example.com
#PBS -m abe
#PBS -q workq
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe


if [ "$PBS_ENVIRONMENT" == "PBS_BATCH" ]; then
   cd $PBS_O_WORKDIR
fi

# ===== Parameters Setting =====
START=0
END=75
DELTA=1


# ===== Run =====
python3 dust_density.py
python3 gas_temp.py -s $START -e $END -d $DELTA
