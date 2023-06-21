# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --mpi=true \
                       --model=HYDRO --mhd=true
