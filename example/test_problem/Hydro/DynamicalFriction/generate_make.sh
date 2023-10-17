# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --model HYDRO --gravity True --particle True --max_patch 1000000 --hdf5 True --gsl True --fftw FFTW3 --gpu True --mpi True
