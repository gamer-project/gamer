# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=HYDRO --particle=true --tracer=true --gravity=true
