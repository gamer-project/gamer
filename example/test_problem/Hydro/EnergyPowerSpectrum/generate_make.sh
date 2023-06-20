# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --fftw=FFTW3 --model=HYDRO --eos=GAMMA --mpi=true
