# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true --gsl=true \
                       --model=HYDRO --particle=true --gravity=true "$@"
