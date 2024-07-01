# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --hdf5=true --gpu=true --mpi=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true "$@"
