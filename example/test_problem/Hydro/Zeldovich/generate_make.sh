# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --gsl=true --double=true --model=HYDRO --gravity=true \
                       --comoving=true --dual=ENPY --particle=true "$@"
