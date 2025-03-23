# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --mpi=true --hdf5=true --fftw=FFTW3 --gravity=true --gpu=true \
                       --particle=true --feedback=true --par_attribute_flt=2 --passive=1 "$@"
