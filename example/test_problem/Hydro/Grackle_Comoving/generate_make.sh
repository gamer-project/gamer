# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --hdf5=true --model=HYDRO --comoving=true --gravity=true --fftw=FFTW3 \
                       --grackle=true --double=true --passive=1 "$@"
