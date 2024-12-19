# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --gpu=true --hdf5=true --fftw=FFTW3 --model=HYDRO --gravity=true "$@"
