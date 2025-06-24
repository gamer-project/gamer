# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=ELBDM --double=True --hdf5=true --fftw=FFTW3 "$@"
