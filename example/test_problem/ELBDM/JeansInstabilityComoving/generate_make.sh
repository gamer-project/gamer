# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --hdf5=true --gpu=true --fftw=FFTW3 --double=true \
                       --model=ELBDM --gravity=true --comoving=true "$@"
