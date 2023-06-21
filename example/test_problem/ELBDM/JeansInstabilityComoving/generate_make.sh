# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 --double=true \
                       --model=ELBDM --gravity=true --comoving=true
