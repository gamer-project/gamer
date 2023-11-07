# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 \
                       --model=ELBDM --elbdm_scheme=ELBDM_WAVE \
                       --gravity=true --comoving=false