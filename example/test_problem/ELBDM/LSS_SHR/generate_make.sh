# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=ELBDM --elbdm_scheme=ELBDM_HYBRID --wave_scheme=WAVE_GRAMFE --gramfe_scheme=GRAMFE_MATMUL \
                       --gravity=true --comoving=true --gsl=true --spectral_interpolation=true
