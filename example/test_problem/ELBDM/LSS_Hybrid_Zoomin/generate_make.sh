# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=ELBDM --elbdm_scheme=HYBRID --wave_scheme=GRAMFE --gramfe_scheme=MATMUL \
                       --gravity=true --comoving=true --gsl=true --spectral_interpolation=true "$@"
