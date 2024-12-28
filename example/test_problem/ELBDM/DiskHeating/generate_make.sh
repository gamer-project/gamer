# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=ELBDM --gravity=true --particle=true --store_par_acc=true \
                       --gsl=true --max_patch=20000000
