# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true --gpu_arch=TURING \
                       --model=ELBDM --gravity=true --particle=true --par_attribute=1 --store_par_acc=true \
                       --gsl=true --max_patch=20000000
