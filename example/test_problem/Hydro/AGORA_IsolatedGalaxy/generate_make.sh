# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=HYDRO --particle=true --gravity=true --flu_scheme=MHM --flux=HLLC \
                       --passive=1 --par_attribute=1 --dual=DE_ENPY --star_formation=true --grackle=true
