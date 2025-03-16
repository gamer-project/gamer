# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true --model=HYDRO \
                       --particle=true --gravity=true --flu_scheme=MHM --flux=HLLC --passive=1 \
                       --par_attribute_flt=1 --dual=ENPY --star_formation=true --grackle=true "$@"
