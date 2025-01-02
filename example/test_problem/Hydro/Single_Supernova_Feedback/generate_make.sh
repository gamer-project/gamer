# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=false --model=HYDRO \
                       --particle=true --gravity=true --flu_scheme=MHM --flux=HLLC \
		       --feedback=true \
                       --passive=0 --par_attribute=0 --dual=DE_ENPY --star_formation=true --grackle=false "$@"
