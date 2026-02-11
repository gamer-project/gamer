# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --flu_scheme=MHM --flux=HLLD --mhd=true --eos=USER \
                       --barotropic=true --gravity=true --particle=true --par_attribute_int=1 \
                       --star_formation=true --feedback=true --double=true --double_par=true \
                       --hdf5=true --gsl=true --fftw=FFTW3 "$@"