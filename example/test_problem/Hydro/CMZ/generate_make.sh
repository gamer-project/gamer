# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --par_attribute_flt=1 "$@"
