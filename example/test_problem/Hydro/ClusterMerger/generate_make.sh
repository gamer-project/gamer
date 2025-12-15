# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=HYDRO --particle=true --gravity=true --passive=2 --par_attribute_int=1 "$@"


# cool-core descruction
# two clusters
# ${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM --flux=HLLC --passive=1 \
#                        --particle=true --store_par_acc=true --par_attribute_int=1 \
#                        --gravity=true --unsplit_gravity=true --fftw=FFTW3 \
#                        --mpi=true --gpu=true \
#                        --hdf5=true --gsl=true \
#                        --exact_cooling=true \
#                        "$@"

# single cluster
# ${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM --flux=HLLC --passive=2 \
#                        --particle=true --store_par_acc=true --par_attribute_int=1 \
#                        --gravity=true --unsplit_gravity=true --fftw=FFTW3 \
#                        --mpi=true --gpu=true \
#                        --hdf5=true --gsl=true \
#                        --exact_cooling=true \
#                        "$@"
