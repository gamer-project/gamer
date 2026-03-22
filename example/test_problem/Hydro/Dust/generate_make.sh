# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --gpu=true --model=HYDRO \
                       --passive=2  --double=true \
                       --grackle=true "$@"

