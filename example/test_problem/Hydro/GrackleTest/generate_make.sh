# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --gsl=true --gpu=true --model=HYDRO --grackle=true --passive=14 "$@"
