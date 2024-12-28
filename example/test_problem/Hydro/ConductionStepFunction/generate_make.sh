# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=mac --double=true --model=HYDRO --eos=GAMMA --hdf5=True --openmp=False --flu_scheme=MHM --conduction=True "$@"
