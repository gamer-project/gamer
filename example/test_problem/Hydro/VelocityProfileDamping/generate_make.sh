# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --hdf5=true --mhd=true --viscosity=true --double=True --flu_scheme=MHM_RP "$@"
