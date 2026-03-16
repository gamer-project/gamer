# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --double=true --hdf5=true --conduction=true --gravity=true --mhd=true --fftw=FFTW3 --flu_scheme=MHM_RP "$@"
