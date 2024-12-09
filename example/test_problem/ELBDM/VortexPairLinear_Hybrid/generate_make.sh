# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --model=ELBDM --elbdm_scheme=ELBDM_HYBRID --hdf5=true --mpi=true
