# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --srhd=true --double=true --eos=TAUBMATHEWS --flux=HLLC --flu_scheme=MHM --slope=PLM "$@"
