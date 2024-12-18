# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --double=true --model=HYDRO --mhd=true --eos=GAMMA "$@"
