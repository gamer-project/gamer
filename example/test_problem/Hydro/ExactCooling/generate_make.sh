# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --exact_cooling=true --gpu=true "$@"
