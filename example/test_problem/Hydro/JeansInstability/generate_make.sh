# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --gpu=true --fftw=FFTW3 --double=true \
                       --model=HYDRO --gravity=true --eos=GAMMA "$@"
