# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --gpu=true --fftw=FFTW3 --double=true \
                       --model=HYDRO --gravity=true --eos=GAMMA
