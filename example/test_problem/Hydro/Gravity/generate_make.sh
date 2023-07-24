# This script should run in the same directory as configure.py

PYTHON=python

${PYTHON} configure.py --machine=eureka_intel --gpu=true --fftw=FFTW3 --model=HYDRO --gravity=true
