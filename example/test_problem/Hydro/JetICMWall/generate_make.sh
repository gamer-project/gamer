# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=pleiades --mpi=true --hdf5=true --gpu=true --model=HYDRO \
	               --srhd=true --eos=TAUBMATHEWS --flux=HLLC --flu_scheme=MHM \
	               --gpu_arch=VOLTA --passive=4
