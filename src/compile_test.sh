set -e
PYTHON=python3

if [[ $HOSTNAME == "eureka00" ]]; then echo "This script should not run on login node!"; exit; fi
if [[ $HOSTNAME == "spock00"  ]]; then echo "This script should not run on login node!"; exit; fi

# no particle
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1
make clean
make -j

# no particle mpi
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --mpi=true
make clean
make -j

# particle
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --particle=true
make clean
make -j

# particle mpi
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --particle=true --mpi=true
make clean
make -j

# flt att
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --particle=true --par_attribute_flt=1
make clean
make -j

# int att
${PYTHON} configure.py --machine=eureka_intel --hdf5=true --gpu=true --fftw=FFTW3 \
                       --model=HYDRO --gravity=true --eos=ISOTHERMAL --barotropic=true \
                       --flu_scheme=MHM --flux=HLLC --passive=1 --particle=true --par_attribute_int=1
make clean
make -j

# feedback
${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=HYDRO --particle=true --gravity=true --feedback=true
make clean
make -j

# star formation
${PYTHON} configure.py --machine=eureka_intel --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true \
                       --model=HYDRO --particle=true --gravity=true --flu_scheme=MHM --flux=HLLC \
                       --passive=1 --par_attribute_flt=1 --dual=DE_ENPY --star_formation=true
make clean
make -j
