# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --flu_scheme=MHM_RP --gravity=true --pot_scheme=SOR --store_pot_ghost=true \
                       --unsplit_gravity=false --slope=PLM --flux=HLLC --srhd=true --cosmic_ray=true \
                       --double=true --hdf5=true --nlevel=14 --fftw=FFTW2 --passive=3 --mpi=true --gpu=true "$@"
