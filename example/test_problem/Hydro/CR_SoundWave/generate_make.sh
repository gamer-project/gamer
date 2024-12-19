# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --machine=eureka_intel --model=HYDRO --flu_scheme=MHM_RP --mhd=false --flux=HLLE \
                       --double=true --cosmic_ray=true --eos=COSMIC_RAY "$@"
#${PYTHON} configure.py --machine=eureka_intel --model=HYDRO --flu_scheme=MHM_RP --mhd=true --flux=HLLD \
#                       --double=true --cosmic_ray=true --eos=COSMIC_RAY "$@"
