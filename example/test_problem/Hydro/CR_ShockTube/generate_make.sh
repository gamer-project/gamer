# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM_RP --mhd=false --flux=HLLE \
                       --cosmic_ray=true --eos=COSMIC_RAY "$@"
#${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM_RP --mhd=true --flux=HLLD \
#                       --cosmic_ray=true --eos=COSMIC_RAY "$@"
