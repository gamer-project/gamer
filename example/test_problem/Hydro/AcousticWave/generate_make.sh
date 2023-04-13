# This script should run in the same directory of configure.py (default directory is `src`)

PYTHON=python

${PYTHON} configure.py --cluster=eureka --flag=intel --model=HYDRO --double --timing
