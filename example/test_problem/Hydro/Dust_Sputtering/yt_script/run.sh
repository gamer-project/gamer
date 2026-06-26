# Example usage of yt analysis scripts

START=0
END=50
DELTA=1

# EDOTMODE has two options:
#   edot_0     -> corresponds to edot(i) = 0
#   edot_const -> corresponds to edot(i) = -k_sup * e
#
# Change the EDOTMODE value as needed

EDOTMODE=edot_0    # or edot_const

python3 metal_density.py -s $START -e $END -d $DELTA
python3 dust_density.py  -option $EDOTMODE
python3 gas_density.py   -s $START -e $END -d $DELTA
python3 gas_temp.py      -s $START -e $END -d $DELTA -option $EDOTMODE
