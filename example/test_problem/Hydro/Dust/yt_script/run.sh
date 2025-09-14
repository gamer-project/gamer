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

python metal_density.py -s $START -e $END -d $DELTA
python dust_density.py  -s $START -e $END -d $DELTA -option $EDOTMODE
python dust_temp.py     -s $START -e $END -d $DELTA
python gas_temp.py      -s $START -e $END -d $DELTA -option $EDOTMODE
