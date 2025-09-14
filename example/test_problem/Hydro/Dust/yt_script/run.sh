START=0
END=50
DELTA=1

python metal_density.py    -s $START -e $END -d $DELTA 
python dust_density.py     -s $START -e $END -d $DELTA 
python dust_temp.py        -s $START -e $END -d $DELTA 
python gas_temp.py         -s $START -e $END -d $DELTA 