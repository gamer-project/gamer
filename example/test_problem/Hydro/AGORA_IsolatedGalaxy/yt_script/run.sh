NNODE=20
START=0
END=50
DELTA=1

mpirun -np $NNODE python plot_gas_profile.py      -s $START -e $END -d $DELTA 1>>log.gas.profile 2>&1
mpirun -np $NNODE python plot_gas_temp-vs-dens.py -s $START -e $END -d $DELTA 1>>log.gas.temp-vs-dens 2>&1
mpirun -np $NNODE python plot_gas_slice.py        -s $START -e $END -d $DELTA 1>>log.gas.slice 2>&1
mpirun -np $NNODE python plot_particle_mass.py    -s $START -e $END -d $DELTA 1>>log.par 2>&1

python plot_star_formation_rate.py
