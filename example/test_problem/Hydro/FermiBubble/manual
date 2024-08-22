# Simulation
#----------------------------------------------------------------------------------------------------
# 1. Generate makefile (double precision)
python3 configure.py --flu_scheme=MHM_RP --gravity=true --pot_scheme=SOR --store_pot_ghost=true --unsplit_gravity=false --slope=PLM --flux=HLLC --srhd=true --cosmic_ray=true --double=true --hdf5=true --nlevel=14 --fftw=FFTW2 --passive=3 --mpi=true --gpu=true

# 2. Compile GAMER
make clean
make -j
mv ../bin/gamer ../bin/gamer_double

# 3. Generate makefile (single precision)
python3 configure.py --flu_scheme=MHM_RP --gravity=true --pot_scheme=SOR --store_pot_ghost=true --unsplit_gravity=false --slope=PLM --flux=HLLC --srhd=true --cosmic_ray=true --double=false --hdf5=true --nlevel=14 --fftw=FFTW2 --passive=3 --mpi=true --gpu=true

# 4. Compile GAMER
make clean
make -j
mv ../bin/gamer ../bin/gamer_single

# 5. Copy simulation files
cd ../bin
cp -r ../example/test_problem/Hydro/FermiBubble ./
cd FermiBubble
cp ../gamer_* ./
sh download_ic.sh

# 6. Run the simulation (stage 1)
# 6-1. Change the parameters to below as follows:
   Input__Parameter:
      END_T        3.14872516940871e+01
      OUTPUT_DT    2
      OPT__INIT    1

   Input__TestProb:
      Jet_Fire     3

# 6-2. Execute gamer
mpirun -map-by ppr:2:socket:pe=8 --report-bindings ./gamer_double 1>>log 2>&1

# 7. Run the simulation (stage 2)
# 7-1. Change the parameters to below as follows:
   Input__Parameter:
      END_T     3.8e3
      OUTPUT_DT 200
      OPT__INIT 2

   Input__TestProb:
      Jet_Fire 0

# 7-2. Link the last file from stage 1
ln -fs Data_000016 RESTART

# 7-3. Execute gamer
mpirun -map-by ppr:2:socket:pe=8 --report-bindings ./gamer_single 1>>log 2>&1



# Analysis
#----------------------------------------------------------------------------------------------------
# Please execute the following commands before the analysis
cd plot_scripts
ln -s ../R12
ln -s ../../../tool/analysis/PerspectiveProjection


## Slice plot
python plot_slice.py -s 0 -e 35


## Profile plot along z axis of the lower jet
python plot_profile.py -s 0 -e 35


## Generate the fixed resolution data (Data_000035)
python AMR2FR.py


## Generate x-ray perspective projection map (Data_000035)
# 1. If you do not have FRB_Data_000035.h5, please execute the following command
python AMR2FR.py

# 2. Compile the projection tool
cd PerspectiveProjection
make clean && make CFLAGS+=-DNUM_THREADS=32 CFLAGS+=-DXRAY_EROSITA
cd ../
cp PerspectiveProjection/bin/Project ./

# 3. Execute the projection tool
./Project FRB_Data_000035.h5

# 4. Plot the map
python plot_map.py -t x_ray


## Generate x-ray profile (Data_000035)
# 1. If you do not have FRB_Data_000035.h5, please execute the following command
python AMR2FR.py

# 2. Plot the x-ray profile
python plot_xray_profile.py

## Generate gamma ray perspective projection map (Data_000035)
# 1. If you do not have FRB_Data_000035.h5, please execute the following command
python AMR2FRB.py

# 2. Compile the projection tool
cd PerspectiveProjection
make clean && make CFLAGS+=-DNUM_THREADS=32 CFLAGS+=-DLEPTONIC_GAMMARAY
cd ../
cp PerspectiveProjection/bin/Project ./

# 3. Execute the projection tool
./Project FRB_Data_000035.h5 100e9 1e6 2.4 R12/robitaille_DL07_PAHISMMix.dat
# NOTE : arg1: name of fixed resolution data in Step1
#        arg2: observed photon energy (eV)
#        arg3: the cut-off Lorentz factor of CR
#        arg4: spectral index of CR
#        arg5: the path of ISRF data

# 4. Plot the map
python plot_map.py -t gamma_ray


## Generate spectrum
# 1. If you do not have FRB_Data_000035.h5, please execute the following command
python AMR2FR.py

# 2. Compile the projection tool
# 2.1 synchrotron emissivities
cd PerspectiveProjection
make clean && make CFLAGS+=-DNUM_THREADS=32 CFLAGS+=-DSYNCHROTRON
cd ../
cp PerspectiveProjection/bin/Project synchrotron_spectrum/

# 2.2 gamma-ray emissivities
cd PerspectiveProjection
make clean && make CFLAGS+=-DNUM_THREADS=32 CFLAGS+=-DLEPTONIC_GAMMARAY
cd ../
cp PerspectiveProjection/bin/Project gamma_ray_spectrum/

# 3. Calculate spectrums
# 3-1. synchrotron emissivities
cd synchrotron_spectrum
ln -s ../FRB_Data_000035.h5
ln -s ../../R12
sh get_spectrum.sh
cd ../

# 3-2. gamma ray emissivities
cd gamma_ray_spectrum
ln -s ../FRB_Data_000035.h5
ln -s ../../R12
sh get_spectrum.sh
cd ../

# 4. Plot the spectrum
python plot_spectrum.py
