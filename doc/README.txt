
#Obtain code
#==============================================================
hg clone https://bitbucket.org/gamer_project/gamer
cd gamer
export GAMER_PATH=$(pwd)


#Set up the "ClusterMerger" test
#==============================================================
#NOTE: Currently it only supports a single cluster in hydrostatic equilibrium.
#      All other test problems haven't been updated for a while and hence they cannot
#      be compiled directly. I can enable them easily if you want to test any of them.
#      In addition, currently there is no way to run different test problems
#      in a flexible way, and the script "Copy.Single.sh" will simply overwrite all
#      files necessary for this test. This is in my TODO list to make it more flexible.
cd ${GAMER_PATH}/test_problem/Model_ParticleOnly/ClusterMerger
sh Copy.Single.sh
cd ${GAMER_PATH}/src


#Edit Makefile
#==============================================================
#NOTE: Please edit the following options and paths. Default is to turn off GPU and MPI
#      (which are controlled by -DGPU, -DGPU_ARCH, -DSERIAL -DLOAD_BALANCE=HILBERT),
#      and hence you don't have to set CUDA_TOOLKIT_PATH and MPI_PATH.
 
#SIMU_OPTION += -DINTEL  # use Intel compiler (otherwise GNU compiler is used)
 
#CUDA_TOOLKIT_PATH := /usr/local/cuda/7.0
#FFTW_PATH         := /projects/ncsa/grav/softwares/fftw/2.1.5
#MPI_PATH          := /usr/local/mpi/openmpi-1.8.4-intel-15.0
#HDF5_PATH         := /projects/ncsa/grav/softwares/miniconda2
#GSL_PATH          := /usr/local/gsl/gsl-1.16


#Compile
#==============================================================
make clean;
make
#(or "make -j NThread" to compile in parallel using NThread threads)


#Edit runtime parameters and run
#==============================================================
#NOTE: Input__Parameter is the runtime parameter file. Currently it is hard coded.
#      So please don't comment out, delete or add any line. Edit OMP_NTHREAD to set
#      the number of OpenMP threads. You could leave all other parameters unchanged.
#      Record__* are different log files for this simulation. It takes ~10
#      minuites to run on the Campus cluster using 20 OpenMP threads. You should get
#      53 snapshots in total (Data_000000 ... Data_000052).
cd ${GAMER_PATH}/bin/run
# edit Input__Parameter if necessary
./Dizzy 1>stdout 2>stderr


#Quick visualization
#==============================================================
#NOTE: Data_XXXXXX are the snapshots that can be loaded by YT.

#plot__conservation.gpt           : plot mass and energy conservations
#plot__particle_single-cluster.gpt: plot particle distribution

#Please don't use plot__density-profile.gpt and plot__temperature-profile.gpt.
#They need aother post-processing tool to calculate profiles first (but you
#can use YT to do everything!)


#Initial condition set-up
#==============================================================
#Gas:
#   Edit "src/Init/Init_TestProb.cpp --> Par_TestProbSol_ClusterMerger".
#   This is actually a function pointer replacing
#   "src/Model_Hydro/Hydro_Init_StartOver_AssignData.cpp --> Init_Function_User"
#   Please set "fluid[DENS/MOMX/MOMY/MOMZ/ENGY]" as a function of the input x/y/z/Time.
 
#Particle:
#   Edit "src/Particle/Par_Init_ByFunction.cpp"
#   Please set "amr->Par->Mass/PosX/PosY/PosZ/VelX/VelY/VelZ" for all particles.
#   The total number of particles at ALL MPI ranks and THIS MPI rank are given by
#   "amr->Par->NPar_Active_AllRank" and "NPar_Active", respectively.

#Code units:
#   Edit "Input__Parameter --> UNIT_L/M/T/V/D" in CGS.
#   Some physical constants and common units are predefined in "include/PhysicalConstant.h" in CGS.

