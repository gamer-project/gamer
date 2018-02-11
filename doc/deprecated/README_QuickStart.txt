
#====================
Download GAMER
#====================

      hg clone https://bitbucket.org/gamer_project/gamer
      cd gamer
      export GAMER=$(pwd)


#====================
Compilation
#====================
1. Please edit at least the following options and paths in ${GAMER}/src/Makefile

      SIMU_OPTION += -DGPU             # use GPU or not
      SIMU_OPTION += -DGPU_ARCH=KEPLER # GPU architecture

      CUDA_TOOLKIT_PATH :=             # useful only if GPU is enabled
      FFTW_PATH         :=             # useful only if GRAVITY is enabled
      MPI_PATH          :=             # useful only if SERIAL is disabled
      HDF5_PATH         :=             # useful only if SUPPORT_HDF5 is enabled

2. Currently GAMER only supports FFTW2. Please configure it with the
   floating-point type prefix enabled. For example,

      # installation path
      export FFTW_PATH=YOUR_FFTW_PATH

      # single precision
      ./configure --enable-mpi --enable-type-prefix --prefix $FFTW_PATH --enable-float
      make
      make install

      # double precision
      make clean
      ./configure --enable-mpi --enable-type-prefix --prefix $FFTW_PATH
      make
      make install

3. You may need to edit other compilation flags for the target test problem.
   Please refer to the README file associated with each test problem.
   For example, to run the AGORA_IsolatedGalaxy and ClusterMerger_vs_Flash,
   test problems, please turn on

      GRAVITY
      POT_SCHEME=SOR
      STORE_POT_GHOST
      UNSPLIT_GRAVITY
      PARTICLE

   For AGORA_IsolatedGalaxy, please also enable

      DUAL_ENERGY=DE_ENPY

4. Compile GAMER by

      cd ${GAMER}/src
      make clean
      make  # (or "make -j NThread" to compile in parallel using NThread threads)


#====================
Run GAMER
#====================
1. Runtime parameters files of various test problems are put at

      ${GAMER}/example/test_problem/

   To run, for example, the BlastWave test problem, please

      cd ${GAMER}/bin/run
      cp ${GAMER}/example/test_problem/Hydro/BlastWave/* .
      ./gamer >& log

      # to quickly visualize the density profile evolution
      gnuplot plot.gpt

2. By default, GAMER will use all CPU cores available if the compilation flag
   OPENMP is enabled
   --> Adjust OMP_NTHREAD in Input__Parameter if you want

3. To run the AGORA_IsolatedGalaxy and ClusterMerger_vs_Flash test problems,
   please first download the initial condition files using the bash script
   "download_ic.sh" associated with each of them.


#====================
Runtime parameters:
#====================
1. Input__Parameter: runtime parameters applied to all test problems
                     (e.g., BOX_SIZE, NX0_TOT[0/1/2], END_T, ... )
2. Input__TestProb:  problem-specific runtime parameters
3. Input__Flag_XXX:  thresholds of various grid refinement criteria


#====================
Log files
#====================
1. Record__Note       : compilation options and runtime parameters
2. Record__PatchCount : number of patches at each level
3. Record__Performance: simulation time-step and performance
4. Record__Dump       : physical times of all data dumps


#====================
Data analysis
#====================
1. Data_XXXXXX are the snapshots that can be loaded by yt.


#====================
Restart
#====================
1. Set "OPT__INIT == 2" in Input__Parameter
2. Choose the snapshot for restart by "ln -s Data_XXXXXX RESTART"
3. run GAMER


#====================
OpenMP optimization
#====================
1. Thread affinity:
   Ref: https://software.intel.com/en-us/node/522691
   Example:

      export KMP_AFFINITY=verbose,granularity=core,compact

   may be faster than

      export KMP_AFFINITY=verbose,granularity=core,scatter
