
       *****************************************************************
       **                                                             **
       **    GAMER : GPU-accelerated Adaptive MEsh Refinement code    **
       **                                                             **
       *****************************************************************

-------------------------------------------------------------------------------


This document contains the instructions for using GAMER. It consists
of five parts :
   I.   Makefile
   II.  Perform Simulation
   III. Data Structure
   IV.  Notice
   V.   New Features and Revision History


==========================================================================
I. Makefile
==========================================================================

1. To compile the full version of GAMER, it requires three libraries to be
   properly installed :
   a. CUDA Toolkit
   b. FFTW (version 2.1.5 with the single/double-precision prefixes)
   c. MPI

   Please set the location of each of these three libraries in the
   corresponding variable in the Makefile :

   --> CUDA_TOOLKIT_PATH, FFTW_PATH, MPI_PATH

   For example, you might want to set "CUDA_TOOLKIT_PATH := /usr/local/cuda".
   You might also want to replace

      LIB := -L$(CUDA_TOOLKIT_PATH)/lib64:

   by

      LIB := -L$(CUDA_TOOLKIT_PATH)/lib32

   if you are using the 32-bit CUDA Toolkit.

**Important** For the CPU-only mode (by turning off the option "GPU"), no CUDA
              libraries are required. For the purely hydrodynamic mode (by
              turning off the option "GRAVITY"), no FFTW libraries are
              required. For the serial mode (by turning on the option
              "SERIAL"), no MPI libraries are required.


2. In the Makefile, currently there are 6 options for choosing different
   numerical schemes:

   a. FLU_SCHEME :
         Choose the hydrodynamic scheme.
            RTVD   -> relaxing TVD scheme
            WAF    -> weighted-average-flux scheme
            MHM    -> MUSCL-Hancock method
            MHM_RP -> MUSCL-Hancock method with Riemann predicion
                      (the VL scheme adopted in Athena)
            CTU    -> corner-transport-upwind scheme

   b. LR_SCHEME :
         Choose the spatial data reconstruction scheme (only useful in
         "MHM/MHM_RP/CTU" schemes).
            PLM -> piecewise-linear method
            PPM -> piecewise-parabolic method

   c. RSOLVER :
         Choose the Riemann solver for WAF/MHM/MHM_RP/CTU schemes (WAF only
         supports EXACT/ROE solvers).
            EXACT -> exact Riemann solver (vacuum state is NOT included yet)
            ROE   -> Roe's Riemann solver
            HLLE  -> Harten, Lax, van Leer and Einfeldt's Riemann solver
            HLLC  -> Extension of HLLE solver (including contact wave) by Toro

   d. GRAVITY :
         Simulate the hydrodynamic + self-gravity system. If this option is
         turned off, GAMER assumes a purely hydrodynamic system.

   e. INDIVIDUAL_TIMESTEP :
         Perform simulations with the individual time-step scheme, in which
         the time-step at level Lv is twice smaller than that at level Lv-1.
         At each global step, the level Lv will contain 2^(Lv+1) sub-steps.
         If this option is turned off, the shared time-step scheme is adopted.

   f. COMOVING :
         Perform cosmological simulations in the comoving coordinates. Please
         refer to the GAMER reference paper for the definitions of physical
         variables in the comoving coordinates.


3. In the Makefile, currently there are 11 options for optimization and
   compilation :

   a. GPU :
         Use the GPU solvers instead of CPU solvers. If this option is turned
         off, GAMER will use CPUs only. Therefore, to measure the performance
         improvement by using GPUs, one can perform two simulations with and
         without this option and compare the elapsed times.

   b. GAMER_OPTIMIZATION :
         Optimize the GAMER performance. Currently this option only has little
         effect on performance.

   c. GAMER_DEBUG :
         Enable more error checks. Currently this option is not fully
         supported. Note that the round-off error in debug mode may be
         different from that in non-debug mode.

   d. TIMING :
         Measure the elapsed wall-clock times of different parts in GAMER.
         The timing results will be recorded in the file "Record__Timing".

   e. TIMING_SOLVER :
         Measure the elapsed wall-clock times of different parts of the
         fluid/Poisson/gravity solvers. Note that this option will disable
         the concurrent execution between CPU and GPU, and hence will
         deteriorate the overall performance. Accordingly, we suggest that
         this option should never be turned on unless detailed timing analyses
         are required.

   f. INTEL :
         Use the intel compiler (default is GNU compiler).

   g. FLOAT8 :
         Use double precision for simulations. If this option is turned off,
         single precision is adopted.

**Important** Currently the double-precision GPU Poisson solver is only
              supported in Fermi GPUs. Therefore, the options "GPU",
              "GRAVITY", and "FLOAT8" can not be turned on at the same time
              without the option "FERMI".

**Important** Even in Fermi GPUs, the double-precision performance is still
              much lower than the single-precision performance in GAMER. This
              issue will be addressed in the future release.

   h. SERIAL :
         Compile GAMER in the serial mode, in which only one GPU can be used.
         No MPI libraries are required in this mode. OpenMP is still supported
         in this mode.

   i. OOC : (NOT SUPPORTED YET)
         Enable the out-of-core computation.

   j. OPENMP :
         Enable the OpenMP parallelization. Use the parameter "OMP_NTHREAD" in
         the parameter file to choose the number of OpenMP threads per MPI
         process.

**Imporant** Turn on this option whenever it is possible. Significant
             performance improvement can be achieved even when the option "GPU"
             is turned on.

   k. FERMI :
         Enable performance optimization in Fermi GPUs. Running GAMER with
         this option in non-Fermi GPUs will produce error messages.


4. Currently there are 7 symbolic constants defined in the Makefile :

   a. NCOMP :
         The number of physical attributes stored in each cell (excluding
         the gravitational potential). This number is designed for making
         GAMER more general-purpose in the future. However, currently it is
         restricted to be 5.

   b. PATCH_SIZE :
         The number of cells in each spatial direction of a single patch.
         Currently it is restricted to be 8.

   c. MAX_PATCH :
         The maximum number of patches at each level (in one MPI process).
         It is used to prevent from running out of CPU memory. If the number
         of patches at any level has exceeded this number, an error message
         will be generated and the program will be terminated.

   d. NLEVEL :
         The maximum number of levels (including the root level and the
         refinement levels). For example, if NLEVEL == 3, we have one root
         level and two refinement levels (level = 0,1,2).

   e. POT_GHOST_SIZE :
         The number of ghost zones on each side in each direction for the
         Poisson solver. It is a free positive integer, which defines the
         number of cells overlapping between adjacent patches when
         evaluating the potential. Larger value will make the potential more
         smooth across the patch boundary. If the option "GPU" is turned on,
         currently this value is restricted to "1 ~ 5". For the CPU-only mode,
         it is restricted to "1 ~ PATCH_SIZE". We suggest to set
         "POT_GHOST_SIZE = 4 or 5" for the balance between performance and
         accuracy.

   f. GRA_GHOST_SIZE :
         The number of ghost zones on each side in each direction for
         evaluating the potential gradient. This value should be set to 1(2)
         if a 3(5)-point stencil is adopted for evaluating the potential
         gradient. The 5-point stencil can be applied by turning on the option
         "OPT__GRA_P5_GRADIENT" in the parameter file "Input__Parameter".

   g. FLU_GHOST_SIZE :
         The number of ghost zones on each side in each direction for the
         hydrodynamic solver. This variable is determined by the hydrodynamic
         scheme and will be set automatically.


5. Other user-defined variables related to the compilation.
   a. CXX :
         The name of the c++ compiler.

   b. LIB :
         The paths and names of the libraries to be linked.

   c. CXXWARN_FLAG :
         The warning flags for the c++ compiler.

   d. CXXFLAG :
         The flags for the c++ compiler.

   e. NVCCFLAG :
         The flags for the nvcc compiler.

**Important** If the option "FERMI" is not turned on, the program will be
              compiled for GPU architectures "sm_10, sm_13, and sm_20".


6. If the program compiled successfully, an executable file named "Dizzy" will
   be created and all object files will be placed in the directory
   "src/Object". A copy of the executable file will put in the directory
   "bin/Run".



==========================================================================
II. Perform Simulation
==========================================================================

1. To perform simulation, GAMER requires at most seven input files, which must
   be put in the same directory as the executable file. All input files have
   the prefix "Input__", Examples of these input files are placed in the
   directory "input_example". Below we give detailed description for each of
   these input files.

-------------------------------------------
   a. Input__Parameter :
-------------------------------------------

      This input file sets the simulation parameters and must be provided
      in all cases. In the following, we explain the usage of each of these
      parameters.

      *****************************************************************
      Parameters related to the simulation scale :
      *****************************************************************

      BOX_SIZE :
         The physical size of simulation box in the longest dimension. For
         example, if NX0_TOT[0/1/2] = 32/64/96, the BOX_SIZE will correspond
         to the simulation box size in the z direction (in this example, the
         box sizes in x and y directions will be set to BOX_SIZE/3 and
         BOX_SIZE*2/3, automatically).

      NX0_TOT[0] :
         The total number of root-level grids in the x direction. It must be a
         multiple of "2*PATCH_SIZE". Moreover, in multi-GPU cases, since each
         MPI process must have at least two root-level patches in each spatial
         direction, NX0_TOT[0] must be a multiple of
         "2*PATCH_SIZE*MPI_NRANK_X[0]", where MPI_NRANK_X[0] is the number of
         MPI processes in the x direction.

         For example, if we set "PATCH_SIZE=8, NLEVEL=3, NX0_TOT[0]=32", the
         spatial resolution in the x direction is "NX0_TOT[0] = 32" in the
         root level (level 0) and "NX0_TOT[0]*2^(NLEVEL-1) = 32*2^(3-1) = 128"
         in the maximum refinement level (level 2). In this case, we can only
         have MPI_NRANK_X[0] = 1 or 2, while MPI_NRANK_X[0] > 2 is prohibited
         since each MPI process must have at least "2*PATCH_SIZE = 16"
         root-level cells in the x direction.

         Also note that in GAMER the refinement ratio between levels (the ratio
         of the grid resolutions between adjacent levels) is limited to 2.

      NX0_TOT[1] :
         The total number of root-level grids in the y direction.

      NX0_TOT[2] :
         The total number of root-level grids in the z direction.

**Important** The grid size at the root level will be
              "BOX_SIZE / max( NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] )".

      MPI_NRANK :
         The total number of MPI ranks. Currently we assume that each MPI
         process will use only one GPU. The option "OPT__GPUID_SELECT" can be
         used to choose different modes of GPU ID selection.

      MPI_NRANK_X[0] :
         The number of MPI ranks in the x direction. For example, if we have
         "NX0_TOT[0]=32", setting "MPI_NRANK_X[0]=2" means that each GPU
         (MPI process) will have "NX0_TOT[0] / MPI_NRANK_X[0] = 16" root-level
         cells in the x direction.

      MPI_NRANK_X[1] :
         The number of MPI ranks in the y direction.

      MPI_NRANK_X[2] :
         The number of MPI ranks in the z direction.

**Important** One must have "MPI_NRANK_X[0]*MPI_NRANK_X[1]*MPI_NRANK_X[2]
              == MPI_NRANK", otherwise an error message will be generated.

      OMP_NTHREAD :
         The number of OpenMP threads (useless if the option "OPENMP" is
         turned off in the Makefile).

**Important** If a negative number is provided here, this parameter will be set
              to the default value [omp_get_max_threads].

      END_T :
         The physical time to terminate the simulation.

      END_STEP :
         The maximum number of the evolution time-steps. Note that either
         "Time > END_T" or "Step > END_STEP" will terminate the simulation.


      *****************************************************************
      Parameters related to the cosmological simulation :
      *****************************************************************

      A_INIT :
         The initial value of the scale factor for cosmological simulations.
         If the option "COMOVING" is not turned on in the Makefile, this
         parameter is useless and the initial time is always set to zero. Note
         that for cosmological simulations, the physical time recorded during
         the run is actually the "scale factor".

      OMEGA_M0 :
         The matter density at the present time for cosmological simulations.
         This parameter is useless if the option "COMOVING" is turned off in
         the Makefile.


      *****************************************************************
      Parameters related to the simulation time-step  :
      *****************************************************************

      DT__CFL :
         The CFL condition parameter to determine the evolution time-step.

**Important** " If a negative number is provided here, this parameter will be
                set to the default value [0.5].

      DT__ACC :
         The parameter used to determine the evolution time-step by estimating
         the maximum potential gradient.

      DT__MAX_DELTA_A :
         The parameter used to determine the evolution time-step by restricting
         the percentage of scale factor updated at one time-step. For example,
         if we set "DT__MAX_DELTA_A=1.e-2, A_INIT=1.e-3", the maximum time-step
         at the first step is restricted to be smaller than
         "A_INIT*DT__MAX_DELTA_A = 1.e-5". This parameter is useless if the
         option "COMOVING" is turned off.

      OPT__GPU_DT_CFL :
         Option : use GPUs to evaluate the CFL condition. Currently this option
         is not supported for the individual time-step scheme.

      OPT__GPU_DT_ACC : (NOT SUPPORTED YET)
         Option : use GPUs to evaluate the maximum potential gradient for the
         time-step estimation.

      OPT__RECORD_DT :
         Record the information of time-step determination in the file
         "Record__TimeStep".


      *****************************************************************
      Parameters related to the domain refinement :
      *****************************************************************

      REGRID_COUNT :
         The frequency of performing the refine operation. The patches will be
         reconstructed/removed every "REGRID_COUNT" step. Note that if the
         individual time-step scheme is enabled, this parameter is refer to the
         "sub-step" at each level rather than the global time-step.

      FLAG_BUFFER_SIZE :
         The size of flag buffer. Please refer to the GAMER reference paper
         for the usage of this parameter. Generally speaking, the larger the
         FLAG_BUFFER_SIZE, the more region will be refined.

      OPT__FLAG_RHO :
         Option : use the mass density as one of the flag criteria. The flag
         threshold at each level will be loaded from the file "Input__Flag_Rho".

      OPT__FLAG_RHO_GRADIENT :
         Option : use the mass density gradient as one of the flag criteria.
         The flag threshold at each level will be loaded from the file
         "Input__Flag_RhoGradient".

      OPT__FLAG_RHO_GRADIENT :
         Option : use the pressure gradient as one of the flag criteria. The
         flag threshold at each level will be loaded from the file
         "Input__Flag_PresGradient".

      OPT__FLAG_USER :
         Option : use the user-defined function as one of the flag criteria.
         The flag threshold at each level will be loaded from the file
         "Input__Flag_User".

**Important** To use this option, please complete the function
              "src/Refine/Flag_UserCriteria.cpp".

      OPT__PATCH_COUNT :
         Option : record the number of patches at each level in the file
         "Record__PatchCount".


      *****************************************************************
      Parameters related to the hydrodynamic simulation :
      *****************************************************************

      GAMMA :
         The ratio of the specific heats.

      FLU_GPU_NPGROUP :
         The number of patch groups for the GPU hydrodynamic solver. Please
         refer to the GAMER reference paper for the definition of patch group.

**Important** If a negative number is provided here, this parameter will be set
              to the default value [2*number of multiprocessors*GPU_NSTREAM].

      GPU_NSTREAM :
         The number of CUDA streams used by the GPU
         hydrodynamic/Poisson/gravity solvers. The CUDA streams are used to
         hide the data transfer time between CPU and GPU. Please refer to the
         GAMER paper for a more detailed description of the usage of CUDA
         streams in GAMER. Also note that for performance consideration,
         "FLU_GPU_NPGROUP/GPU_NSTREAM" and "POT_GPU_NPGROUP/GPU_NSTREAM" should
         be multiples of the number of multiprocessors in one GPU. For example,
         for T10 GPU which has 30 multiprocessors, one might want to set
         "FLU_GPU_NPGROUP=240, GPU_NSTREAM=4", so that each multiprocessor will
         work on 2 patch groups at a time.

**Important** If a negative number is provided here, this parameter will be set
              to the default value, which is (1 / 4) if the asynchronous memory
              copy is (not supported / supported) in the adopted GPU.

      MINMOD_COEFF :
         The coefficient of the generalized MinMod limiter for the data
         reconstruction in MHM/MHM_RP/CTU schemes. This parameter must be
         within the range [1.0~2.0] to ensure the TVD condition. Note that this
         parameter will be useless if "OPT__LR_LIMITER != 1 && != 3".

      EP_COEFF :
         The coefficient of the extrema-preserving limiter. This parameter
         will be useless if "OPT__LR_LIMITER != 4".

      OPT__LR_LIMITER :
         Option : the slope limiter for the data reconstruction in
         MHM/MHM_RP/CTU schemes.

            0 -> van Leer limiter
            1 -> generalized MinMod limiter
                 (please set the parameter "MINMOD_COEFF" for this limiter)
            2 -> van Albada limiter
            3 -> generalized MinMod + van Leer limiter
                 (please set the parameter "MINMOD_COEFF" for this limiter)
            4 -> extrema-preserving limiter (not well tested yet)

      OPT__WAF_LIMITER :
         The flux limiter in WAF scheme.

            0 -> SuperBee limiter
            1 -> van Leer limiter
            2 -> van Albada limiter
            3 -> MinBee limiter

      OPT__FIXUP_FLUX :
         Option : perform the flux correction. It ensures that the conservative
         variables are conserved even for patches adjacent to the coarse-fine
         boundaries (assuming no self-gravity).

      OPT__FIXUP_RESTRICT :
         Option : perform the restrict correction. It ensures that the
         conservative variables of a parent cell are equal to the averages of
         the values of its eight child cells.


      *****************************************************************
      Parameters related to the self-gravity simulation :
      *****************************************************************

      NEWTON_G :
         The gravitational constant. Note that this value will be overwritten
         in cosmological simulations, in which the gravitation constant is
         determined by other pre-defined constants.

      SOR_OMEGA :
         The over-relaxation parameter for the SOR scheme. This parameter will
         be set to the default value if a negative value is found in this file.

      SOR_MAX_ITER :
         The maximum number of iterations for the SOR scheme. This parameter
         will be set to the default value if a negative value is found in this
         file.

      SOR_MIN_ITER :
         The minimum number of iterations for the SOR scheme. This parameter
         will be set to the default value if a negative value is found in this
         file.

      POT_GPU_NPGROUP :
         The number of patch groups for the GPU Poisson/gravity solver. Please
         refer to the GAMER paper for the definition of the patch group.

**Important** If a negative number is provided here, this parameter will be set
              to the default value, which is
              "2 * number of multiprocessors * GPU_NSTREAM".

      OPT__GRA_P5_GRADIENT :
         Option : use a five-point stencil to evaluate the potential gradient.
         To turn on this option, one must set "GRA_GHOST_SIZE=2" in the
         Makefile.


      *****************************************************************
      Parameters related to the initialization :
      *****************************************************************

      OPT__INIT :
         Option : initialization mode. Currently there are three modes
         supported in GAMER.
         (i  ) Mode 0 : Construct the initial condition by the function
                        "Init_AssignData".
         (ii ) Mode 1 : Load any previous data dump as a restart file. The
                        restart file must be named as "RESTART". The typical
                        method is to make a symbolic link, for example,
                        "ln -s Data_000010 RESTART"
         (iii) Mode 2 : Construct the initial condition by loading a
                        uniform-mesh data named "UM_START", which should be a
                        binary file. Two parameters "OPT__UM_START_LEVEL" and
                        "OPT__UM_START_NVAR" are used to further specify the
                        input uniform-mesh data. See below for a more detailed
                        description.

**Important** For mode 0, one must complete the function "Init_AssignData".

      OPT__RESTART_FORMAT :
         Option : the data format of the restart file. Please set it to 1.
         Options 0 and 2 are out-of-date and will be removed in the future
         release.

      OPT__UM_START_LEVEL :
         The refinement level of the UM_START file. For example, if
         "NX0_TOT[0]=NX0_TOT[1]=NX0_TOT[2]=64, and OPT__UM_START_LEVEL=2",
         the UM_START file should contain
         "   OPT__UM_START_NVAR*( (64*2^(OPT__UM_START_LEVEL) )^3
           = OPT__UM_START_NVAR*256^3 " data.
         The function "Init_UM" will construct the initial condition from
         level OPT__UM_START_LEVEL down to level 0, and patches failed to
         fulfill the refinement criteria will be removed. Accordingly,
         depending on the refinement criteria, the simulation domain may not
         be fully refined to level "OPT__UM_START_LEVEL". No patch will be
         constructed at level > OPT__UM_START_LEVEL.

**Important** Currently, only the density (option OPT__FLAG_RHO) is adopted as
              the refinement criterion for constructing the initial condition
              in this mode (OPT__INIT == 2).

      OPT__UM_START_NVAR :
         The number of physical variables in each cell stored in the file
         "UM_START". Currently, this value must be set to 1 or NCOMP (5).
         (i ) OPT__UM_START_NVAR == 1 :
              The input uniform-mesh data only store the mass density. Velocity
              is initialized as zero and sound speed is initialized as 1.
         (ii) OPT__UM_START_NVAR == NCOMP (5):
              The input uniform-mesh data stores NCOMP data at each cell. The
              data should be organized in the following order (in binary form).

              Variable 0       at cell (0,0,0) <-- (x,y,z)
              Variable 1       at cell (0,0,0)
                  ...   ...   ...   ...   ...
              Variable NCOMP-1 at cell (0,0,0)
              Variable 0       at cell (1,0,0)
              Variable 1       at cell (1,0,0)
                  ...   ...   ...   ...   ...
              Variable NCOMP-1 at cell (1,0,0)
                  ...   ...   ...   ...   ...
                  ...   ...   ...   ...   ...
              Variable NCOMP-1 at cell ( NX0_TOT[0]*2^OPT__UM_START_LEVEL - 1,
                                         NX0_TOT[1]*2^OPT__UM_START_LEVEL - 1,
                                         NX0_TOT[2]*2^OPT__UM_START_LEVEL - 1 )

**Important** Variable (0,1,2,3,4)
              = (mass density, momentum density x,y,z, total energy density).
              The pseudo code below gives an example of producing the
              UM_START file:

              for (z=0; z<NX0_TOT[2]*2^OPT__UM_START_LEVEL - 1; z++)
              for (y=0; y<NX0_TOT[1]*2^OPT__UM_START_LEVEL - 1; y++)
              for (x=0; x<NX0_TOT[0]*2^OPT__UM_START_LEVEL - 1; x++)
              for (v=0; v<NCOMP; v++)
                 fwrite( &Array[z][y][x][v], sizeof(float), 1, File );

      OPT__INIT_RESTRICT :
         Option : after constructing the initial condition from level 0 to
         level NLEVEL-1, perform the restrict (average) operation to ensure
         that the data of any cell at level "lv" is equal to the spatial
         average of the data of its eight child cells at level "lv+1".

      OPT__GPUID_SELECT :
         Option : the GPU ID selection mode. Currently three modes are
         supported in GAMER.
         (i  ) Mode -2 : Choose the GPU ID by the CUDA libraries automatically.

**Important** One should turn on the "exclusive" compute mode for
              "OPT__GPUID_SELECT = -2" to ensure that different MPI processes
              will use different GPUs.

         (ii ) Mode -1 : Choose the GPU ID according to the MPI rank.
                         --> GPU ID = MPI_Rank % ( number of GPUs reported by
                                                   cudaGetDeviceCount )
                         For example, if there are 2 GPUs in each node, we
                         will have
                         MPI_Rank = (0,1,2,3...) <-> GPU ID = (0,1,0,1...)

         (iii) Mode >= 0 : Set the GPU ID equal to the input number
                           (GPU ID = Mode).


      *****************************************************************
      Parameters related to the spatial and temporal interpolation :
      *****************************************************************

      OPT__INT_TIME :
         Option : perform the temporal interpolation in the individual
         time-step scheme.

      OPT__FLU_INT_SCHEME :
         Option : the interpolation scheme to set the ghost-zone fluid data
         for the hydrodynamic solver.

      OPT__POT_INT_SCHEME :
         Option : the interpolation scheme to set the ghost-zone potential data
         for the Poisson solver.

      OPT__RHO_INT_SCHEME :
         Option : the interpolation scheme to set the ghost-zone density data
         for the Poisson solver.

      OPT__GRA_INT_SCHEME :
         Option : the interpolation scheme to set the ghost-zone potential data
         for the gravity solver (evaluating the potential gradient).

      OPT__REFINE_FLU_INT_SCHEME :
         Option : the interpolation scheme to set the fluid data for the
         newly-allocated patches during the refine operation.

      OPT__REFINE_POT_INT_SCHEME :
         Option : the interpolation scheme to set the potential data for the
         newly-allocated patches during the refine operation.

         Interpolation schemes:
            Mode 0 : Interpolation with the central limiter
            Mode 1 : Interpolation with the Min-Mod limiter
            Mode 2 : Interpolation with the van-Leer limiter
            Mode 3 : Conservative quadratic interpolation
            Mode 4 : Quadratic interpolation

**Important** For "OPT__FLU_INT_SCHEME" and "OPT__REFINE_FLU_INT_SCHEME",
              currently we suggest to use the linear interpolation with the
              Min-Mod scheme (set options = 1) to maintain both the
              conservation and monotonicity of data. More accurate
              interpolation schemes will be supported in the future version.


      *****************************************************************
      Parameters related to the data output :
      *****************************************************************

      OPT__OUTPUT_TOTAL :
         Option : output all simulation data (in binary form). The output files
         are named as "Data_XXXXXX" and can be used as restart files for
         "OPT__INIT == 1". The file "Record__Dump" will be created, which will
         record the physical time and simulation step of each data dump.

**Important** Please set "OPT__OUTPUT_TOTAL = 2" if you want to output the
              binary data. "OPT__OUTPUT_TOTAL = 1" is out-of-date and will be
              removed in the future release.

      OPT__OUTPUT_PART :
         Option : output a portion of simulation data (in text form). Currently
         GAMER supports 6 modes,
            Mode 0 : no output
            Mode 1 : output a single xy slice
            Mode 2 : output a single yz slice
            Mode 3 : output a single xz slice
            Mode 4 : output a line in the x direction
            Mode 5 : output a line in the y direction
            Mode 6 : output a line in the z direction
         No interpolation is performed, so that the output data may not be
         uniformly distributed. The output format is as following.
         "X, Y, Z, Density, Momentum X/Y/Z, Energy, Pressure, Potential"

      OPT__OUTPUT_BASE :
         Option : only output the root-level data for the option
         "OPT__OUTPUT_PART". Data output by this option can be plotted using
         the gnuplot tool "GAMER_2DBaseData2gnuplot.1.0".

      OPT__OUTPUT_POT :
         Option : output the potential data. Since that we can always evaluate
         potential from the density distribution, this option is necessary only
         when one wants to analyze the potential distribution. A data dump
         without storing the potential data can still be used as the restart
         file.

      OPT__OUTPUT_MODE :
         Option : the criterion of data dump. Currently three modes are
         supported :
         (i  ) Mode 0 : output data every "OUTPUT_STEP" step.
         (ii ) Mode 1 : output data every "OUTPUT_DT" time interval.
         (iii) Mode 2 : output data according to the physical time recorded
                        in the file "Input__DumpTable". Please refer to the
                        section "Input__DumpTable" for a more detailed
                        description of the usage of the dump table.

      OUTPUT_STEP :
         Output data every "OUTPUT_STEP" step.
         This parameter is useful only when "OPT__OUTPUT_MODE == 0"

      OUTPUT_DT :
         Output data every "OUTPUT_DT" time interval.
         This parameter is useful only when "OPT__OUTPUT_MODE == 1"

      OUTPUT_PART_X :
         The x coordinate for the option "OPT__OUTPUT_PART == 2/5/6".
         It must be smaller than the simulation box size in the x direction.

      OUTPUT_PART_Y :
         The y coordinate for the option "OPT__OUTPUT_PART == 3/4/6".
         It must be smaller than the simulation box size in the y direction.

      OUTPUT_PART_Z :
         The z coordinate for the option "OPT__OUTPUT_PART == 1/4/5".
         It must be smaller than the simulation box size in the z direction.

      OPT__VERBOSE :
         Option : output the detail of simulation progress.


      *****************************************************************
      Parameters related to the simulation check routines :
      *****************************************************************

      OPT__CK_REFINE :
         Option : verify that all patches fulfilling the refinement criteria
         are properly refined. Note that this check may fail due to the
         proper-nesting constraint. Also note that currently only the density
         flag criterion is checked in this option (OPT__FLAG_RHO).

      OPT__CK_PROPER_NESTING :
         Option : verify that the proper-nesting constraint is satisfied. This
         check should ALWAYS pass.

      OPT__CK_CONSERVATION :
         Option : check the conservation of fluid conservative variables. Note
         that in a hydrodynamic + self-gravity system, the integration scheme
         currently adopted in GAMER only ensures the conservation of "mass density".

      OPT__CK_NEGATIVE :
         Option : check whether there are cells with negative density or pressure.

      OPT__CK_RESTRICT :
         Option : verify that the values stored in each cell at level "lv" are
         equal to the average values of its eight child cells at level "lv+1".

      OPT__CK_FLUX_ALLOCATE :
         Option : verify that the flux arrays are properly allocated for
         patches near the coarse-fine boundaries. This check should ALWAYS pass.

      OPT__CK_FINITE :
         Option : check whether there are cells with INF of NAN values.


-------------------------------------------
   b. Input__DumpTable :
-------------------------------------------

      This file provides the physical time for each data dump. The input values
      should be organized as following.

      Dump ID              Dump Time
      0           0.0909090920000000
      1           0.0965026241318910
      ...         ...
      ...         ...
      ***************END LINE***************

      The "Dump ID" recorded in the first column are useless and will not be
      loaded anyway. The second column records the physical time to dump data,
      which must be monotonically increasing. Note that the line
      "***************END LINE***************" must be provided in the end of
      the file. The program will stop loading data from this file once it has
      detected this line. Also note that the parameter "END_T" will be reset
      to the physical time of the last data dump. The original value of "END_T"
      in the parameter file will be overwritten.


-------------------------------------------
   c. Input__Flag_Rho :
   d. Input__Flag_RhoGradient :
   e. Input__Flag_PresGradient :
   f. Input__Flag_User :
-------------------------------------------

      These files provide the density, density gradient, pressure gradient,
      and user-defined thresholds for the domain refinement at each level.
      Below, we give an example of the data format in the file
      "Input__Flag_Rho".

      level                         Density
          0                             8.0
          1                            64.0
          2                           512.0
        ...                             ...
        ...                             ...

      The "level" recorded in the first column is useless and will not be
      loaded anyway. The second column records the threshold at each level.
      At least "NLEVEL-1" rows (excluding the header) must be provided. In the
      example above, patches at level 0 with density larger than 8.0 will be
      refined to level 1. Likewise, patches at level 2 with density larger than
      512.0 will be refined to level 3. Note that patches satisfying the
      refinement criteria may still not be refined if they do not satisfy the
      proper-nesting constraint.


-------------------------------------------
   g. Input__NoteScript
-------------------------------------------

      This file is used only to take notes for each simulation run. All text
      in this file will be copied to the file "Record__Note". The file
      "Record__Note" will also record all simulation parameters, symbolic
      constants, and GPU information in each run.



==========================================================================
III. Data Structure
==========================================================================

There are three main "structures" used by GAMER :

(1). AMR_t     : data structure for the AMR implementation
(2). Patch_t   : data structure of a single patch
(3). ParaVar_t : data structure for the parallelization

Please refer to the file "include/AMR.h" "include/Patch_t", and
"include/ParaVar_t" for detailed descriptions of these structures.


In the following, we clarify some definitions adopted in GAMER.

1. In order to perform the temporal interpolation for the individual time-step
   scheme, as well as to prevent the hydrodynamic solver from overwriting the
   previous data (which might be still useful for constructing the ghost-zone
   values of nearby patches), we allocate TWO data arrays (0 and 1) for each
   patch to store the hydrodynamic data (in the array "fluid") and the
   potential data (in the array "pot") at two different physical times.
   Conceptually, in the first step, the hydrodynamic solver will update the
   data stored in array 0 and store the updated results at array 1. In the next
   step, the hydrodynamic solver will update the data stored in array 1 and
   store the updated results at array 0. We use the variables "FluSg[NLEVEL]"
   and "PotSg[NLEVEL]" to record the array index (0 or 1) of the latest fluid
   and potential data at each level.
   ( "Sg" is a short name of "Sandglass" that implies the operation
     0->1->0->1->0->1 ... ).


2. Definition of the sibling index (SibID) :
   For each patch, the array "sibling[26]" is used to store the patch indices
   of sibling patches in 26 directions. They are ordered as following :

   Sibling  Direction (x,y,z)
         0            (-,0,0)
         1            (+,0,0)
         2            (0,-,0)
         3            (0,+,0)
         4            (0,0,-)
         5            (0,0,+)
         6            (-,-,0)
         7            (+,-,0)
         8            (-,+,0)
         9            (+,+,0)
        10            (0,-,-)
        11            (0,+,-)
        12            (0,-,+)
        13            (0,+,+)
        14            (-,0,-)
        15            (-,0,+)
        16            (+,0,-)
        17            (+,0,+)
        18            (-,-,-)
        19            (+,-,-)
        20            (-,+,-)
        21            (+,+,-)
        22            (-,-,+)
        23            (+,-,+)
        24            (-,+,+)
        25            (+,+,+)

   A schematic diagram would look like :

            24  13  25
            15  05  17     z+1 plane
            22  12  23

            08  03  09
            00  XX  01     z   plane
            06  02  07
   y
   ^        20  11  21
   |        14  04  16     z-1 plane
   --->x    18  10  19


3. Definition of the local index (LocalID) :
   For each patch, its local patch index in one patch group is defined as
   following :

   LocalID  Direction (x,y,z)
         0            (0,0,0)
         1            (1,0,0)
         2            (0,1,0)
         3            (0,0,1)
         4            (1,1,0)
         5            (0,1,1)
         6            (1,0,1)
         7            (1,1,1)

   A schematic diagram would look like :


            5  7       z+1 plane
   y        3  6
   ^
   |        2  4       z   plane
   --->x    0  1


==========================================================================
IV. Notice
==========================================================================
1. Currently only the periodic boundary condition is supported.

2. All global variables are defined in the file "src/GAMER/Main.cpp" and
   declared as "extern" in the file "include/Global.h".

3. Function prototypes are declared in the file "include/Prototype.h".

4. All symbolic constants (except for those defined in the Makefile) are
   defined in the files "include/Macro.h", "include/CUPOT.h", and
   "include/CUFLU.h". They should not be modified.

5. A simple structure "Timer_t" is defined in the file "include/Timer.h",
   which is used to measure the elapsed times of different parts in GAMER.

6. If one wants to add a new input parameter in the file "Input__Parameter",
   the following procedure should be proceeded.

   a. Add a new line in the file "Input_Parameter" with the following format.

      71.0           HUBBLE_PARAMETER          # DESCRIPTION

   b. Add the following lines in the correct position in the file
      "src/Init/Init_Load_Parameter.cpp". Since that all input parameters are
      loaded sequentially from the input file, these two lines must appear in
      the same order as they are put in the file "Input__Parameter".

      getline( &input_line, &len, File );
      sscanf( input_line, "%f%s", &HUBBLE_PARAMETER, string );

   c. Add the following line in the file "src/Auxiliary/Aux_TakeNote.cpp" to
      record the value of this new input.

      fprintf( Note, "HUBBLE_PARAMETER %13.7e\n", HUBBLE_PARAMETER );

   ***A more flexible input format will be provided in the future release***

7. The multi-level Poisson solver adopted in this version is different from
   the scheme described in the GAMER paper. In this version, when solving the
   Poisson equation, we let each patch to overlap between the neighboring
   patches by "POT_GHOST_SIZE" to make the potential smoothly across the patch
   interfaces. In other words, we solve the Poisson equation using

      "( PATCH_SIZE + 2*POT_GHOST_SIZE )^3"

   grids for each patch. We recommend setting "POT_GHOST_SIZE = 4 or 5" in
   order to achieve good accuracy.

8. To use Valgrind for debugging memory misusage, one can
   (1) Turn on the option "GAMER_DEBUG" in the makefile
   (2) Run GAMER by

         mpirun -np 8 valgrind -v --leak-check=full ./Dizzy 1>record 2>&1 &



==========================================================================
Revision History :
==========================================================================

Version     Date (D/M/Y)      Note
--------------------------------------------------------------------------
4.0         08/17/2010
--------------------------------------------------------------------------
4.0.1       09/01/2010
--------------------------------------------------------------------------
4.0.2       09/28/2010
--------------------------------------------------------------------------
4.1         12/21/2010
--------------------------------------------------------------------------
4.1.1       01/09/2011
--------------------------------------------------------------------------
4.1.2       02/07/2011        Latest version
--------------------------------------------------------------------------



==========================================================================
New Feature in Version 4.1.2
==========================================================================
1. Support Fermi GPUs

2. CUTIL library is no longer required

3. Support several hydro schemes, including
   a. RTVD/WAF/MHM/MHM_RP/CTU hydro schemes
   b. PLM/PPM data reconstruction in MHM/MHM_RP/CTU schemes
   c. EXACT/ROE/HLLE/HLLC Riemann solvers in WAF/MHM/MHM_RP/CTU
      schemes (WAF only supports EXACT/ROE solvers)

4. Support OpenMP parallelization in most of the CPU calculation
   --> overall performance is highly improved

5. Add the new parameter "BOX_SIZE", which set the physical size of the
   simulation box
   --> Now the spatial unit is determined by this parameter
   --> The physical grid size at each level is stored in "amr->dh[lv]", which
       is different from "amr->scale[lv]". The latter represents the grid
       size by the equivalent number of finest grids.

6. The output format version is updated to 1102

7. Support modifying the parameter "NLEVEL" during the restart process
   --> The NLEVEL in "Makefile" and "RESTART" files can be different
   --> In this case, the grid scale "amr->scale[lv] will be rescaled

8. Add the new parameter "OPT__OUTPUT_BASE"
   --> Only output the base-level data for the option "OPT__OUTPUT_PART"
   --> Data output with this option can be visualized by the tool
       "GAMER_2DBaseData2gnuplot.1.0"

9. Several visualization tools are put in "tool/Visualization"

10. For CUDA 4.0 + InfiniBand, it is recommended to set
    "export CUDA_NIC_INTEROP=1" to remove the incompatibility between CUDA and
    InfiniBand driver



