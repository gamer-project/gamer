This file records all the simulation configurations, including
* [[Simulation notes | Running-the-Code#taking-notes]]
* [[Configuration options | Installation:-Option-List#Option-List]]
* Symbolic constants defined in various headers
* [[Runtime parameters | Runtime Parameters]] set by various `Input__*` tables
* Simulation resolutions (i.e., cell width) on different AMR levels
* Compilation time
* OpenMP diagnosis (see also [[MPI Binding and Thread Affinity | MPI-and-OpenMP#mpi-binding-and-thread-affinity]])
* CPU and GPU specifications
* Total simulation time

Example:
``` markdown
Simulation Notes
***********************************************************************************
My notes
***********************************************************************************


Makefile Options (numerical schemes)
***********************************************************************************
MODEL                           HYDRO
GRAVITY                         ON
POT_SCHEME                      SOR
STORE_POT_GHOST                 ON
UNSPLIT_GRAVITY                 ON
COMOVING                        OFF
PARTICLE                        ON
SUPPORT_GRACKLE                 OFF
FLU_SCHEME                      CTU
LR_SCHEME                       PPM
RSOLVER                         ROE
DUAL_ENERGY                     DE_ENPY
STORE_PAR_ACC                   ON
STAR_FORMATION                  OFF
***********************************************************************************


Makefile Options (optimization and compilation)
***********************************************************************************
GPU                             ON
GAMER_DEBUG                     OFF
BITWISE_REPRODUCIBILITY         OFF
TIMING                          ON
TIMING_SOLVER                   OFF
FLOAT8                          OFF
SERIAL                          OFF
LOAD_BALANCE                    HILBERT
OVERLAP_MPI                     OFF
OPENMP                          ON
GPU_ARCH                        MAXWELL
LAOHU                           OFF
SUPPORT_HDF5                    ON
SUPPORT_GSL                     OFF
RANDOM_NUMBER                   RNG_GNU_EXT
***********************************************************************************


Other Options (in CUFLU.h and CUPOT.h)
***********************************************************************************
CHECK_NEGATIVE_IN_FLUID         OFF
CHAR_RECONSTRUCTION             OFF
CHECK_INTERMEDIATE              HLLE
HLL_NO_REF_STATE                OFF
HLL_INCLUDE_ALL_WAVES           OFF
WAF_DISSIPATE                   OFF
USE_PSOLVER_10TO14              ON
SOR_RHO_SHARED                  ON
SOR_CPOT_SHARED                 OFF
SOR_USE_SHUFFLE                 ON
SOR_USE_PADDING                 ON
SOR_MOD_REDUCTION               2
DT_FLU_USE_SHUFFLE              ON
DT_GRA_USE_SHUFFLE              ON
EXT_POT_NAUX_MAX                10
EXT_ACC_NAUX_MAX                10
***********************************************************************************


Symbolic Constants
***********************************************************************************
#define NCOMP_FLUID             5
#define NCOMP_PASSIVE           1
#define FLU_NIN                 6
#define FLU_NOUT                6
#define NFLUX_FLUID             5
#define NFLUX_PASSIVE           1
#define GRA_NIN                 5
#define PATCH_SIZE              8
#define MAX_PATCH               1000000
#define NLEVEL                  10

#define FLU_GHOST_SIZE          3
#define POT_GHOST_SIZE          5
#define RHO_GHOST_SIZE          4
#define GRA_GHOST_SIZE          2
#define USG_GHOST_SIZE          1
#define RHOEXT_GHOST_SIZE       2
#define FLU_NXT                 22
#define POT_NXT                 12
#define RHO_NXT                 16
#define GRA_NXT                 12
#define USG_NXT_F               18
#define USG_NXT_G               10
#define RHOEXT_NXT              12
#define FLU_BLOCK_SIZE_X        512
#define FLU_BLOCK_SIZE_Y        1
#define POT_BLOCK_SIZE_Z        4
#define GRA_BLOCK_SIZE_Z        4
#define DT_FLU_BLOCK_SIZE       512
#define DT_GRA_BLOCK_SIZE_ Z    4
#define PAR_NATT_TOTAL          11
#define PAR_NATT_USER           0
#define PAR_NATT_STORED         7
#define MAX_STRING              512
#define TINY_NUMBER             1.17549435082229e-38
#define HUGE_NUMBER             3.40282346638529e+38
***********************************************************************************


Parameters of Simulation Scale
***********************************************************************************
BOX_SIZE (input)                3.00000000000000e+00
BOX_SIZE_X                      3.00000000000000e+00
BOX_SIZE_Y                      3.00000000000000e+00
BOX_SIZE_Z                      3.00000000000000e+00
BOX_SCALE_X                     32768
BOX_SCALE_Y                     32768
BOX_SCALE_Z                     32768
NX0_TOT[0]                      64
NX0_TOT[1]                      64
NX0_TOT[2]                      64
MPI_NRank                       2
MPI_NRank_X[0]                  -1
MPI_NRank_X[1]                  -1
MPI_NRank_X[2]                  -1
OMP_NTHREAD                     6
END_T                           1.00000000000000e+00
END_STEP                        2147483647
***********************************************************************************


Parameters of Test Problems
***********************************************************************************
TESTPROB_ID                     11
***********************************************************************************


Parameters of Code Units
***********************************************************************************
OPT__UNIT                       0
***********************************************************************************


Parameters of Boundary Condition
***********************************************************************************
OPT__BC_FLU[0] (-x)             2
OPT__BC_FLU[1] (+x)             2
OPT__BC_FLU[2] (-y)             2
OPT__BC_FLU[3] (+y)             2
OPT__BC_FLU[4] (-z)             2
OPT__BC_FLU[5] (+z)             2
OPT__BC_POT                     2
GFUNC_COEFF0                    4.8000000e+00
***********************************************************************************


Parameters of Particle
***********************************************************************************
DEBUG_PARTICLE                  OFF
Par->NPar_Active_AllRank        20000
Par->Init                       1
Par->ParICFormat                1
Par->ParICMass                 -1.0000000e+00
Par->Interp                     3
Par->Integ                      2
Par->GhostSize                  1
Par->ImproveAcc                 1
Par->PredictPos                 1
Par->RemoveCell                 2.0000000e+00
***********************************************************************************


Parameters of Time-step Determination
***********************************************************************************
DT__FLUID                       5.0000000e-01
DT__FLUID_INIT                  5.0000000e-01
DT__GRAVITY                     5.0000000e-01
DT__PARVEL                      5.0000000e-01
DT__PARVEL_MAX                 -1.0000000e+00
DT__PARACC                      5.0000000e-01
DT__SYNC_PARENT_LV              1.0000000e-01
DT__SYNC_CHILDREN_LV            1.0000000e-01
OPT__DT_USER                    0
OPT__DT_LEVEL                   3
AUTO_REDUCE_DT                  1
AUTO_REDUCE_DT_FACTOR           8.0000000e-01
AUTO_REDUCE_DT_FACTOR_MIN       1.0000000e-01
OPT__RECORD_DT                  1
***********************************************************************************


Parameters of Domain Refinement
***********************************************************************************
REGRID_COUNT                    4
FLAG_BUFFER_SIZE                8
FLAG_BUFFER_SIZE_MAXM1_LV       4
FLAG_BUFFER_SIZE_MAXM2_LV       6
MAX_LEVEL                       2
OPT__FLAG_RHO                   0
OPT__FLAG_RHO_GRADIENT          0
OPT__FLAG_PRES_GRADIENT         0
OPT__FLAG_VORTICITY             0
OPT__FLAG_JEANS                 0
OPT__FLAG_LOHNER_DENS           0
OPT__FLAG_LOHNER_ENGY           0
OPT__FLAG_LOHNER_PRES           0
OPT__FLAG_LOHNER_TEMP           0
OPT__FLAG_LOHNER_FORM           LOHNER_FLASH2
OPT__FLAG_USER                  0
OPT__FLAG_REGION                0
OPT__FLAG_NPAR_PATCH            2
OPT__FLAG_NPAR_CELL             0
OPT__FLAG_PAR_MASS_CELL         0
OPT__NO_FLAG_NEAR_BOUNDARY      0
OPT__PATCH_COUNT                1
OPT__PARTICLE_COUNT             1
OPT__REUSE_MEMORY               2
OPT__MEMORY_POOL                0
***********************************************************************************


Parameters of Parallelization
***********************************************************************************
Flu_ParaBuf                     3
Pot_ParaBuf                     4
Rho_ParaBuf                     4
LB_WLI_MAX                      1.0000000e-01
LB_PAR_WEIGHT                   0.0000000e+00
OPT__RECORD_LOAD_BALANCE        1
OPT__MINIMIZE_MPI_BARRIER       1
***********************************************************************************


Parameters of Fluid Solver (in different models)
***********************************************************************************
GAMMA                           1.6666667e+00
MOLECULAR_WEIGHT                6.0000000e-01
MINMOD_COEFF                    1.5000000e+00
EP_COEFF                        1.2500000e+00
OPT__LR_LIMITER                 VL_GMINMOD
OPT__WAF_LIMITER                NONE
OPT__1ST_FLUX_CORR              3D1D
OPT__1ST_FLUX_CORR_SCHEME       RSOLVER_1ST_ROE
***********************************************************************************


Parameters of Fluid Solver (in all models)
***********************************************************************************
FLU_GPU_NPGROUP                 768
GPU_NSTREAM                     32
OPT__FIXUP_FLUX                 1
OPT__FIXUP_RESTRICT             1
OPT__CORR_AFTER_ALL_SYNC        2
OPT__NORMALIZE_PASSIVE          0
   Number of fields             0
OPT__OVERLAP_MPI                0
OPT__RESET_FLUID                0
MIN_DENS                        0.0000000e+00
MIN_PRES                        1.0000000e-15
JEANS_MIN_PRES                  0
DUAL_ENERGY_SWITCH              2.0000000e-02
WITH_COARSE_FINE_FLUX           1
MPI Thread Level                MPI_THREAD_SINGLE
***********************************************************************************


Parameters of Poisson and Gravity Solvers
***********************************************************************************
NEWTON_G                        1.0000000e+00
SOR_OMEGA                       1.6900000e+00
SOR_MAX_ITER                    60
SOR_MIN_ITER                    10
POT_GPU_NPGROUP                 768
OPT__GRA_P5_GRADIENT            0
OPT__SELF_GRAVITY               1
OPT__EXT_ACC                    0
OPT__EXT_POT                    0
AveDensity_Init                 2.9407660e-04
***********************************************************************************


Parameters of Initialization
***********************************************************************************
OPT__INIT                       1
RESTART_LOAD_NRANK              1
OPT__RESTART_RESET              0
OPT__UM_IC_LEVEL                0
OPT__UM_IC_NVAR                 -1
OPT__UM_IC_FORMAT               1
OPT__UM_IC_DOWNGRADE            1
OPT__UM_IC_REFINE               1
OPT__UM_IC_LOAD_NRANK           1
OPT__INIT_RESTRICT              1
OPT__INIT_GRID_WITH_OMP         1
OPT__GPUID_SELECT               -1
INIT_SUBSAMPLING_NCELL          0
***********************************************************************************


Parameters of Interpolation Schemes
***********************************************************************************
OPT__INT_TIME                   1
OPT__FLU_INT_SCHEME             CQUAD
OPT__POT_INT_SCHEME             QUAD
OPT__RHO_INT_SCHEME             CQUAD
OPT__GRA_INT_SCHEME             QUAD
OPT__REF_FLU_INT_SCHEME         CQUAD
OPT__REF_POT_INT_SCHEME         QUAD
INT_MONO_COEFF                  2.0000000e+00
***********************************************************************************


Parameters of Data Dump
***********************************************************************************
OPT__OUTPUT_TOTAL               1
OPT__OUTPUT_PART                0
OPT__OUTPUT_USER                0
OPT__OUTPUT_PAR_TEXT            1
OPT__OUTPUT_BASEPS              0
OPT__OUTPUT_BASE                0
OPT__OUTPUT_POT                 0
OPT__OUTPUT_PAR_DENS            1
OPT__OUTPUT_MODE                2
OUTPUT_STEP                     5
OUTPUT_DT                       5.00000000000000e-01
OUTPUT_PART_X                   -1.00000000000000e+00
OUTPUT_PART_Y                   -1.00000000000000e+00
OUTPUT_PART_Z                   -1.00000000000000e+00
INIT_DUMPID                     -1
***********************************************************************************


Parameters of Miscellaneous Purposes
***********************************************************************************
OPT__VERBOSE                    0
OPT__TIMING_BARRIER             0
OPT__TIMING_BALANCE             0
OPT__TIMING_MPI                 0
OPT__RECORD_MEMORY              1
OPT__RECORD_PERFORMANCE         1
OPT__MANUAL_CONTROL             1
OPT__RECORD_USER                0
OPT__OPTIMIZE_AGGRESSIVE        0
***********************************************************************************


Parameters of Simulation Checks
***********************************************************************************
OPT__CK_REFINE                  0
OPT__CK_PROPER_NESTING          0
OPT__CK_CONSERVATION            1
OPT__CK_NORMALIZE_PASSIVE       0
OPT__CK_RESTRICT                0
OPT__CK_FINITE                  0
OPT__CK_PATCH_ALLOCATE          0
OPT__CK_FLUX_ALLOCATE           0
OPT__CK_NEGATIVE                0
OPT__CK_MEMFREE                 1.0000000e+00
OPT__CK_PARTICLE                0
***********************************************************************************


Flag Criterion (# of Particles per Patch)
***********************************************************************************
  Level      # of Particles
      0                3000
      1                3000
***********************************************************************************


Cell Size and Scale (scale = number of cells at the finest level)
***********************************************************************************
  Level                           Cell Size                Cell Scale
      0              0.04687500000000000000                       512
      1              0.02343750000000000000                       256
      2              0.01171875000000000000                       128
      3              0.00585937500000000000                        64
      4              0.00292968750000000000                        32
      5              0.00146484375000000000                        16
      6              0.00073242187500000000                         8
      7              0.00036621093750000000                         4
      8              0.00018310546875000000                         2
      9              0.00009155273437500000                         1
***********************************************************************************


Compilation Time
***********************************************************************************
Sep 21 2018 21:04:55
***********************************************************************************


OpenMP Diagnosis
***********************************************************************************
OMP__SCHEDULE                   DYNAMIC
OMP__SCHEDULE_CHUNK_SIZE        1
OMP__NESTED                     OFF

CPU core IDs of all OpenMP threads (tid == thread ID):
------------------------------------------------------------------------
 Rank        Host  NThread  tid-00  tid-01  tid-02  tid-03  tid-04  tid-05
    0      hulk32        6       0       5       2       1       1       3
    1      hulk32        6       0       4       3       4       4       4
***********************************************************************************


Device Diagnosis
***********************************************************************************
MPI_Rank =   0, hostname =     hulk32, PID = 24703

CPU Info :
CPU Type        : Intel(R) Core(TM) i7-3930K CPU @ 3.20GHz
CPU MHz         : 3201.000
Cache Size      : 12288 KB
CPU Cores       : 6
Total Memory    : 63.0 GB

GPU Info :
Number of GPUs                    : 1
GPU ID                            : 0
GPU Name                          : GeForce GTX TITAN X
CUDA Driver Version               : 7.50
CUDA Runtime Version              : 7.50
CUDA Major Revision Number        : 5
CUDA Minor Revision Number        : 2
Clock Rate                        : 1.076000 GHz
Global Memory Size                : 12287 MB
Constant Memory Size              : 64 KB
Shared Memory Size per Block      : 48 KB
Number of Registers per Block     : 65536
Warp Size                         : 32
Number of Multiprocessors:        : 24
Number of Cores per Multiprocessor: 128
Total Number of Cores:            : 3072
Max Number of Threads per Block   : 1024
Max Size of the Block X-Dimension : 1024
Max Size of the Grid X-Dimension  : 2147483647
Concurrent Copy and Execution     : Yes
Concurrent Up/Downstream Copies   : Yes
Concurrent Kernel Execution       : Yes
GPU has ECC Support Enabled       : No


MPI_Rank =   1, hostname =     hulk32, PID = 24704

CPU Info :
CPU Type        : Intel(R) Core(TM) i7-3930K CPU @ 3.20GHz
CPU MHz         : 3201.000
Cache Size      : 12288 KB
CPU Cores       : 6
Total Memory    : 63.0 GB

GPU Info :
Number of GPUs                    : 1
GPU ID                            : 0
GPU Name                          : GeForce GTX TITAN X
CUDA Driver Version               : 7.50
CUDA Runtime Version              : 7.50
CUDA Major Revision Number        : 5
CUDA Minor Revision Number        : 2
Clock Rate                        : 1.076000 GHz
Global Memory Size                : 12287 MB
Constant Memory Size              : 64 KB
Shared Memory Size per Block      : 48 KB
Number of Registers per Block     : 65536
Warp Size                         : 32
Number of Multiprocessors:        : 24
Number of Cores per Multiprocessor: 128
Total Number of Cores:            : 3072
Max Number of Threads per Block   : 1024
Max Size of the Block X-Dimension : 1024
Max Size of the Grid X-Dimension  : 2147483647
Concurrent Copy and Execution     : Yes
Concurrent Up/Downstream Copies   : Yes
Concurrent Kernel Execution       : Yes
GPU has ECC Support Enabled       : No
***********************************************************************************

Total Processing Time : 21.777708 s
```

<br>

## Links
* [[Simulation Logs]]
