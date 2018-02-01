


# ===========================================================================
# INSTRUCTIONS
#
# To compile GAMER, please set the following configurations properly:
#
# (1) simulation options
# (2) compiler and flags
# (3) library paths
# ===========================================================================



# executable
#######################################################################################################
EXECUTABLE := gamer



# output detailed compilation commands (0/1 = off/on)
#######################################################################################################
COMPILE_VERBOSE := 0



# simulation options
#######################################################################################################

# (a) physical modules
# -------------------------------------------------------------------------------
# physical model: HYDRO/MHD/ELBDM/PAR_ONLY (ELBDM: wave dark matter; PAR_ONLY: particle-only)
# --> must be set in any cases; PAR_ONLY is not supported yet;
#     MHD and ELBDM are not publicly available yet but will be released soon
SIMU_OPTION += -DMODEL=HYDRO

# enable gravity
#SIMU_OPTION += -DGRAVITY

# enable particles
# --> must enable GRAVITY
#SIMU_OPTION += -DPARTICLE

# support Grackle, a chemistry and radiative cooling library
# --> please enable OpenMP when compiling Grackle (by "make omp-on")
#SIMU_OPTION += -DSUPPORT_GRACKLE


# (b) options of different physical modules
# (b-1) hydro options
# ------------------------------------------------------------------------------------
ifeq "$(filter -DMODEL=HYDRO, $(SIMU_OPTION))" "-DMODEL=HYDRO"
# hydrodynamic scheme: RTVD/WAF/MHM/MHM_RP/CTU
# --> must be set when MODEL=HYDRO; WAF is deprecated
SIMU_OPTION += -DFLU_SCHEME=CTU

# scheme of spatial data reconstruction: PLM/PPM (piecewise-linear/piecewise-parabolic)
# --> useless for RTVD/WAF
SIMU_OPTION += -DLR_SCHEME=PPM

# Riemann solver: EXACT/ROE/HLLE/HLLC
# --> all are useless for RTVD; HLLE/HLLC are useless for WAF
SIMU_OPTION += -DRSOLVER=ROE

# dual energy formalism: DE_ENPY/DE_EINT (evolve entropy or internal energy)
# --> DE_EINT is not supported yet; useless for RTVD/WAF
#SIMU_OPTION += -DDUAL_ENERGY=DE_ENPY

# number of passively advected scalars (**excluding** entropy or internal energy used by DUAL_ENERGY)
# --> useless for RTVD/WAF
# --> for simulations without any passive scalar, please set this variable to 0 and do NOT comment it out
SIMU_OPTION += -DNCOMP_PASSIVE_USER=0


# (b-2) MHD options
# ------------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=MHD, $(SIMU_OPTION))" "-DMODEL=MHD"


# (b-3) ELBDM options
# ------------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=ELBDM, $(SIMU_OPTION))" "-DMODEL=ELBDM"
# enforce the mass conservation
SIMU_OPTION += -DCONSERVE_MASS

# 4-th order Laplacian
SIMU_OPTION += -DLAPLACIAN_4TH

# include the quartic self-interaction potential
# --> must turn on GRAVITY; does not support COMOVING
#SIMU_OPTION += -DQUARTIC_SELF_INTERACTION

# number of passively advected scalars
# --> not supported yet and must be set to 0
SIMU_OPTION += -DNCOMP_PASSIVE_USER=0

endif # MODEL


# (b-4) gravity options
# -------------------------------------------------------------------------------
ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
# Poisson solver: SOR/MG (successive-overrelaxation (recommended)/multigrid)
# --> must be set when GRAVITY is enabled
SIMU_OPTION += -DPOT_SCHEME=SOR

# store GRA_GHOST_SIZE ghost-zone potential for each patch on each side
# --> recommended when PARTICLE is enabled for improving accuracy for particles around the patch boundaries
# --> must be enabled for STAR_FORMATION + STORE_PAR_ACC
SIMU_OPTION += -DSTORE_POT_GHOST

# use unsplitting method to couple gravity to the target model (recommended)
# --> for HYDRO and MHD only; ELBDM is not supported yet
SIMU_OPTION += -DUNSPLIT_GRAVITY

# comoving frame for cosmological simulations
# --> not supported for PARTICLE
#SIMU_OPTION += -DCOMOVING
endif # GRAVITY


# (b-5) particle options
# ------------------------------------------------------------------------------------
ifeq "$(filter -DPARTICLE, $(SIMU_OPTION))" "-DPARTICLE"
# store particle acceleration (recommended)
SIMU_OPTION += -DSTORE_PAR_ACC

# allow for creating new particles after initialization
# --> must turn on STORE_POT_GHOST when STORE_PAR_ACC is adopted
#SIMU_OPTION += -DSTAR_FORMATION

# number of passive particle attributes (**excluding** PAR_CREATION_TIME used by STAR_FORMATION)
# --> for simulations without any passive attributes, please set this variable to 0 and do NOT comment it out
SIMU_OPTION += -DPAR_NPASSIVE_USER=0
endif # PARTICLE


# (d) miscellaneous options (AMR, GPU, parallelization, optimizations, ...)
# ------------------------------------------------------------------------------------
# maximum number of AMR levels including the root level (level id = 0 ... NLEVEL-1)
# --> must be set in any cases
SIMU_OPTION += -DNLEVEL=10

# maximum number of patches on each AMR level
# --> must be set in any cases
SIMU_OPTION += -DMAX_PATCH=1000000

# GPU acceleration
# --> must set GPU_ARCH as well
#SIMU_OPTION += -DGPU

# GPU architecture: FERMI/KEPLER/MAXWELL/PASCAL
SIMU_OPTION += -DGPU_ARCH=KEPLER

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG

# bitwise reproducibility
#SIMU_OPTION += -DBITWISE_REPRODUCIBILITY

# measure the wall-clock time of GAMER
SIMU_OPTION += -DTIMING

# measure the wall-clock time of GPU solvers (will disable CPU/GPU overlapping)
# --> must enable TIMING
#SIMU_OPTION += -DTIMING_SOLVER

# double precision
#SIMU_OPTION += -DFLOAT8

# serial mode (in which no MPI libraries are required)
# --> must disable LOAD_BALANCE
SIMU_OPTION += -DSERIAL

# load-balance parallelization: HILBERT
# --> must disable SERIAL
#SIMU_OPTION += -DLOAD_BALANCE=HILBERT

# overlap MPI communication with computation
# --> NOT supported yet; must enable LOAD_BALANCE
#SIMU_OPTION += -DOVERLAP_MPI

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP

# work on the NAOC Laohu GPU cluster
#SIMU_OPTION += -DLAOHU

# support HDF5 format
#SIMU_OPTION += -DSUPPORT_HDF5

# support GNU scientific library
#SIMU_OPTION += -DSUPPORT_GSL

# support yt inline analysis
#SIMU_OPTION += -DSUPPORT_LIBYT





# compiler and flags
#######################################################################################################

# intel
# -------------------------------------------------------------------------------
# CXX         = icpc                                      # serial compiler
##CXX         = $(MPI_PATH)/bin/mpicxx                    # MPI compiler
##CXX         = CC                                        # cray wrapper script
# CXXFLAG     = -g -O3                                    # optimization flags
##CXXFLAG     = -g -fast
# CXXFLAG    += -w1                                       # warning flags
# OPENMPFLAG  = -fopenmp                                  # openmp flag
# LIB         = -limf                                     # libraries and linker flags
#
## for debug only
#ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
##CXXFLAG    += -fstack-protector-all
##LIB        += -lssp
#endif
#
## suppress warning when OpenMP is disabled
#ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
# CXXFLAG    += -Wno-unknown-pragmas -diag-disable 3180
#endif


# gnu
# -------------------------------------------------------------------------------
 CXX         = g++                                       # serial compiler
#CXX         = $(MPI_PATH)/bin/mpicxx                    # MPI compiler
#CXX         = CC                                        # cray wrapper script
 CXXFLAG     = -g -O3                                    # optimization flags
#CXXFLAG     = -g -Ofast
 CXXFLAG    += -Wall -Wextra                             # warning flags
 CXXFLAG    += -Wno-unused-variable -Wno-unused-parameter \
               -Wno-maybe-uninitialized -Wno-unused-but-set-variable \
               -Wno-unused-result -Wno-unused-function
 OPENMPFLAG  = -fopenmp                                  # openmp flag
 LIB         =                                           # libraries and linker flags

# for debug only
ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
#CXXFLAG    += -fstack-protector-all
endif

# suppress warning when OpenMP is disabled
ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
 CXXFLAG    += -Wno-unknown-pragmas
endif


# CUDA
# -------------------------------------------------------------------------------
 NVCC = $(CUDA_PATH)/bin/nvcc                            # CUDA compiler
#NVCC = nvcc -ccbin CC                                   # CUDA compiler used in the Cray environment





# library paths
#######################################################################################################

# template
CUDA_PATH    :=
FFTW_PATH    :=
MPI_PATH     :=
HDF5_PATH    :=
GRACKLE_PATH :=
GSL_PATH     :=

# NCSA-bluewaters (only need to specify the path of user-installed packages)
#GRACKLE_PATH :=

# UIUC-campus (intel-15.0 + openmpi-1.8.4)
#CUDA_PATH    := /usr/local/cuda/7.0
#FFTW_PATH    := /projects/ncsa/grav/softwares/fftw/2.1.5
#MPI_PATH     := /usr/local/mpi/openmpi-1.8.4-intel-15.0
#HDF5_PATH    := /projects/ncsa/grav/softwares/miniconda2
#GRACKLE_PATH := /home/hyschive/local
#GSL_PATH     := /usr/local/gsl/gsl-1.16

# UIUC-campus (gcc-4.9.2 + openmpi-1.8.4)
#CUDA_PATH    := /usr/local/cuda/7.0
#FFTW_PATH    := /projects/ncsa/grav/softwares/fftw/2.1.5-openmpi-1.8.4-gcc-4.9.2
#MPI_PATH     := /usr/local/mpi/openmpi-1.8.4-gcc-4.9.2
#HDF5_PATH    := /projects/ncsa/grav/softwares/miniconda2
#GRACKLE_PATH := /home/hyschive/local
#GSL_PATH     := /usr/local/gsl/gsl-1.16

# NCSA-hybrid (gcc)
#CUDA_PATH    := /usr/local/cuda-8.0
#MPI_PATH     := /home/hyschive/project/software/openmpi/2.0.2
#FFTW_PATH    := /home/hyschive/project/software/fftw/2.1.5
#HDF5_PATH    := /home/hyschive/project/software/miniconda2
#GSL_PATH     :=

# NCSA-nano (intel 2016)
#CUDA_PATH    := /usr/local/cuda
#MPI_PATH     := /mnt/lustre/allen/hyschive/software/openmpi/2.0.2
#FFTW_PATH    := /mnt/lustre/allen/hyschive/software/fftw/2.1.5
#HDF5_PATH    := /mnt/lustre/allen/hyschive/software/miniconda2
#GSL_PATH     :=

# NCSA-dxl (gcc, serial)
#CUDA_PATH    :=
#FFTW_PATH    := /home/hyschive/software/fftw/2.1.5
#MPI_PATH     := /usr
#HDF5_PATH    := /home/hyschive/software/miniconda2
#GSL_PATH     := /usr/lib
#LIBYT_PATH   := /home/hyschive/project/libyt/libyt.dev/libyt

# NTU-hulk (openmpi-intel-qlc)
#CUDA_PATH    := /opt/gpu/cuda/4.2
#FFTW_PATH    := /opt/math/fftw/2.1.5-intel-qlc
#MPI_PATH     := /usr/mpi/intel/openmpi-1.4.3-qlc
#HDF5_PATH    := /opt/hdf5/1.8.9
#MPI_PATH     := /opt/mpi/openmpi/1.6.0-intel
#CUDA_PATH    := /usr/local/cuda-7.5

# NTU-hulk (openmpi-gcc-qlc)
#CUDA_PATH    := /opt/gpu/cuda/4.2
#FFTW_PATH    := /opt/math/fftw/2.1.5-gcc-qlc
#MPI_PATH     := /opt/mpi/openmpi/1.8.1-gcc-qlc
#HDF5_PATH    := /opt/hdf5/1.8.9

# NAOC-laohu (OpenMPI-1.3.2-intel)
#CUDA_PATH    := $(CUDA_INSTALL_PATH)
#FFTW_PATH    := /ifs/data/hsiyu/Software/fftw/2.1.5-openmpi-1.3.2-intel
#MPI_PATH     := /usr/mpi/intel/openmpi-1.3.2-qlc
#GPUID_PATH   := /ifs/data/hsiyu/Software/gpudevmgr





# source files
#######################################################################################################

# common source files
# -------------------------------------------------------------------------------
# Cuda source files (compiled with nvcc)
CUDA_FILE   := CUAPI_Asyn_FluidSolver.cu  CUAPI_DiagnoseDevice.cu  CUAPI_MemAllocate_Fluid.cu \
               CUAPI_MemFree_Fluid.cu  CUAPI_Set_Default_GPU_Parameter.cu  CUAPI_SetDevice.cu \
               CUAPI_Synchronize.cu  CUAPI_Asyn_dtSolver.cu


# C/C++ source files (compiled with c++ compiler)
CC_FILE     := Main.cpp  EvolveLevel.cpp  InvokeSolver.cpp  Prepare_PatchData.cpp \
               InterpolateGhostZone.cpp

CC_FILE     += Aux_Check_Parameter.cpp  Aux_Check_Conservation.cpp  Aux_Check.cpp  Aux_Check_Finite.cpp \
               Aux_Check_FluxAllocate.cpp  Aux_Check_PatchAllocate.cpp  Aux_Check_ProperNesting.cpp \
               Aux_Check_Refinement.cpp  Aux_Check_Restrict.cpp  Aux_Error.cpp  Aux_GetCPUInfo.cpp \
               Aux_GetMemInfo.cpp  Aux_Message.cpp  Aux_Record_PatchCount.cpp  Aux_TakeNote.cpp  Aux_Timing.cpp \
               Aux_Check_MemFree.cpp  Aux_Record_Performance.cpp  Aux_CheckFileExist.cpp  Aux_Array.cpp \
               Aux_Record_User.cpp  Aux_Record_CorrUnphy.cpp  Aux_SwapPointer.cpp  Aux_Check_NormalizePassive.cpp \
               Aux_LoadTable.cpp

CC_FILE     += CPU_FluidSolver.cpp  Flu_AdvanceDt.cpp  Flu_Prepare.cpp  Flu_Close.cpp  Flu_FixUp.cpp \
               Flu_Restrict.cpp  Flu_AllocateFluxArray.cpp  Flu_BoundaryCondition_User.cpp  Flu_ResetByUser.cpp \
               Flu_CorrAfterAllSync.cpp  Flu_SwapFluxPointer.cpp

CC_FILE     += End_GAMER.cpp  End_MemFree.cpp  End_MemFree_Fluid.cpp  End_StopManually.cpp  End_User.cpp \
               Init_BaseLevel.cpp  Init_GAMER.cpp  Init_Load_DumpTable.cpp \
               Init_Load_FlagCriteria.cpp  Init_Load_Parameter.cpp  Init_MemAllocate.cpp \
               Init_MemAllocate_Fluid.cpp  Init_Parallelization.cpp  Init_RecordBasePatch.cpp  Init_Refine.cpp \
               Init_ByRestart_v1.cpp  Init_ByFunction.cpp  Init_TestProb.cpp  Init_ByFile.cpp  Init_OpenMP.cpp \
               Init_ByRestart_HDF5.cpp  Init_ResetParameter.cpp  Init_ByRestart_v2.cpp  Init_MemoryPool.cpp \
               Init_Unit.cpp  Init_PassiveVariable.cpp

CC_FILE     += Interpolate.cpp  Int_CQuadratic.cpp  Int_MinMod1D.cpp  Int_MinMod3D.cpp  Int_vanLeer.cpp \
               Int_Quadratic.cpp  Int_Table.cpp  Int_CQuartic.cpp  Int_Quartic.cpp

CC_FILE     += Mis_CompareRealValue.cpp  Mis_GetTotalPatchNumber.cpp  Mis_GetTimeStep.cpp  Mis_Heapsort.cpp \
               Mis_BinarySearch.cpp  Mis_1D3DIdx.cpp  Mis_Matching.cpp  Mis_GetTimeStep_User.cpp \
               Mis_dTime2dt.cpp  Mis_CoordinateTransform.cpp  Mis_BinarySearch_Real.cpp  Mis_InterpolateFromTable.cpp \
               CPU_dtSolver.cpp  dt_Prepare_Flu.cpp  dt_Prepare_Pot.cpp  dt_Close.cpp  dt_InvokeSolver.cpp

CC_FILE     += Output_DumpData_Total.cpp  Output_DumpData.cpp  Output_DumpManually.cpp  Output_PatchMap.cpp \
               Output_DumpData_Part.cpp  Output_FlagMap.cpp  Output_Patch.cpp  Output_PreparedPatch_Fluid.cpp \
               Output_PatchCorner.cpp  Output_Flux.cpp  Output_User.cpp  Output_BasePowerSpectrum.cpp \
               Output_DumpData_Total_HDF5.cpp  Output_L1Error.cpp

CC_FILE     += Flag_Real.cpp  Refine.cpp   SiblingSearch.cpp  SiblingSearch_Base.cpp  FindFather.cpp \
               Flag_User.cpp  Flag_Check.cpp  Flag_Lohner.cpp  Flag_Region.cpp

CC_FILE     += Table_01.cpp  Table_02.cpp  Table_03.cpp  Table_04.cpp  Table_05.cpp  Table_06.cpp \
               Table_07.cpp

vpath %.cu     GPU_API
vpath %.cpp    GAMER  Init  Refine  Fluid  Interpolation  Tables  Output  Miscellaneous  Auxiliary


# hydrodynamic source files (included only if "MODEL=HYDRO")
# ------------------------------------------------------------------------------------
ifeq "$(filter -DMODEL=HYDRO, $(SIMU_OPTION))" "-DMODEL=HYDRO"
CUDA_FILE   += CUFLU_dtSolver_HydroCFL.cu  CUFLU_FluidSolver_RTVD.cu  CUFLU_FluidSolver_WAF.cu  CUFLU_FluidSolver_MHM.cu \
               CUFLU_FluidSolver_CTU.cu

CC_FILE     += CPU_FluidSolver_RTVD.cpp  CPU_FluidSolver_WAF.cpp  CPU_FluidSolver_MHM.cpp  CPU_FluidSolver_CTU.cpp \
               CPU_Shared_DataReconstruction.cpp  CPU_Shared_FluUtility.cpp  CPU_Shared_ComputeFlux.cpp \
               CPU_Shared_FullStepUpdate.cpp  CPU_Shared_RiemannSolver_Exact.cpp  CPU_Shared_RiemannSolver_Roe.cpp \
               CPU_Shared_RiemannSolver_HLLE.cpp  CPU_Shared_RiemannSolver_HLLC.cpp  CPU_Shared_DualEnergy.cpp \
               CPU_dtSolver_HydroCFL.cpp

CC_FILE     += Hydro_Init_ByFunction_AssignData.cpp  Hydro_Aux_Check_Negative.cpp \
               Hydro_Init_ByFile_AssignData.cpp  Hydro_BoundaryCondition_Outflow.cpp  Hydro_BoundaryCondition_Reflecting.cpp \
               Hydro_Flag_Vorticity.cpp

vpath %.cu     Model_Hydro/GPU_Hydro
vpath %.cpp    Model_Hydro/CPU_Hydro  Model_Hydro

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
CUDA_FILE   += CUPOT_HydroGravitySolver.cu  CUPOT_dtSolver_HydroGravity.cu

CC_FILE     += CPU_HydroGravitySolver.cpp  CPU_dtSolver_HydroGravity.cpp

vpath %.cu     Model_Hydro/GPU_HydroGravity
vpath %.cpp    Model_Hydro/CPU_HydroGravity
endif


# MHD source files (included only if "MODEL=MHD")
# -------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=MHD, $(SIMU_OPTION))" "-DMODEL=MHD"

vpath %.cu     Model_MHD/GPU_MHD
vpath %.cpp    Model_MHD/CPU_MHD  Model_MHD

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"

vpath %.cu     Model_MHD/GPU_MHDGravity
vpath %.cpp    Model_MHD/CPU_MHDGravity
endif


# ELBDM source files (included only inf "MODEL=ELBDM")
# -------------------------------------------------------------------------------
else ifeq "$(filter -DMODEL=ELBDM, $(SIMU_OPTION))" "-DMODEL=ELBDM"
CUDA_FILE   += CUFLU_ELBDMSolver.cu

CC_FILE     += CPU_ELBDMSolver.cpp  ELBDM_Init_ByFunction_AssignData.cpp  ELBDM_Init_ByFile_AssignData.cpp \
               ELBDM_GetTimeStep_Fluid.cpp  ELBDM_Flag_EngyDensity.cpp  ELBDM_UnwrapPhase.cpp \
               ELBDM_GetTimeStep_Phase.cpp  ELBDM_SetTaylor3Coeff.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDM
vpath %.cpp    Model_ELBDM/CPU_ELBDM  Model_ELBDM

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
CUDA_FILE   += CUPOT_ELBDMGravitySolver.cu

CC_FILE     += CPU_ELBDMGravitySolver.cpp  ELBDM_GetTimeStep_Gravity.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDMGravity
vpath %.cpp    Model_ELBDM/CPU_ELBDMGravity
endif

endif # MODEL


# self-gravity source files (included only if "GRAVITY" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
CUDA_FILE   += CUAPI_MemAllocate_PoissonGravity.cu  CUAPI_MemFree_PoissonGravity.cu \
               CUAPI_Asyn_PoissonGravitySolver.cu  CUAPI_Init_ExternalAccPot.cu

CUDA_FILE   += CUPOT_PoissonSolver_SOR_10to14cube.cu  CUPOT_PoissonSolver_SOR_16to18cube.cu \
               CUPOT_PoissonSolver_MG.cu

CC_FILE     += CPU_PoissonGravitySolver.cpp  CPU_PoissonSolver_SOR.cpp  CPU_PoissonSolver_FFT.cpp \
               CPU_PoissonSolver_MG.cpp

CC_FILE     += Init_FFTW.cpp  Gra_Close.cpp  Gra_Prepare_Flu.cpp  Gra_Prepare_Pot.cpp  Gra_Prepare_Corner.cpp \
               Gra_AdvanceDt.cpp  Poi_Close.cpp  Poi_Prepare_Pot.cpp  Poi_Prepare_Rho.cpp \
               Output_PreparedPatch_Poisson.cpp  Init_MemAllocate_PoissonGravity.cpp  CPU_ExternalPot.cpp \
               End_MemFree_PoissonGravity.cpp  Init_Set_Default_SOR_Parameter.cpp \
               Init_Set_Default_MG_Parameter.cpp  Poi_GetAverageDensity.cpp  Init_GreenFuncK.cpp \
               Init_ExternalPot.cpp  Poi_BoundaryCondition_Extrapolation.cpp  CPU_ExternalAcc.cpp \
               Gra_Prepare_USG.cpp  Init_ExternalAcc.cpp  Poi_StorePotWithGhostZone.cpp  Init_ExternalAccPot.cpp

vpath %.cu     SelfGravity/GPU_Poisson  SelfGravity/GPU_Gravity
vpath %.cpp    SelfGravity/CPU_Poisson  SelfGravity/CPU_Gravity  SelfGravity
endif # GRAVITY


# particle source files (included only if "PARTICLE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DPARTICLE, $(SIMU_OPTION))" "-DPARTICLE"
CUDA_FILE   +=

CC_FILE     += Par_Init_ByFunction.cpp  Par_Output_TextFile.cpp  Par_FindHomePatch_Base.cpp  Par_PassParticle2Son.cpp \
               Par_Aux_Check_Particle.cpp  Par_PassParticle2Father.cpp  Par_CollectParticle2OneLevel.cpp \
               Par_MassAssignment.cpp  Par_UpdateParticle.cpp  Par_GetTimeStep_VelAcc.cpp \
               Par_PassParticle2Sibling.cpp  Par_CountParticleInDescendant.cpp  Par_Aux_GetConservedQuantity.cpp \
               Par_Aux_InitCheck.cpp  Par_Aux_Record_ParticleCount.cpp  Par_PassParticle2Son_AllPatch.cpp \
               Par_Synchronize.cpp  Par_PredictPos.cpp  Par_Init_ByFile.cpp

vpath %.cu     Particle/GPU
vpath %.cpp    Particle/CPU  Particle

ifeq "$(findstring -DLOAD_BALANCE, $(SIMU_OPTION))" "-DLOAD_BALANCE"
CC_FILE     += Par_LB_Init_RedistributeByRectangular.cpp  Par_LB_SendParticleData.cpp  Par_LB_CollectParticle2OneLevel.cpp \
               Par_LB_CollectParticleFromRealPatch.cpp  Par_LB_RecordExchangeParticlePatchID.cpp \
               Par_LB_MapBuffer2RealPatch.cpp  Par_LB_ExchangeParticleBetweenPatch.cpp \
               Par_LB_Refine_SendParticle2Father.cpp  Par_LB_Refine_SendParticle2Son.cpp

vpath %.cpp    Particle/LoadBalance
endif # LOAD_BALANCE

endif # PARTICLE


# parallelization source files (included only if "SERIAL" is turned off)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
CC_FILE     += Flu_AllocateFluxArray_Buffer.cpp

CC_FILE     += Flag_Buffer.cpp  Refine_Buffer.cpp

CC_FILE     += Buf_AllocateBufferPatch.cpp  Buf_AllocateBufferPatch_Base.cpp  Buf_GetBufferData.cpp \
               Buf_RecordExchangeDataPatchID.cpp  Buf_RecordExchangeFluxPatchID.cpp Buf_SortBoundaryPatch.cpp \
               Buf_RecordBoundaryFlag.cpp  Buf_RecordBoundaryPatch.cpp  Buf_RecordBoundaryPatch_Base.cpp \
               Buf_ResetBufferFlux.cpp

CC_FILE     += MPI_ExchangeBoundaryFlag.cpp  MPI_ExchangeBufferPosition.cpp  MPI_ExchangeData.cpp \
               Init_MPI.cpp  MPI_Exit.cpp

CC_FILE     += Output_BoundaryFlagList.cpp  Output_ExchangeDataPatchList.cpp  Output_ExchangeFluxPatchList.cpp \
               Output_ExchangePatchMap.cpp

CC_FILE     += Aux_Record_BoundaryPatch.cpp

vpath %.cpp    Buffer  MPI
endif # !SERIAL


# load-balance source files (included only if "LOAD_BALANCE" is turned on)
# ------------------------------------------------------------------------------------
CC_FILE     += LB_HilbertCurve.cpp  LB_Utility.cpp

ifeq "$(findstring -DLOAD_BALANCE, $(SIMU_OPTION))" "-DLOAD_BALANCE"
CC_FILE     += LB_Init_LoadBalance.cpp  LB_AllocateBufferPatch_Sibling.cpp  LB_RecordOvelapMPIPatchID.cpp \
               LB_Output_LBIdx.cpp  LB_AllocateBufferPatch_Father.cpp  LB_FindFather.cpp  LB_SiblingSearch.cpp \
               LB_RecordExchangeDataPatchID.cpp  LB_GetBufferData.cpp  LB_AllocateFluxArray.cpp \
               LB_RecordExchangeRestrictDataPatchID.cpp  LB_GrandsonCheck.cpp  LB_ExchangeFlaggedBuffer.cpp \
               LB_Refine.cpp  LB_Refine_GetNewRealPatchList.cpp  LB_Refine_AllocateNewPatch.cpp \
               LB_FindSonNotHome.cpp  LB_Refine_AllocateBufferPatch_Sibling.cpp \
               LB_AllocateBufferPatch_Sibling_Base.cpp  LB_RecordExchangeFixUpDataPatchID.cpp \
               LB_EstimateWorkload_AllPatchGroup.cpp  LB_EstimateLoadImbalance.cpp

endif # LOAD_BALANCE

vpath %.cpp    LoadBalance


# yt inline analysis source files (included only if "SUPPORT_LIBYT" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
CC_FILE     += YT_Init.cpp  YT_End.cpp  YT_SetParameter.cpp  YT_AddAllGrid.cpp  YT_Inline.cpp

vpath %.cpp    YT
endif # SUPPORT_LIBYT


# Grackle source files (included only if "SUPPORT_GRACKLE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
CUDA_FILE   += CUAPI_MemAllocate_Grackle.cu  CUAPI_MemFree_Grackle.cu

CC_FILE     += CPU_GrackleSolver_Original.cpp

CC_FILE     += Grackle_Init.cpp  Grackle_End.cpp  Init_MemAllocate_Grackle.cpp  End_MemFree_Grackle.cpp \
               Grackle_Prepare.cpp  Grackle_Close.cpp  Grackle_Init_FieldData.cpp  Grackle_AdvanceDt.cpp

vpath %.cpp    Grackle  Grackle/CPU_Grackle
endif # SUPPORT_GRACKLE


# star formation source files (included only if "STAR_FORMATION" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(filter -DSTAR_FORMATION, $(SIMU_OPTION))" "-DSTAR_FORMATION"
CC_FILE     += SF_CreateStar.cpp  SF_CreateStar_AGORA.cpp

vpath %.cpp    StarFormation
endif # STAR_FORMATION


# test problem source files
# --> note that here we simply compile all .cpp files under TestProblem/*/*/
# ------------------------------------------------------------------------------------
CC_FILE     += $(notdir $(wildcard TestProblem/*/*/*.cpp))

VPATH := $(dir $(wildcard TestProblem/*/*/))



# rules and targets
#######################################################################################################

# object files
# -------------------------------------------------------------------------------
OBJ_PATH := Object
OBJ := $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(CC_FILE))
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
OBJ += $(patsubst %.cu,  $(OBJ_PATH)/%.o, $(CUDA_FILE))
endif


# libraries
# -------------------------------------------------------------------------------
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
LIB += -L$(CUDA_PATH)/lib64
#LIB += -L$(CUDA_PATH)/lib
LIB += -lcudart
endif

ifeq "$(filter -DLAOHU, $(SIMU_OPTION))" "-DLAOHU"
LIB += -L$(GPUID_PATH) -lgpudevmgr
endif

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
   LIB += -L$(FFTW_PATH)/lib
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
      ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" "-DSERIAL"
         LIB += -ldrfftw -ldfftw
      else
         LIB += -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
      endif
   else
      ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" "-DSERIAL"
         LIB += -lsrfftw -lsfftw
      else
         LIB += -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
      endif
   endif
endif

ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
LIB += -L$(GRACKLE_PATH)/lib -lgrackle
endif

ifeq "$(filter -DSUPPORT_HDF5, $(SIMU_OPTION))" "-DSUPPORT_HDF5"
LIB += -L$(HDF5_PATH)/lib -lhdf5
endif

ifeq "$(filter -DSUPPORT_GSL, $(SIMU_OPTION))" "-DSUPPORT_GSL"
LIB += -L$(GSL_PATH)/lib -lgsl -lgslcblas
endif

ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
LIB += -L$(LIBYT_PATH)/lib -lyt
endif


# headers
# -------------------------------------------------------------------------------
INCLUDE := -I../include

ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
INCLUDE += -I$(MPI_PATH)/include
endif

ifeq "$(filter -DGRAVITY, $(SIMU_OPTION))" "-DGRAVITY"
INCLUDE += -I$(FFTW_PATH)/include
endif

ifeq "$(filter -DSUPPORT_GRACKLE, $(SIMU_OPTION))" "-DSUPPORT_GRACKLE"
INCLUDE += -I$(GRACKLE_PATH)/include
endif

ifeq "$(filter -DSUPPORT_HDF5, $(SIMU_OPTION))" "-DSUPPORT_HDF5"
INCLUDE += -I$(HDF5_PATH)/include
endif

ifeq "$(filter -DSUPPORT_GSL, $(SIMU_OPTION))" "-DSUPPORT_GSL"
INCLUDE += -I$(GSL_PATH)/include
endif

ifeq "$(filter -DSUPPORT_LIBYT, $(SIMU_OPTION))" "-DSUPPORT_LIBYT"
INCLUDE += -I$(LIBYT_PATH)/include
endif


# CXX flags
# -------------------------------------------------------------------------------
# remove the OpenMP flag if OPENMP is disabled
ifeq "$(filter -DOPENMP, $(SIMU_OPTION))" ""
   OPENMPFLAG =
endif

# fixes compilation issues on Intel MPI
ifeq "$(filter -DSERIAL, $(SIMU_OPTION))" ""
   CXXFLAG += -DMPICH_IGNORE_CXX_SEEK
endif

COMMONFLAG := $(INCLUDE) $(SIMU_OPTION)
CXXFLAG    += $(COMMONFLAG) $(OPENMPFLAG)


# NVCC flags
# -------------------------------------------------------------------------------
# common flags
NVCCFLAG_COM := -Xcompiler $(COMMONFLAG) -O3 #-Xopencc -OPT:Olimit=0

ifeq      "$(filter -DGPU_ARCH=FERMI,   $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
   NVCCFLAG_COM += -gencode arch=compute_20,code=\"compute_20,sm_20\"
#  NVCCFLAG_COM += -gencode arch=compute_20,code=\"compute_20,sm_21\"
else ifeq "$(filter -DGPU_ARCH=KEPLER,  $(SIMU_OPTION))" "-DGPU_ARCH=KEPLER"
   NVCCFLAG_COM += -gencode arch=compute_30,code=\"compute_30,sm_30\"
   NVCCFLAG_COM += -gencode arch=compute_35,code=\"compute_35,sm_35\"
   NVCCFLAG_COM += -gencode arch=compute_37,code=\"compute_37,sm_37\"
else ifeq "$(filter -DGPU_ARCH=MAXWELL, $(SIMU_OPTION))" "-DGPU_ARCH=MAXWELL"
   NVCCFLAG_COM += -gencode arch=compute_50,code=\"compute_50,sm_50\"
   NVCCFLAG_COM += -gencode arch=compute_52,code=\"compute_52,sm_52\"
else ifeq "$(filter -DGPU_ARCH=PASCAL,  $(SIMU_OPTION))" "-DGPU_ARCH=PASCAL"
   NVCCFLAG_COM += -gencode arch=compute_60,code=\"compute_60,sm_60\"
   NVCCFLAG_COM += -gencode arch=compute_61,code=\"compute_61,sm_61\"
else ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
   $(error unknown GPU_ARCH (please set it in the Makefile))
endif

ifeq "$(filter -DGAMER_DEBUG, $(SIMU_OPTION))" "-DGAMER_DEBUG"
   NVCCFLAG_COM += -g -Xptxas -v   # "-G" may cause the GPU Poisson solver to fail
endif

# fluid solver flags
NVCCFLAG_FLU += -Xptxas -dlcm=ca -prec-div=false -ftz=true #-prec-sqrt=false

ifeq      "$(filter -DGPU_ARCH=FERMI, $(SIMU_OPTION))" "-DGPU_ARCH=FERMI"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
#    NVCCFLAG_FLU += --maxrregcount=XX
   else
#    NVCCFLAG_FLU += --maxrregcount=XX
   endif
else ifeq "$(filter -DGPU_ARCH=KEPLER, $(SIMU_OPTION))" "-DGPU_ARCH=KEPLER"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=128
   else
     NVCCFLAG_FLU += --maxrregcount=64
   endif
else ifeq "$(filter -DGPU_ARCH=MAXWELL, $(SIMU_OPTION))" "-DGPU_ARCH=MAXWELL"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
else ifeq "$(filter -DGPU_ARCH=PASCAL, $(SIMU_OPTION))" "-DGPU_ARCH=PASCAL"
   ifeq "$(filter -DFLOAT8, $(SIMU_OPTION))" "-DFLOAT8"
     NVCCFLAG_FLU += --maxrregcount=192
   else
     NVCCFLAG_FLU += --maxrregcount=128
   endif
endif

# Poisson/gravity solvers flags
NVCCFLAG_POT += -Xptxas -dlcm=ca


# remove extra whitespaces
# -------------------------------------------------------------------------------
CXX          := $(strip $(CXX))
CXXFLAG      := $(strip $(CXXFLAG))
OPENMPFLAG   := $(strip $(OPENMPFLAG))
LIB          := $(strip $(LIB))
NVCC         := $(strip $(NVCC))
NVCCFLAG_COM := $(strip $(NVCCFLAG_COM))
NVCCFLAG_FLU := $(strip $(NVCCFLAG_FLU))
NVCCFLAG_POT := $(strip $(NVCCFLAG_POT))


# implicit rules (do NOT modify the order of the following rules)
# -------------------------------------------------------------------------------
# output detailed compilation commands or not
ifeq "$(COMPILE_VERBOSE)" "0"
ECHO = @
else
ECHO =
endif

# GPU codes
ifeq "$(filter -DGPU, $(SIMU_OPTION))" "-DGPU"
$(OBJ_PATH)/CUAPI_%.o : CUAPI_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) -o $@ -c $<

$(OBJ_PATH)/CUFLU_%.o : CUFLU_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_FLU) -o $@ -c $<

$(OBJ_PATH)/CUPOT_%.o : CUPOT_%.cu
	@echo "Compiling $<"
	$(ECHO)$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_POT) -o $@ -c $<
endif # GPU

# CPU codes
$(OBJ_PATH)/%.o : %.cpp
	@echo "Compiling $<"
	$(ECHO)$(CXX) $(CXXFLAG) -o $@ -c $<


# linking
# -------------------------------------------------------------------------------
$(EXECUTABLE) : $(OBJ)
	@echo "Linking"
ifeq "$(COMPILE_VERBOSE)" "1"
	$(CXX) -o $@ $^ $(LIB) $(OPENMPFLAG)
	@printf "\nCompiling GAMER --> Successful!\n\n"; \
	cp $(EXECUTABLE) ../bin/
else
	@$(CXX) -o $@ $^ $(LIB) $(OPENMPFLAG); \
	(if [ -e $@ ]; then \
		printf "\nCompiling GAMER --> Successful!\n\n"; \
		cp $(EXECUTABLE) ../bin/; \
	else \
		printf "\nCompiling GAMER --> Failed!\n\n"; \
	fi)
endif
	@rm ./*.linkinfo -f


# clean
# -------------------------------------------------------------------------------
.PHONY: clean
clean :
	@rm -f $(OBJ_PATH)/*
	@rm -f $(EXECUTABLE)
	@rm ./*.linkinfo -f



