


# name of the executable file
#######################################################################################################
EXECUTABLE := Dizzy



# simulation options
#######################################################################################################

# (a) common options (applied to all models)
# -------------------------------------------------------------------------------
# model: HYDRO/MHD/ELBDM/PAR_ONLY (PAR_ONLY must work with PARTICLE)
SIMU_OPTION += -DMODEL=ELBDM

# self-gravity
SIMU_OPTION += -DGRAVITY

# Poisson solver: SOR/MG (successive-overrelaxation/multigrid)
SIMU_OPTION += -DPOT_SCHEME=SOR

# individual time-step (default: shared time-step)
SIMU_OPTION += -DINDIVIDUAL_TIMESTEP

# comoving frame
#SIMU_OPTION += -DCOMOVING

# particles (must work with GRAVITY; set MODEL=PAR_ONLY for particle-only simulations)
#SIMU_OPTION += -DPARTICLE


# (b) hydro options
# ------------------------------------------------------------------------------------
ifeq "$(findstring MODEL=HYDRO, $(SIMU_OPTION))" "MODEL=HYDRO"
# hydrodynamic scheme: RTVD/WAF/MHM/MHM_RP/CTU
SIMU_OPTION += -DFLU_SCHEME=CTU

# scheme of spatial data reconstruction: PLM/PPM (piecewise-linear/piecewise-parabolic) ##USELESS IN RTVD/WAF##
SIMU_OPTION += -DLR_SCHEME=PPM

# Riemann solver: EXACT/ROE/HLLE/HLLC ##ALL ARE USELESS IN RTVD, HLLE/HLLC ARE USELESS IN WAF##
SIMU_OPTION += -DRSOLVER=ROE

# number of passively advected scalars
SIMU_OPTION += -DNPASSIVE=0

# star formation
#SIMU_OPTION += -DGALAXY


# (c) MHD options
# ------------------------------------------------------------------------------------
else ifeq "$(findstring MODEL=MHD, $(SIMU_OPTION))" "MODEL=MHD"


# (d) ELBDM options
# ------------------------------------------------------------------------------------
else ifeq "$(findstring MODEL=ELBDM, $(SIMU_OPTION))" "MODEL=ELBDM"
# enforce the mass conservation
SIMU_OPTION += -DCONSERVE_MASS

# 4-th order Laplacian
SIMU_OPTION += -DLAPLACIAN_4TH

# include the quartic self-interaction potential (must turn on GRAVITY as well, not work with COMOVING yet)
#SIMU_OPTION += -DQUARTIC_SELF_INTERACTION

endif # MODEL


# (e) PARTICLE options
# ------------------------------------------------------------------------------------


# (f) optimization and compilation options
# ------------------------------------------------------------------------------------
# GPU acceleration
SIMU_OPTION += -DGPU

# optimization mode
SIMU_OPTION += -DGAMER_OPTIMIZATION

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG

# measure the elapsing wall-clock time of different parts of the program
SIMU_OPTION += -DTIMING

# measure the elapsing wall-clock time of different parts of the GPU solvers (will disable CPU/GPU overlapping)
#SIMU_OPTION += -DTIMING_SOLVER

# intel compiler (default: GNU compiler)
SIMU_OPTION += -DINTEL

# double precision (not supported for the GPU + self-gravity mode in non-Fermi GPUs)
SIMU_OPTION += -DFLOAT8

# serial mode (in which no MPI libraries are required)
SIMU_OPTION += -DSERIAL

# out-of-core computation (deprecated)
#SIMU_OPTION += -DOOC

# load-balance parallelization (one must not turn on the option "SERIAL") (HILBERT)
#SIMU_OPTION += -DLOAD_BALANCE=HILBERT

# overlap MPI communication with computation (experimental, only supported in LOAD_BALANCE)
#SIMU_OPTION += -DOVERLAP_MPI

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP

# enable performance optimization in Fermi GPUs
SIMU_OPTION += -DFERMI

# work on the NAOC Laohu GPU cluster
#SIMU_OPTION += -DLAOHU

# support HDF5 format
#SIMU_OPTION += -DSUPPORT_HDF5



# simulation parameters
#######################################################################################################

NLEVEL    := 6        # maximum number of grid levels (including the base level)
MAX_PATCH := 1000000  # maximum number of patches at each level

NLEVEL    := $(strip $(NLEVEL))
MAX_PATCH := $(strip $(MAX_PATCH))

SIMU_PARA := -DNLEVEL=$(NLEVEL)  -DMAX_PATCH=$(MAX_PATCH)



# source files
#######################################################################################################

# common source files
# -------------------------------------------------------------------------------
# Cuda source files (compiled with nvcc)
CUDA_FILE   := CUAPI_Asyn_FluidSolver.cu  CUAPI_DiagnoseDevice.cu  CUAPI_MemAllocate_Fluid.cu \
               CUAPI_MemFree_Fluid.cu  CUAPI_Set_Default_GPU_Parameter.cu  CUAPI_SetDevice.cu \
               CUAPI_Synchronize.cu


# C/C++ source files (compiled with c++ compiler)
CC_FILE     := Main.cpp  EvolveLevel.cpp  InvokeSolver.cpp  Prepare_PatchData.cpp \
               InterpolateGhostZone.cpp

CC_FILE     += Aux_Check_Parameter.cpp  Aux_Check_Conservation.cpp  Aux_Check.cpp  Aux_Check_Finite.cpp \
               Aux_Check_FluxAllocate.cpp  Aux_Check_PatchAllocate.cpp  Aux_Check_ProperNesting.cpp \
               Aux_Check_Refinement.cpp  Aux_Check_Restrict.cpp  Aux_Error.cpp  Aux_GetCPUInfo.cpp \
               Aux_GetMemInfo.cpp  Aux_Message.cpp  Aux_PatchCount.cpp  Aux_TakeNote.cpp  Aux_Timing.cpp \
               Aux_Check_MemFree.cpp  Aux_RecordPerformance.cpp  Aux_CheckFileExist.cpp  Aux_Array.cpp

CC_FILE     += CPU_FluidSolver.cpp  Flu_AdvanceDt.cpp  Flu_Prepare.cpp  Flu_Close.cpp  Flu_FixUp.cpp \
               Flu_Restrict.cpp  Flu_AllocateFluxArray.cpp  Flu_BoundaryCondition_User.cpp

CC_FILE     += End_GAMER.cpp  End_MemFree.cpp  End_MemFree_Fluid.cpp  End_StopManually.cpp  End_TestProb.cpp \
               Init_BaseLevel.cpp  Init_GAMER.cpp  Init_Load_DumpTable.cpp \
               Init_Load_FlagCriteria.cpp  Init_Load_Parameter.cpp  Init_MemAllocate.cpp \
               Init_MemAllocate_Fluid.cpp  Init_Parallelization.cpp  Init_RecordBasePatch.cpp  Init_Refine.cpp \
               Init_Restart.cpp  Init_StartOver.cpp  Init_TestProb.cpp  Init_UM.cpp  Init_OpenMP.cpp \
               Init_Restart_HDF5.cpp \

CC_FILE     += Interpolate.cpp  Int_CQuadratic.cpp  Int_MinMod1D.cpp  Int_MinMod3D.cpp  Int_vanLeer.cpp \
               Int_Quadratic.cpp  Int_Table.cpp  Int_CQuartic.cpp  Int_Quartic.cpp

CC_FILE     += Mis_Check_Synchronization.cpp  Mis_GetTotalPatchNumber.cpp  Mis_GetTimeStep.cpp  Mis_Heapsort.cpp \
               Mis_BinarySearch.cpp  Mis_1D3DIdx.cpp  Mis_Matching.cpp  Mis_GetTimeStep_UserCriteria.cpp \
               Mis_dTime2dt.cpp  Mis_CoordinateTransform.cpp  Mis_BinarySearch_Real.cpp

CC_FILE     += Output_DumpData_Total.cpp  Output_DumpData.cpp  Output_DumpManually.cpp  Output_PatchMap.cpp \
               Output_DumpData_Part.cpp  Output_FlagMap.cpp  Output_Patch.cpp  Output_PreparedPatch_Fluid.cpp \
               Output_PatchCorner.cpp  Output_Flux.cpp  Output_TestProbErr.cpp  Output_BasePowerSpectrum.cpp \
               Output_DumpData_Total_HDF5.cpp

CC_FILE     += Flag_Real.cpp  Refine.cpp   SiblingSearch.cpp  SiblingSearch_Base.cpp  FindFather.cpp \
               Flag_UserCriteria.cpp  Flag_Check.cpp  Flag_Lohner.cpp  Flag_Region.cpp

CC_FILE     += Table_01.cpp  Table_02.cpp  Table_03.cpp  Table_04.cpp  Table_05.cpp  Table_06.cpp \
               Table_07.cpp

vpath %.cu     GPU_API
vpath %.cpp    GAMER  Init  Refine  Fluid  Interpolation  Tables  Output  Miscellaneous  Auxiliary


# hydrodynamic source files (included only if "MODEL=HYDRO")
# ------------------------------------------------------------------------------------
ifeq "$(findstring MODEL=HYDRO, $(SIMU_OPTION))" "MODEL=HYDRO"
CUDA_FILE   += CUFLU_GetMaxCFL.cu  CUFLU_FluidSolver_RTVD.cu  CUFLU_FluidSolver_WAF.cu  CUFLU_FluidSolver_MHM.cu \
               CUFLU_FluidSolver_CTU.cu

CC_FILE     += CPU_FluidSolver_RTVD.cpp  CPU_FluidSolver_WAF.cpp  CPU_FluidSolver_MHM.cpp \
               CPU_FluidSolver_CTU.cpp  CPU_Shared_DataReconstruction.cpp  CPU_Shared_FluUtility.cpp \
               CPU_Shared_ComputeFlux.cpp  CPU_Shared_FullStepUpdate.cpp \
               CPU_Shared_RiemannSolver_Exact.cpp  CPU_Shared_RiemannSolver_Roe.cpp \
               CPU_Shared_RiemannSolver_HLLE.cpp  CPU_Shared_RiemannSolver_HLLC.cpp

CC_FILE     += Hydro_Init_StartOver_AssignData.cpp  Hydro_Aux_Check_Negative.cpp  Hydro_GetTimeStep_Fluid.cpp \
               Hydro_Init_UM_AssignData.cpp  Hydro_BoundaryCondition_Outflow.cpp  Hydro_BoundaryCondition_Reflecting.cpp

vpath %.cu     Model_Hydro/GPU_Hydro
vpath %.cpp    Model_Hydro/CPU_Hydro  Model_Hydro

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
CUDA_FILE   += CUPOT_HydroGravitySolver.cu

CC_FILE     += CPU_HydroGravitySolver.cpp  Hydro_GetTimeStep_Gravity.cpp

vpath %.cu     Model_Hydro/GPU_HydroGravity
vpath %.cpp    Model_Hydro/CPU_HydroGravity
endif


# MHD source files (included only if "MODEL=MHD")
# -------------------------------------------------------------------------------
else ifeq "$(findstring MODEL=MHD, $(SIMU_OPTION))" "MODEL=MHD"

vpath %.cu     Model_MHD/GPU_MHD
vpath %.cpp    Model_MHD/CPU_MHD  Model_MHD

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"

vpath %.cu     Model_MHD/GPU_MHDGravity
vpath %.cpp    Model_MHD/CPU_MHDGravity
endif


# ELBDM source files (included only inf "MODEL=ELBDM")
# -------------------------------------------------------------------------------
else ifeq "$(findstring MODEL=ELBDM, $(SIMU_OPTION))" "MODEL=ELBDM"
CUDA_FILE   += CUFLU_ELBDMSolver.cu

CC_FILE     += CPU_ELBDMSolver.cpp  ELBDM_Init_StartOver_AssignData.cpp  ELBDM_Init_UM_AssignData.cpp \
               ELBDM_GetTimeStep_Fluid.cpp  ELBDM_Flag_EngyDensity.cpp  ELBDM_UnwrapPhase.cpp \
               ELBDM_GetTimeStep_Phase.cpp  ELBDM_SetTaylor3Coeff.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDM
vpath %.cpp    Model_ELBDM/CPU_ELBDM  Model_ELBDM

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
CUDA_FILE   += CUPOT_ELBDMGravitySolver.cu

CC_FILE     += CPU_ELBDMGravitySolver.cpp  ELBDM_GetTimeStep_Gravity.cpp

vpath %.cu     Model_ELBDM/GPU_ELBDMGravity
vpath %.cpp    Model_ELBDM/CPU_ELBDMGravity
endif

endif # MODEL


# self-gravity source files (included only if "GRAVITY" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
CUDA_FILE   += CUAPI_MemAllocate_PoissonGravity.cu  CUAPI_MemFree_PoissonGravity.cu \
               CUAPI_Asyn_PoissonGravitySolver.cu

CUDA_FILE   += CUPOT_PoissonSolver_SOR_10to14cube.cu  CUPOT_PoissonSolver_SOR_16to18cube.cu \
               CUPOT_PoissonSolver_MG.cu

CC_FILE     += CPU_PoissonGravitySolver.cpp  CPU_PoissonSolver_SOR.cpp  CPU_PoissonSolver_FFT.cpp \
               CPU_PoissonSolver_MG.cpp

CC_FILE     += Init_FFTW.cpp  Gra_Close.cpp  Gra_Prepare_Flu.cpp  Gra_Prepare_Pot.cpp  Gra_Prepare_Corner.cpp \
               Gra_AdvanceDt.cpp  Poi_Close.cpp  Poi_Prepare_Pot.cpp  Poi_Prepare_Rho.cpp \
               Output_PreparedPatch_Poisson.cpp  Init_MemAllocate_PoissonGravity.cpp  CPU_ExternalPot.cpp \
               End_MemFree_PoissonGravity.cpp  Init_Set_Default_SOR_Parameter.cpp \
               Init_Set_Default_MG_Parameter.cpp  Poi_GetAverageDensity.cpp  Init_GreenFuncK.cpp \
               Init_ExternalPot.cpp  Poi_BoundaryCondition_Extrapolation.cpp  CPU_ExternalGravity.cpp

vpath %.cu     SelfGravity/GPU_Poisson  SelfGravity/GPU_Gravity
vpath %.cpp    SelfGravity/CPU_Poisson  SelfGravity/CPU_Gravity  SelfGravity
endif # GRAVITY


# particle source files (included only if "PARTICLE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(findstring PARTICLE, $(SIMU_OPTION))" "PARTICLE"
CUDA_FILE   +=

CC_FILE     += Par_Init_Function.cpp  Par_Output_Particle.cpp  Par_FindHomePatch_Base.cpp  Par_PassParticle2Son.cpp \
               Par_Aux_Check_Particle.cpp  Par_PassParticle2Father.cpp  Par_CollectParticleFromDescendant.cpp \
               Par_MassAssignment.cpp  Par_UpdateParticle.cpp  Par_GetTimeStep_Velocity.cpp \
               Par_PassParticle2Sibling.cpp  Par_CountParticleInDescendant.cpp

vpath %.cu     Particle/GPU
vpath %.cpp    Particle/CPU  Particle
endif # PARTICLE


# galaxy source files (included only if "GALAXY" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(findstring GALAXY, $(SIMU_OPTION))" "GALAXY"
#CUDA_FILE   +=
CC_FILE     += Galaxy.cpp  Gal_Init_StartOver_AssignData.cpp

ifeq "$(findstring GPU, $(SIMU_OPTION))" ""
CC_FILE     += Gal_ExternalGravity_Burkert.cpp
endif

#vpath %.cu
vpath %.cpp    Galaxy
endif # GALAXY


# parallelization source files (included only if "SERIAL" is turned off)
# ------------------------------------------------------------------------------------
ifeq "$(findstring SERIAL, $(SIMU_OPTION))" ""
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

CC_FILE     += Aux_RecordBoundaryPatch.cpp

vpath %.cpp    Buffer  MPI
endif # !SERIAL


# load-balance source files (included only if "LOAD_BALANCE" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(findstring LOAD_BALANCE, $(SIMU_OPTION))" "LOAD_BALANCE"
CC_FILE     += LB_HilbertCurve.cpp  LB_Utility.cpp

CC_FILE     += LB_Init_LoadBalance.cpp  LB_AllocateBufferPatch_Sibling.cpp  LB_RecordOvelapMPIPatchID.cpp \
               LB_Output_LBIdx.cpp  LB_AllocateBufferPatch_Father.cpp  LB_FindFather.cpp  LB_SiblingSearch.cpp \
               LB_RecordExchangeDataPatchID.cpp  LB_GetBufferData.cpp  LB_AllocateFluxArray.cpp \
               LB_RecordExchangeRestrictDataPatchID.cpp  LB_GrandsonCheck.cpp  LB_ExchangeFlaggedBuffer.cpp \
               LB_Refine.cpp  LB_Refine_GetNewRealPatchList.cpp  LB_Refine_AllocateNewPatch.cpp \
               LB_FindSonNotHome.cpp  LB_Refine_AllocateBufferPatch_Sibling.cpp \
               LB_AllocateBufferPatch_Sibling_Base.cpp  LB_RecordExchangeFixUpDataPatchID.cpp \

endif # LOAD_BALANCE

vpath %.cpp    LoadBalance


# out-of-core source files (included only if "OOC" is turned on)
# ------------------------------------------------------------------------------------
ifeq "$(findstring OOC, $(SIMU_OPTION))" "OOC"
CC_FILE     += Init_OOC.cpp  OOC_SwapPointer.cpp  OOC_Dump.cpp  OOC_Load.cpp  OOC_Send.cpp  OOC_Recv.cpp \
               OOC_ExchangeBoundaryFlag.cpp  OOC_MPI_ExchangeBoundaryFlag.cpp  OOC_MPI_GetBufferData.cpp \
               OOC_GetAverageDensity.cpp  OOC_UpdateBuffer.cpp  OOC_SortPatch.cpp  OOC_SetTargetP.cpp \
               OOC_Implement.cpp  OOC_Init_SortPatch.cpp

CC_FILE     += OOC_Init_StartOver_ConstructAllLevels.cpp  OOC_Init_StartOver_Restrict.cpp \
               OOC_Init_Reload_LoadData.cpp  OOC_Init_Reload_FindFather.cpp  OOC_Output_DumpData_Part.cpp \
               OOC_Init_Reload_ConstructAllLevels.cpp  OOC_Output_DumpData_Total.cpp  OOC_Mis_GetMaxCFL.cpp \
               OOC_Integration_IndiviTimeStep_FluAdvanceDt.cpp  OOC_Init_Reload_Restrict.cpp \
               OOC_EvolveLevel.cpp  OOC_Integration_IndiviTimeStep_FixUp.cpp \
               OOC_Integration_IndiviTimeStep_Flag.cpp  OOC_Integration_IndiviTimeStep_Refine.cpp \
               OOC_GetBufferData.cpp

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
CC_FILE     += OOC_Init_GAMER_GetPot.cpp  OOC_Integration_IndiviTimeStep_GraAdvanceDt.cpp \
               OOC_Gra_AdvanceDt_InvokeSolver.cpp  OOC_CPU_PoissonSolver_FFT_Cube_to_Slice.cpp \
               OOC_CPU_PoissonSolver_FFT_Slice_to_Cube.cpp  OOC_Mis_GetMaxAcc.cpp
endif

vpath %.cpp    OOC_API  OOC_Alternative
endif # OOC



# rules and targets
#######################################################################################################

# location of libraries
# -------------------------------------------------------------------------------
# hulk (openmpi-intel-qlc)
CUDA_TOOLKIT_PATH := /opt/gpu/cuda/default
FFTW_PATH         := /opt/math/fftw/2.1.5-intel
MPI_PATH          := /opt/mpi/openmpi/1.6.0-intel
HDF5_PATH         := /opt/hdf5/1.8.9

# hulk (openmpi-gcc-qlc)
#CUDA_TOOLKIT_PATH := /opt/gpu/cuda/4.2
#FFTW_PATH         := /opt/math/fftw/2.1.5-gcc-qlc
#MPI_PATH          := /opt/mpi/openmpi/1.8.1-gcc-qlc
#HDF5_PATH         := /opt/hdf5/1.8.9

# geisha (intel)
#CUDA_TOOLKIT_PATH := /usr/local/cuda
#FFTW_PATH         := /opt/math/fftw_openmpi_intel
#MPI_PATH          := /opt/mpi/openmpi_intel

# geisha (gcc)
#CUDA_TOOLKIT_PATH := /usr/local/cuda
#FFTW_PATH         := /opt/math/fftw_openmpi
#MPI_PATH          := /opt/mpi/openmpi

# dirac at NERSC (OpenMPI-1.6.0)
#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/4.1
#FFTW_PATH         := /project/projectdirs/gamer/hsiyu/Software/fftw/2.1.5-openmpi-1.6.0-intel
#MPI_PATH          := /usr/common/usg/openmpi/1.6/intel

# dirac at NERSC (without thread support)
#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/4.1
#FFTW_PATH         := /project/projectdirs/gamer/hsiyu/Software/fftw/2.1.5-openmpi-1.4.2-intel
#MPI_PATH          := /usr/common/usg/openmpi/1.4.2/intel

# dirac at NERSC (with thread support)
#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/4.1
#FFTW_PATH         := /global/scratch/sd/hsiyu/fftw.2.1.5-intel-thread
#MPI_PATH          := /global/scratch/sd/hsiyu/OpenMPI.1.5.3_Intel_Thread

# dirac at NERSC (OpenMPI-1.5.4)
#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/4.1
#FFTW_PATH         := /project/projectdirs/gamer/hsiyu/Software/fftw/2.1.5-openmpi-1.5.4-intel
#MPI_PATH          := /project/projectdirs/gamer/hsiyu/Software/mpi/openmpi/1.5.4-intel

# dirac at NERSC (MVAPICH2-1.5.1p1)
#CUDA_TOOLKIT_PATH := /usr/common/usg/cuda/4.1
#FFTW_PATH         := /project/projectdirs/gamer/hsiyu/Software/fftw/2.1.5-mvapich2-1.5.1p1-gcc
#MPI_PATH          := /usr/common/usg/mvapich2/1.5.1p1/gcc

# laohu at NAOC (OpenMPI-1.3.2-intel)
#CUDA_TOOLKIT_PATH := $(CUDA_INSTALL_PATH)
#FFTW_PATH         := /ifs/data/hsiyu/Software/fftw/2.1.5-openmpi-1.3.2-intel
#MPI_PATH          := /usr/mpi/intel/openmpi-1.3.2-qlc
#GPUID_PATH        := /ifs/data/hsiyu/Software/gpudevmgr

# formosa4 (fish:mvapich2-1.8a2-intel)
#CUDA_TOOLKIT_PATH := /usr/local/cuda
#FFTW_PATH         := /home/u11thc00/fish/Software/fftw/2.1.5-mvapich2-1.8a2-intel
#MPI_PATH          := /home/u11thc00/fish/Software/mpi/mvapich2/1.8a2-intel

# coxhydrae at ARI
#CUDA_TOOLKIT_PATH := /usr/local/cuda
#FFTW_PATH         := /data/hsiyu/Software/fftw/2.1.5
#FFTW_PATH         := /home/liulei/local
#MPI_PATH          := /usr


# compilers
# -------------------------------------------------------------------------------
NVCC := $(CUDA_TOOLKIT_PATH)/bin/nvcc
ifeq "$(findstring SERIAL, $(SIMU_OPTION))" "SERIAL"
   ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
   CXX := icpc
   else
   CXX := g++
   endif
else
#  CXX := $(MPI_PATH)/bin/mpiCC
   CXX := $(MPI_PATH)/bin/mpicxx
endif


# object files
# -------------------------------------------------------------------------------
OBJ_PATH := Object
OBJ := $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(CC_FILE))
ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
OBJ += $(patsubst %.cu,  $(OBJ_PATH)/%.o, $(CUDA_FILE))
endif


# libraries
# -------------------------------------------------------------------------------
LIB :=
ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
LIB += -L$(CUDA_TOOLKIT_PATH)/lib64
#LIB += -L$(CUDA_TOOLKIT_PATH)/lib
LIB += -lcudart
endif

ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
LIB += -limf
endif

ifeq "$(findstring LAOHU, $(SIMU_OPTION))" "LAOHU"
LIB += -L$(GPUID_PATH) -lgpudevmgr
endif

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
   LIB += -L$(FFTW_PATH)/lib
   ifeq "$(findstring FLOAT8, $(SIMU_OPTION))" "FLOAT8"
      ifeq "$(findstring SERIAL, $(SIMU_OPTION))" "SERIAL"
         LIB += -ldrfftw -ldfftw
      else
         LIB += -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
      endif
   else
      ifeq "$(findstring SERIAL, $(SIMU_OPTION))" "SERIAL"
         LIB += -lsrfftw -lsfftw
      else
         LIB += -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
      endif
   endif
endif

ifeq "$(findstring OOC, $(SIMU_OPTION))" "OOC"
LIB += -laio
endif

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" "OPENMP"
   ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
      OPENMP := -openmp
   else
      OPENMP := -fopenmp
   endif
endif

ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
LIB += -L$(HDF5_PATH)/lib -lhdf5
endif


# headers
# -------------------------------------------------------------------------------
INCLUDE := -I../include

ifeq "$(findstring SERIAL, $(SIMU_OPTION))" ""
INCLUDE += -I$(MPI_PATH)/include
endif

ifeq "$(findstring GRAVITY, $(SIMU_OPTION))" "GRAVITY"
INCLUDE += -I$(FFTW_PATH)/include
endif

ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
INCLUDE += -I$(HDF5_PATH)/include
endif


# CXX flags
# -------------------------------------------------------------------------------
COMMONFLAG := $(INCLUDE) $(SIMU_OPTION) $(SIMU_PARA)

ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
CXXWARN_FLAG := -w1
else
CXXWARN_FLAG := -Wextra -Wall -Wimplicit -Wswitch -Wformat -Wchar-subscripts -Wparentheses -Wmultichar \
                -Wtrigraphs -Wpointer-arith -Wcast-align -Wreturn-type -Wno-unused-function
endif

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" ""
CXXWARN_FLAG += -Wno-unknown-pragmas
endif

ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
CXXFLAG  := $(CXXWARN_FLAG) $(COMMONFLAG) $(OPENMP) -O3 -xSSE4.2 -fp-model precise #-ipo
else
CXXFLAG  := $(CXXWARN_FLAG) $(COMMONFLAG) $(OPENMP) -O3
endif

ifeq "$(findstring GAMER_DEBUG, $(SIMU_OPTION))" "GAMER_DEBUG"
   ifeq "$(findstring INTEL, $(SIMU_OPTION))" "INTEL"
      CXXFLAG += -g -debug
   else
      CXXFLAG += -g
   endif
endif


# NVCC flags
# -------------------------------------------------------------------------------
NVCCFLAG_COM := -Xcompiler $(COMMONFLAG) -O3 #-Xopencc -OPT:Olimit=0

NVCCFLAG_COM += -gencode arch=compute_20,code=\"compute_20,sm_20\"
ifeq "$(findstring FERMI, $(SIMU_OPTION))" ""
   NVCCFLAG_COM += -gencode arch=compute_13,code=\"compute_13,sm_13\"
   ifeq "$(findstring FLOAT8, $(SIMU_OPTION))" ""
      NVCCFLAG_COM += -gencode arch=compute_10,code=\"compute_10,sm_10\"
   endif
endif

ifeq "$(findstring FERMI, $(SIMU_OPTION))" "FERMI"
   NVCCFLAG_FLU += -Xptxas -dlcm=ca -prec-div=false -ftz=true #-prec-sqrt=false
   NVCCFLAG_POT += -Xptxas -dlcm=ca
   ifeq "$(findstring FLOAT8, $(SIMU_OPTION))" "FLOAT8"
#     NVCCFLAG_FLU += --maxrregcount=64
   else
#     NVCCFLAG_FLU += --maxrregcount=XX
   endif
else # non-Fermi GPUs
   ifeq "$(findstring FLOAT8, $(SIMU_OPTION))" "FLOAT8"
#     NVCCFLAG_FLU += --maxrregcount=120
   else
      NVCCFLAG_FLU += --maxrregcount=64
   endif
endif

ifeq "$(findstring GAMER_DEBUG, $(SIMU_OPTION))" "GAMER_DEBUG"
NVCCFLAG_COM += -D_DEBUG -g -G -Xptxas -v #-deviceemu
endif


# implicit rules (do NOT modify the order of the following rules)
# -------------------------------------------------------------------------------
ifeq "$(findstring GPU, $(SIMU_OPTION))" "GPU"
$(OBJ_PATH)/CUAPI_%.o : CUAPI_%.cu
	$(NVCC) $(NVCCFLAG_COM) -o $@ -c $<

$(OBJ_PATH)/CUFLU_%.o : CUFLU_%.cu
	$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_FLU) -o $@ -c $<

$(OBJ_PATH)/CUPOT_%.o : CUPOT_%.cu
	$(NVCC) $(NVCCFLAG_COM) $(NVCCFLAG_POT) -o $@ -c $<
endif

$(OBJ_PATH)/%.o : %.cpp
	$(CXX) $(CXXFLAG) -o $@ -c $<


# linking
# -------------------------------------------------------------------------------
$(EXECUTABLE) : $(OBJ)
	$(CXX) -o $@ $^ $(LIB) $(OPENMP)
	rm ./*.linkinfo -f
	cp $(EXECUTABLE) ../bin/Run/


# clean
# -------------------------------------------------------------------------------
clean :
	rm -f $(OBJ_PATH)/*
	rm -f $(EXECUTABLE)
	rm ./*.linkinfo -f



