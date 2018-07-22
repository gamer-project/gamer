


# name of the executable file 
#######################################################################################################
EXECUTABLE  = GAMER_ExtractUniform



# siimulation options
#######################################################################################################
# model: HYDRO/MHD/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8

# enable OpenMP parallelization
SIMU_OPTION += -DOPENMP

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG

# support HDF5 format
SIMU_OPTION += -DSUPPORT_HDF5

# serial version (no MPI is required)
#SIMU_OPTION += -DSERIAL



# siimulation parameters
#######################################################################################################
NLEVEL        = 10        # level : 0 ~ NLEVEL-1 
NCOMP_PASSIVE = 0         # number of passive scalars
MAX_PATCH     = 1000000   # the maximum number of patches in each level

NLEVEL        := $(strip $(NLEVEL))
NCOMP_PASSIVE := $(strip $(NCOMP_PASSIVE))
MAX_PATCH     := $(strip $(MAX_PATCH))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH) -DNCOMP_PASSIVE=$(NCOMP_PASSIVE)



# sources
#######################################################################################################
SOURCE  = Main.cpp  LoadData.cpp  End_MemFree.cpp  FindFather.cpp  Flu_Restrict.cpp  Init_MemAllocate.cpp \
          Init_RecordBasePatch.cpp  SiblingSearch_Base.cpp  SiblingSearch.cpp  Table_01.cpp  Table_02.cpp \
          Table_03.cpp  Table_04.cpp  Table_05.cpp  Table_07.cpp  Aux_Message.cpp  Aux_Error.cpp \
          Hydro_BoundaryCondition_Outflow.cpp LoadData_HDF5.cpp  Aux_CheckFileExist.cpp

SOURCE += Interpolate.cpp  Int_CQuadratic.cpp  Int_MinMod1D.cpp  Int_MinMod3D.cpp  Int_vanLeer.cpp \
          Int_Quadratic.cpp  Int_Table.cpp  Int_CQuartic.cpp  Int_Quartic.cpp

ifeq "$(findstring SERIAL, $(SIMU_OPTION))" ""
SOURCE += Buf_AllocateBufferPatch_Base.cpp  Buf_AllocateBufferPatch.cpp  Buf_GetBufferData.cpp \
          Buf_RecordBoundaryPatch_Base.cpp  Buf_RecordBoundaryPatch.cpp  Buf_RecordExchangeDataPatchID.cpp \
          Buf_SortBoundaryPatch.cpp  MPI_Exit.cpp  MPI_ExchangeBufferPosition.cpp  MPI_ExchangeInfo.cpp
endif

ifeq "$(findstring MODEL=ELBDM, $(SIMU_OPTION))" "MODEL=ELBDM"
SOURCE += ELBDM_UnwrapPhase.cpp
endif

vpath %.cpp ./ GAMER_Functions GAMER_Functions/Interpolation



# rules and targets
#######################################################################################################
#HDF5_PATH := /opt/hdf5/1.8.9
#MPI_PATH  := /usr/mpi/intel/openmpi-1.4.3-qlc
#MPI_PATH  := /opt/mpi/openmpi/1.6.0-intel
MPI_PATH  := /usr/local/mpi/openmpi-1.8.4-intel-15.0
HDF5_PATH := /projects/ncsa/grav/softwares/miniconda2


ifeq "$(findstring SERIAL, $(SIMU_OPTION))" ""
CXX   := $(MPI_PATH)/bin/mpicxx
else
CXX   := icpc
endif

CFLAG := -O3 -w1 -xSSE4.2 -fp-model precise
#CFLAG := -O3 -Wall

ifeq "$(findstring GAMER_DEBUG, $(SIMU_OPTION))" "GAMER_DEBUG"
CFLAG += -g
endif

ifeq "$(findstring OPENMP, $(SIMU_OPTION))" "OPENMP"
CFLAG += -fopenmp
else
CFLAG += -Wno-unknown-pragmas
endif

INCLUDE := -I./Header
ifeq "$(findstring SERIAL, $(SIMU_OPTION))" ""
INCLUDE += -I$(MPI_PATH)/include
endif
ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
INCLUDE += -I$(HDF5_PATH)/include
endif

LIB := -limf
ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
LIB += -L$(HDF5_PATH)/lib -lhdf5
endif

OBJ_PATH := ./Object
OBJ 	 := $(patsubst %.cpp, $(OBJ_PATH)/%.o, $(SOURCE))


$(OBJ_PATH)/%.o : %.cpp
	$(CXX) $(CFLAG) $(SIMU_PARA) $(SIMU_OPTION) $(INCLUDE) -o $@ -c $<

$(EXECUTABLE): $(OBJ)
	$(CXX) $(CFLAG) -o $@ $^ $(LIB)
	cp $(EXECUTABLE) ./Run/

clean:
	rm -f $(OBJ)
	rm -f $(EXECUTABLE)
