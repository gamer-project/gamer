


# name of the executable file
#######################################################################################################
EXECUTABLE := GAMER_ExtractProfile



# siimulation options
#######################################################################################################
# model: HYDRO/MHD/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8

# double precision for particle floating-point attributes
#SIMU_OPTION += -DFLOAT8_PAR

# long integer for particle integer attributes
SIMU_OPTION += -DINT8_PAR

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG

# support HDF5 format
SIMU_OPTION += -DSUPPORT_HDF5



# siimulation parameters
#######################################################################################################
NLEVEL        := 10        # level : 0 ~ NLEVEL-1
NCOMP_PASSIVE := 0         # number of passive scalars
MAX_PATCH     := 2000000   # the maximum number of patches in each level

NLEVEL        := $(strip $(NLEVEL))
NCOMP_PASSIVE := $(strip $(NCOMP_PASSIVE))
MAX_PATCH     := $(strip $(MAX_PATCH))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH) -DNCOMP_PASSIVE=$(NCOMP_PASSIVE)



# sources
# #######################################################################################################
SOURCE  = Main.cpp

SOURCE += LoadData.cpp  Init_RecordBasePatch.cpp  FindFather.cpp  Table_01.cpp  Table_02.cpp  SiblingSearch.cpp \
          SiblingSearch_Base.cpp  Flu_Restrict.cpp  Aux_Error.cpp  Aux_Message.cpp  Prepare_PatchData.cpp \
          InterpolateGhostZone.cpp  Prepare_PatchData.cpp  InterpolateGhostZone.cpp  Table_03.cpp  Table_04.cpp \
          Output_Patch.cpp  LoadData_HDF5.cpp  Aux_CheckFileExist.cpp

SOURCE += Interpolate.cpp  Int_CQuadratic.cpp  Int_MinMod1D.cpp  Int_MinMod3D.cpp  Int_vanLeer.cpp \
          Int_Quadratic.cpp  Int_Table.cpp  Int_CQuartic.cpp  Int_Quartic.cpp

ifeq "$(findstring MODEL=ELBDM, $(SIMU_OPTION))" "MODEL=ELBDM"
SOURCE += ELBDM_UnwrapPhase.cpp
endif

vpath %.cpp ./ GAMER_Functions GAMER_Functions/Interpolation



# rules and targets
#######################################################################################################
HDF5_PATH := /software/hdf5/default

CC    := icpc
CFLAG := -O3 -w1 -g #-xSSE4.2 -fp-model precise
#CFLAG := -O3 -w1 -mp1

#CC    := g++
#CFLAG := -O3 -Wall

#LIB   := -limf

ifeq "$(findstring GAMER_DEBUG, $(SIMU_OPTION))" "GAMER_DEBUG"
CFLAG += -g
endif

INCLUDE  := -I./Header
ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
INCLUDE += -I$(HDF5_PATH)/include
endif

LIB := $(strip $(LIB))
ifeq "$(findstring SUPPORT_HDF5, $(SIMU_OPTION))" "SUPPORT_HDF5"
LIB += -L$(HDF5_PATH)/lib -lhdf5
endif

OBJ_DIR := ./Object
OBJ     := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(SOURCE))


$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CFLAG) $(SIMU_PARA) $(SIMU_OPTION) $(INCLUDE) -o $@ -c $<

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAG) -o $@ $^ $(LIB)
	cp $(EXECUTABLE) ./Run/

clean:
	rm -f $(OBJ)
	rm -f $(EXECUTABLE)
