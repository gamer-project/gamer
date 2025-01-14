


# executable file
#######################################################################################################
EXECUTABLE := GAMER_CompareData



# simulation options
#######################################################################################################
# model: HYDRO/ELBDM
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
MAX_PATCH     := 2000000   # maximum number of patches in each level

NLEVEL        := $(strip $(NLEVEL))
MAX_PATCH     := $(strip $(MAX_PATCH))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH)



# sources
# #######################################################################################################
SOURCE  = Main.cpp

SOURCE += LoadData.cpp  LoadData_HDF5.cpp  Aux_AllocateArray2D.cpp  Aux_CheckFileExist.cpp  Aux_Error.cpp \
          Aux_Message.cpp  SortParticle.cpp

vpath %.cpp ./ GAMER_Functions



# rules and targets
#######################################################################################################
HDF5_PATH := /software/hdf5/default
#HDF5_PATH := /software/hdf5/1.10.6-intel-2023.1.0

CC    := icpc
CFLAG := -O3 -w1 -g -diag-disable=10441

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
