
# file names
#######################################################################################################
EXECUTABLE := GAMER_Data2Visit



# siimulation options
#######################################################################################################
# model: HYDRO/MHD/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8

# debug mode
#SIMU_OPTION += -DGAMER_DEBUG



# siimulation parameters
#######################################################################################################
NLEVEL     := 6         # level : 0 ~ NLEVEL-1
MAX_PATCH  := 1000000   # maximum number of patches at each level

NLEVEL     := $(strip $(NLEVEL))
MAX_PATCH  := $(strip $(MAX_PATCH))

SIMU_PARA := -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH)



# sources
# #######################################################################################################
SOURCE := Main.cpp

SOURCE += LoadData.cpp  LoadData_HDF5.cpp  Aux_CheckFileExist.cpp  Aux_Message.cpp  Aux_Error.cpp

vpath %.cpp ./ GAMER_Functions



# rules and targets
#######################################################################################################
#HDF5_PATH := /opt/hdf5/1.8.9
#SILO_PATH := /opt/visualization/silo/default
HDF5_PATH := /projects/ncsa/grav/softwares/miniconda3
SILO_PATH := /projects/ncsa/grav/softwares/silo/4.10.2-hdf51.8.15


CC    := icpc
CFLAG := -O3 -w1 -xSSE4.2 -fp-model precise

#CC    := g++
#CFLAG := -O3 -Wall

ifeq "$(findstring GAMER_DEBUG, $(SIMU_OPTION))" "GAMER_DEBUG"
CFLAG += -g -debug
#CFLAG += -g
endif


INCLUDE := -I./Header -I$(SILO_PATH)/include -I$(HDF5_PATH)/include
LIB     := -L$(SILO_PATH)/lib -lsiloh5 -lz -L$(HDF5_PATH)/lib -lhdf5

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
