#MACRO += -DXRAY_ROSAT
#MACRO += -DXRAY_EROSITA
#MACRO += -DHADRONIC_GAMMARAY
#MACRO += -DLEPTONIC_GAMMARAY
#MACRO += -DSYNCHROTRON


#MACRO += -DNUM_THREADS=16

#MACRO += -DDEBUG

#MACRO += -DFLOAT8


CC := mpic++

CFLAGS     = -Wall -O3 -std=c++11
DEBUG_FLAG = -g
OPENMPFLAG = -fopenmp


##########
# TW3
##########
#INC := -I/work/tseng0929/software/hdf5/hdf5-gnu-7.3/include \
#       -I/work/tseng0929/software/cfitsio/include \
#       -I/work/tseng0929/software/CCfits/include \
#       -I/work/tseng0929/software/gsl/2.7.1/include
#
#LIB := -L/work/tseng0929/software/hdf5/hdf5-gnu-7.3/lib -lhdf5 \
#       -L/work/tseng0929/software/cfitsio/lib -lcfitsio \
#       -L/work/tseng0929/software/CCfits/lib  -lCCfits \
#       -L/work/tseng0929/software/gsl/2.7.1/lib -lgsl -lgslcblas

###########
# Eureka
###########
#export PATH=/software/hdf5/default/bin:$PATH
#export PATH=/software/gsl/default/bin:$PATH
#export PATH=/software/gcc/9.3.0/bin:$PATH
#export PATH=/software/openmpi/4.1.1-ucx_mt-gnu-9.3.0/bin:$PATH
#
#export LD_LIBRARY_PATH=/software/hdf5/default/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/projectW/tseng/opt/cfitsio/cfitsio-4.0.0/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/projectW/tseng/opt/CCfits/CCfits-2.6/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/projectW/tseng/opt/HEALPix/HEALPix-3.80/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/software/intel/oneapi/compiler/2021.1.1/linux/compiler/lib/intel64_lin:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/projectW/tseng/opt/openmpi-gcc-9.3.0/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/software/gsl/default/lib:$LD_LIBRARY_PATH

INC := -I /software/hdf5/default/include \
       -I /projectW/tseng/opt/cfitsio/cfitsio-4.0.0/include \
       -I /projectW/tseng/opt/CCfits/CCfits-2.6/include \
       -I /software/gsl/default/include

LIB := -L /software/hdf5/default/lib -lhdf5 \
       -L /projectW/tseng/opt/cfitsio/cfitsio-4.0.0/lib -lcfitsio \
       -L /projectW/tseng/opt/CCfits/CCfits-2.6/lib  -lCCfits \
       -L /software/gsl/default/lib -lgsl -lgslcblas

EXECUTABLE := Project



## source files
CC_FILE  = Main.c Array.c BinarySearch.c Cooling.c Interpolation.c ReadTable.c \
           Make_Projection.c Make_Slice.c XRayMap.c HaloDensFun.c Leptonic.c Hadronic.c \
           Synchrotron.c Utilities.c
SRC_PATH = src


## object files
OBJ_FILE = $(CC_FILE:.c=.o)
OBJ_PATH = objective
OBJ      = $(patsubst %,$(OBJ_PATH)/%,$(OBJ_FILE))


# header files
HEADER_FILE = Macro.h Prototypes.h General.h
INC_PATH    = include
HEADER      = $(patsubst %,$(INC_PATH)/%,$(HEADER_FILE))

ALL_FLAG = $(DEBUG_FLAG) $(OPENMPFLAG) $(LIB) $(INC) $(CFLAGS) $(MACRO)

# Compiling
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c $(HEADER)
	$(CC) $(ALL_FLAG) -o $@ -c $<

# Linking
$(EXECUTABLE): $(OBJ)
	@echo "Linking ..."
	$(CC) -o $@ $^ $(LIB) $(OPENMPFLAG)
	@echo "Done!"
	@mv $(EXECUTABLE) bin/

.PHONY:
clean:
	@rm -rf ${OBJ_PATH}/*.o
	@rm -rf ${EXECUTABLE}
	@rm -rf *.o rm ${EXECUTABLE}
