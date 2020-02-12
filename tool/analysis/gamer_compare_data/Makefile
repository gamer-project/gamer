


# file names
#######################################################################################################
PROGRAM    = GAMER_CompareData
EXECUTABLE = GAMER_CompareData



# siimulation options
#######################################################################################################
# model: HYDRO/ELBDM
SIMU_OPTION += -DMODEL=HYDRO

# double precision
#SIMU_OPTION += -DFLOAT8



# siimulation parameters
#######################################################################################################
NLEVEL        = 10        # level : 0 ~ NLEVEL-1
NCOMP_PASSIVE = 0         # number of passive scalars
MAX_PATCH     = 1000000   # maximum number of patches in each level

NLEVEL        := $(strip $(NLEVEL))
NCOMP_PASSIVE := $(strip $(NCOMP_PASSIVE))
MAX_PATCH     := $(strip $(MAX_PATCH))

SIMU_PARA = -DNLEVEL=$(NLEVEL) -DMAX_PATCH=$(MAX_PATCH) -DNCOMP_PASSIVE=$(NCOMP_PASSIVE)



# rules and targets
#######################################################################################################
CC    := icpc
CFLAG := -O3 -w1
CFLAG += -g


$(EXECUTABLE): $(PROGRAM).o
	$(CC) $(CFLAG) -o $@ $<
	cp $(EXECUTABLE) ./Run/

$(PROGRAM).o: $(PROGRAM).cpp
	$(CC) $(CFLAG) $(SIMU_PARA) $(SIMU_OPTION) -o $@ -c $<

clean:
	rm -f *.o
	rm -f $(EXECUTABLE)

