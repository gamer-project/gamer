
# file names
#######################################################################################################
PROGRAM    := GAMER_MakeDataInOrder
EXECUTABLE := GAMER_MakeDataInOrder



# siimulation options
#######################################################################################################
# double precision
#SIMU_OPTION += -DFLOAT8



# rules and targets
#######################################################################################################
CC    := icpc
CFLAG := -O3

$(EXECUTABLE): $(PROGRAM).o
	 $(CC) $(CFLAG) -o $@ $<

$(PROGRAM).o: $(PROGRAM).cpp
	 $(CC) $(CFLAG) $(SIMU_OPTION) -o $@ -c $<

clean:
	 rm -f *.o
	 rm -f $(EXECUTABLE)

