GRACKLE_DIR = /usr
CC          = gcc
EXE         = coolingrate_proper

all: $(EXE)

$(EXE): coolingrate_proper.c
	$(CC) -I$(GRACKLE_DIR)/include -L$(GRACKLE_DIR)/lib -fopenmp -O2 -o $(EXE) coolingrate_proper.c -lgrackle -lm

run: $(EXE)
	./$(EXE) 	0  4  8  0.02
	./$(EXE) 	1  4  8  0.02
	./$(EXE) 	2  4  8  0.02
	./$(EXE) 	3  4  8  0.02
	./$(EXE) 	4  4  8  0.02
	./$(EXE) 	5  4  8  0.02
	./$(EXE) 	6  4  8  0.02
	./$(EXE) 	7  4  8  0.02
	./$(EXE) 	8  4  8  0.02
	./$(EXE) 	9  4  8  0.02

plot: 
	python Compute_CoolingCurve.py

clean:
	rm -f $(EXE) coolingrate*.dat