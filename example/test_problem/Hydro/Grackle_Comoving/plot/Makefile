GRACKLE_DIR = /usr
CC = gcc

all: coolingrate_proper

coolingrate_proper: coolingrate_proper.c
	$(CC) -I$(GRACKLE_DIR)/include -L$(GRACKLE_DIR)/lib -fopenmp -O2 -o coolingrate_proper coolingrate_proper.c -lgrackle -lm

run: coolingrate_proper
	./coolingrate_proper 	0  4  8  0.02
	./coolingrate_proper 	1  4  8  0.02
	./coolingrate_proper 	2  4  8  0.02
	./coolingrate_proper 	3  4  8  0.02
	./coolingrate_proper 	4  4  8  0.02
	./coolingrate_proper 	5  4  8  0.02
	./coolingrate_proper 	6  4  8  0.02
	./coolingrate_proper 	7  4  8  0.02
	./coolingrate_proper 	8  4  8  0.02
	./coolingrate_proper 	9  4  8  0.02

plot: 
	python Compute_CoolingCurve.py

clean:
	rm -f coolingrate_proper coolingrate*.dat