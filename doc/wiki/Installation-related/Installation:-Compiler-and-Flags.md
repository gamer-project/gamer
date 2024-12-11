To choose a compiler and compilation flags, set the following
variables in the `Makefile`:
```Makefile
CXX        =    # C++ compiler
CXXFLAG    =    # compilation flags
OPENMPFLAG =    # openmp flag
LIB        =    # libraries and linker flags
NVCC       =    # CUDA compiler
```

Example: Intel compiler
```Makefile
CXX        = $(MPI_PATH)/bin/mpicxx   # replace by "icpc" in the serial mode
CXXFLAG    = -g -O3 -w1
OPENMPFLAG = -fopenmp
LIB        = -limf
NVCC       = $(CUDA_PATH)/bin/nvcc
```

Example: GNU compiler
```Makefile
CXX        = $(MPI_PATH)/bin/mpicxx   # replace by "g++" in the serial mode
CXXFLAG    = -g -O3 -Wall -Wextra
CXXFLAG   += -Wno-unused-variable -Wno-unused-parameter \
             -Wno-maybe-uninitialized -Wno-unused-but-set-variable \
             -Wno-unused-result -Wno-unused-function
OPENMPFLAG = -fopenmp
LIB        =
NVCC       = $(CUDA_PATH)/bin/nvcc
```

On a Cray system (even when adopting the Intel or GNU compiler), set `CXX` and `NVCC` as
```Makefile
CXX        = CC
NVCC       = nvcc -ccbin CC
```

The prefixes `$(MPI_PATH)/bin/` and `$(CUDA_PATH)/bin/` in the
above examples are optional to force the Makefile to use the
correct compiler, where `MPI_PATH` and `CUDA_PATH` are the library
paths described in [[External Libraries|Installation: External Libraries]].


<br>

## Links
* [[Makefile configuration -- Simulation Options | Installation: Simulation Options ]]
* [[Makefile configuration -- External Libraries | Installation: External Libraries ]]
* [[Back to the main page of Installation | Installation]]