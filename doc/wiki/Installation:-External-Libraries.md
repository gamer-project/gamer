## Library Paths

Set the following library paths in the `Makefile` to help
compiler locate them (if necessary):

``` Makefile
CUDA_PATH    :=
FFTW2_PATH   :=
FFTW3_PATH   :=
MPI_PATH     :=
HDF5_PATH    :=
GRACKLE_PATH :=
GSL_PATH     :=
LIBYT_PATH   :=
```
Only the paths of libraries being used need to be set. In addition,
it is usually unnecessary to set the paths that have been embedded
into the compiling command (e.g., when using `CC` and `module load`
in a Cray computer system).

## Library Configurations

### FFTW
GAMER supports both FFTW2 and FFTW3 for various calculations (e.g., the root-level Poisson solver).
Follow the installation instructions on the [FFTW website](http://www.fftw.org/download.html).
Note that it must be configured with
floating-point type prefix `--enable-type-prefix` for FFTW2
and MPI support `--enable-mpi` for both FFTW2 and FFTW3.
Here are example installation scripts using the GNU compiler for FFTW2 and FFTW3, respectively:

#### FFTW2
``` bash
export FFTW_PATH=PATH_TO_INSTALL_YOUR_FFTW
export CC=gcc
export F77=gfortran
export CFLAGS=-O3
export FFLAGS=-O3

make clean
./configure --enable-mpi --enable-type-prefix --prefix $FFTW_PATH --enable-float
make
make install

make clean
./configure --enable-mpi --enable-type-prefix --prefix $FFTW_PATH
make
make install
```

#### FFTW3
``` bash
export FFTW_PATH=PATH_TO_INSTALL_YOUR_FFTW
export CC=gcc
export F77=gfortran
export CFLAGS=-O3
export FFLAGS=-O3

make distclean
./configure --enable-openmp --enable-mpi --enable-shared --prefix $FFTW_PATH --enable-float
make
make install

make distclean
./configure --enable-openmp --enable-mpi --enable-shared --prefix $FFTW_PATH
make
make install

```

### GRACKLE
GAMER uses GRACKLE for the chemistry and radiative processes.
Follow the installation instructions in the
[GRACKLE website](http://grackle.readthedocs.io/en/latest/index.html).
Note that it must be configured with a
consistent floating-point accuracy as GAMER using
```
make precision-{32,64}
```
Specifically, configure GRACKLE with `make precision-64/32` when
compiling GAMER with/without the option
[[FLOAT8 | Installation: Simulation-Options#FLOAT8]], respectively.

In addition, when enabling OpenMP in GAMER (i.e., with the
compile-time option [[OPENMP | Installation: Simulation-Options#OPENMP]]),
GRACKLE must be configured with OpenMP
support as well using
```
make omp-on
```
### HDF5
GAMER uses [HDF5](https://support.hdfgroup.org/HDF5/) for storing snapshots.
It is not necessary to enable either `--enable-cxx` or `--enable-parallel` when
configuring HDF5 since GAMER currently adopts the C interface with serial I/O.
Here is an example installation script using the `icc` compiler
(remember to replace `PATH_TO_INSTALL_YOUR_HDF5` with your installation path):

``` bash
CC=icc ./configure --prefix=PATH_TO_INSTALL_YOUR_HDF5
make
make install
```

### LIBYT
GAMER uses [libyt](https://github.com/yt-project/libyt) for in situ Python analysis.
See [[In Situ Python Analysis | In-Situ-Python-Analysis]] for details.

libyt has two modes, normal mode and interactive mode.
Please refer to [libyt -- How to Install](https://yt-project.github.io/libyt/HowToInstall.html#libyt).

Set `LIBYT_PATH` to the folder that contains subfolders `include` and `lib`.


<br>

## Links
* [[Makefile configuration -- Simulation Options | Installation: Simulation Options]]
* [[Makefile configuration -- Compiler and Flags | Installation: Compiler and Flags]]
* [[Back to the main page of Installation | Installation]]
