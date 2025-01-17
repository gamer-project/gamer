* [FFTW](#FFTW)
* [GRACKLE](#GRACKLE)
* [HDF5](#HDF5)
* [LIBYT](#LIBYT)

## FFTW
GAMER supports both FFTW2 and FFTW3 for various calculations (e.g., the root-level Poisson solver).
Follow the installation instructions on the [FFTW website](http://www.fftw.org/download.html).
Note that it must be configured with
floating-point type prefix `--enable-type-prefix` for FFTW2
and MPI support `--enable-mpi` for both FFTW2 and FFTW3.
Here are example installation scripts using the GNU compiler for FFTW2 and FFTW3, respectively:

### FFTW2 (`FFTW2_PATH`) <a name="FFTW2"></a>
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

### FFTW3 (`FFTW3_PATH`) <a name="FFTW3"></a>
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

## GRACKLE (`GRACKLE_PATH`) <a name="GRACKLE"></a>
GAMER uses GRACKLE for the chemistry and radiative processes.
Follow the installation instructions in the
[GRACKLE website](http://grackle.readthedocs.io/en/latest/index.html).

Grackle and GAMER can be compiled with different floating-point precisions.
It is [recommended](https://grackle.readthedocs.io/en/latest/Installation.html#compiler-settings)
to compile Grackle in double precision:

    > make precision-64

When enabling OpenMP in GAMER (i.e., with the
compile-time option [[--openmp | Installation:-Option-List#--openmp]]),
GRACKLE must be configured with OpenMP
support as well using

    > make omp-on

## HDF5 (`HDF5_PATH`) <a name="HDF5"></a>
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

## LIBYT (`LIBYT_PATH`) <a name="LIBYT"></a>
GAMER uses [libyt](https://github.com/yt-project/libyt) for in situ Python analysis.
See [[In Situ Python Analysis | In-Situ-Python-Analysis]] for details.

libyt has two modes, normal mode and interactive mode.
Please refer to [libyt -- How to Install](https://libyt.readthedocs.io/en/latest/how-to-install/how-to-install.html#how-to-install).

Set `LIBYT_PATH` to the folder that contains subfolders `include` and `lib`.


<br>

## Links
* [[Machine Configuration File | Installation:-Machine-Configuration-File]]
* [[Option List | Installation:-Option-List]]
* [[Back to the main page of Installation | Installation]]
