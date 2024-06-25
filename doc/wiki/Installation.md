To compile GAMER, go to the source directory:

    > cd src

We recommend using the Python script
[[configure.py|Installation:-Configure.py]]
to tailor the `Makefile` for your simulation and machine. Alternatively, you can
edit `Makefile` and set the following configurations directly:
1. [[Simulation Options|Installation:-Simulation-Options]]
2. [[Compiler and Flags|Installation:-Compiler-and-Flags]]
3. [[External Libraries|Installation:-External-Libraries]]

**Caution**: on macOS, we recommend using the GNU compiler and set
[[RANDOM_NUMBER | Installation:-Simulation-Options#RANDOM_NUMBER]] to `RNG_CPP11`
in the `Makefile` (or via `--rng=RNG_CPP11` in `configure.py`).

Compile the code by

    > make clean
    > make

To reduce the compilation time, you can perform a parallel
compilation by `make -j N`, where N is the number of compilation
jobs to run in parallel. For example, the following command will
invoke 4 compilation jobs simultaneously:

    > make -j 4

However, please consult the documentation of your system to avoid
violating the usage policy.

If the compilation succeeds, you will see the following message

    > Compiling GAMER --> Successful!

and get an executable `gamer`, which will be automatically copied
to `../bin/gamer`.

**Caution**: the Makefile structure will likely be changed in the
near future to become more extensible and maintainable.