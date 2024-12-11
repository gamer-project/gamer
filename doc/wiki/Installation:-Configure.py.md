This page shows how to use the Python script `configure.py` to tailor the `Makefile`
for your simulation and machine. The script supports both Python2 and Python3.
- [User Guide](#user-guide)
  - [Setup machine configuration file](#setup-machine-configuration-file)
  - [Simulation options](#simulation-options)
  - [Running the script](#running-the-script)
- [Developer Guide](#developer-guide)
  - [Adding new source files](#adding-new-source-files)
  - [Adding new library paths](#adding-new-library-paths)
  - [Adding new compiler flag types](#adding-new-compiler-flag-types)
  - [Rules of Makefile_base](#rules-of-makefile_base)

## User Guide
### Setup machine configuration file
> [!TIP]
> To set up your machine configuration file, make a copy of `template.config` and modify it.

The machine configuration file is under `gamer/configs`. The configuration file contains the library paths, the compiler types, the compilation flags, and the GPU compute capability.
1. Library paths

   For example, `MPI_PATH` can be set by:
   ```
   MPI_PATH /usr/local/mpich-3.2
   ```

2. Compilers

   There are two compilers in this section: C++ (`CXX`) and MPI (`CXX_MPI`).
> [!NOTE]
> `MPI_PATH/bin/` will be combined with `CXX_MPI` automatically

3. Compilation flags

   For example, `CXXFLAG` can be set by:
   ```
   CXXFLAG -g -O2
   ```
   or
   ```
   CXXFLAG -g
   CXXFLAG -O2
   ```
   Here is a table of all the available flag variables:
   | Flag name | Description |
   |---|---|
   | `CXXFLAG`      | Flags for compiler `CXX` and `CXX_MPI` |
   | `OPENMPFLAG`   | Flags for OpenMP |
   | `LIBFLAG`      | Flags for all libraries |
   | `NVCCFLAG_COM` | Flags for `nvcc` compiler |
   | `NVCCFLAG_FLU` | Flags for fluid solver files |
   | `NVCCFLAG_POT` | Flags for Poisson/gravity solvers files |

4. GPU compute capability

   The GPU compute capability can be calculated by `major_verison*100 + minor_version*10`. For example, for `GeForce RTX 4090`, set `GPU_COMPUTE_CAPABILITY 890` (8*100 + 9*10).
> [!TIP]
> * You can also set `GPU_COMPUTE_CAPABILITY` to `-1` to determine the value automatically using `get_gpu_compute_capability()` in `configure.py`.
> * Check your GPU compute capability:
>   1. https://developer.nvidia.com/cuda-gpus
>   1. https://en.wikipedia.org/wiki/CUDA#GPUs_supported

### Simulation options
The following commands list all available [[simulation options|Installation:-Simulation-Options]].
```bash
python configure.py -h    # show a short help message
python configure.py -lh   # show a detailed help message

```

### Running the script
Run the following command to generate a new `Makefile`.
```bash
python configure.py --machine=your_configuration_file [--your_arguments]
```

For example, the following command sets the compiler, flags, and library paths
based on `gamer/configs/pleiades.config`, uses `FFTW2`, and enables both gravity and GPU.

``` bash
python configure.py --machine=pleiades --fftw=FFTW2 --gravity=true --gpu=true
```

An example script `generate_make.sh` can be found in each test problem folder
(e.g., `example/test_problem/Hydro/AcousticWave/generate_make.sh`).

***

## Developer Guide
This script consists of five parts: `Packages`, `Global variables`, `Classes`, `Functions`, and `Main execution`.

### Adding new source files
Edit `Makefile_base` to add new source files.

### Adding new simulation options
1. Add a Python argument reader for the new simulation option under `load_arguments()`. Here is a simple example of the argument reader:
   ```python
   parser.add_argument( "--new_argument", type=int, metavar="INTEGER", gamer_name="NAME_IN_GAMER",
                        default=0,
                        help="Your help message.\n"
                      )
   ```
   * Please check out the available options at [argparse document](https://docs.python.org/3/library/argparse.html#quick-links-for-add-argument).
   * `gamer_name` is the simulation option name in GAMER.
   * `default` is the default value of the argument. If the argument default depends on other arguments
     (e.g., the default of `bitwise_reproducibility` is `True` when enabling `--debug` but otherwise is `False`),
     set `default=None` and assign the default value under `set_conditional_defaults()`.
     ```python
     def set_conditional_defaults( args ):
         ...
         if args["new_argument"] == None:
             args["new_argument"] = default_value_of_true if args["other_argument"] else default_value_of_false
         ...
         return args
     ```

2. [Optional] If the argument depends on other arguments, add `depend={"depend_arg1":depend_value1, "depend_arg2":depend_value2}` so the argument will be loaded only if `depend_arg1==depend_value1` and `depend_arg2==depend_value2`.
   ```python
   parser.add_argument( "--new_argument", type=int, metavar="INTEGER", gamer_name="NEW_SIMUALTION_OPTION",
                        default=0,
                        depend={"depend_arg1":depend_value1, "depend_arg2":depend_value2},
                        help="Your help message.\n"
                      )
   ```

3. [Optional] To validate the input values, add `constraint={ val1:{"arg1":["a", "b"], val2:{"arg2":"c"} }`,
   which will check whether the argument `arg1` is either `a` or `b` when the input value is `val1`
   and whether the argument `arg2` is `c` when the input value is `val2`.
   An error will be raised if any constraints are violated. For example, the following code
   asserts `--eos=GAMMA` when adopting either `--flux=ROE` or `--flux=EXACT`.
   ```python
   parser.add_argument( "--flux", type=str, metavar="TYPE", gamer_name="RSOLVER",
                        choices=["EXACT", "ROE", "HLLE", "HLLC", "HLLD"],
                        constraint={ "ROE":{"eos":"GAMMA"},
                                     "EXACT":{"eos":"GAMMA"} },
                        ...
                      )
   ```

4. [Optional] Add additional checks in `validation()` and warning messages in `warning()` under `Functions`.
   * `validation()`
   ```python
   def validation( paths, depends, constraints, **kwargs ):
       success = True
       ...

       if kwargs["new_argument"] < -1:
           LOGGER.error("Your error message.")
           success = False
       ...

       if not success: raise BaseException("The above validation failed.")
       return
   ```
   * `warning()`
   ```python
   def warning( paths, **kwargs ):
       ...
       if kwargs["new_argument"] == 0:
           LOGGER.warning("Your warning message.")
       ...
       return
   ```

### Adding new library paths
1. Add `NEW_PATH := @@@NEW_PATH@@@` in `Makefile_base` under `# library paths`.
   ```makefile
   # library paths
   #######################################################################################################
   ... other paths ...
   NEW_PATH     := @@@NEW_PATH@@@
   ...
   ```
2. Add `NEW_PATH /path/of/new` in your machine configuration file `configs/YOUR.config`.
   ```
   # 1. Paths
   ... other paths ...
   NEW_PATH        /path/of/new
   ...
   ```

### Adding new compiler flag types
1. Add `NEW_FLAG := @@@NEW_FLAG@@@` in `Makefile_base` under `# compilers and flags`.
   ```makefile
   # compilers and flags
   #######################################################################################################
   ... other flags ...
   NEW_FLAG := @@@NEW_FLAG@@@
   ...
   ```

2. Add `["NEW_FLAG":""]` in the dictionary variable `flags` of `load_config()` in `configure.py`.
   ```python
   def load_config( config ):
       LOGGER.info("Using %s as the config."%(config))
       paths, compilers = {}, {"CXX":"", "CXX_MPI":""}
       flags = {"CXXFLAG":"", "OPENMPFLAG":"", "LIBFLAG":"", "NVCCFLAG_COM":"", "NVCCFLAG_FLU":"", "NVCCFLAG_POT":""}
       gpus  = {"GPU_COMPUTE_CAPABILITY":""}
       ...
       return paths, compilers, flags
   ```
> [!IMPORTANT]
> All flags must be set in `flags`; otherwise, they will be interpreted as library paths.

3. Add `NEW_FLAG -new_flag` in your machine configuration file `configs/YOUR.config`.
   ```
   # 2. Compiler flags
   ... other flags ...
   NEW_FLAG    -new_flag
   ...
   ```

### Rules of `Makefile_base`
* The strings to be replaced by `configure.py` must be sandwiched by `@@@`.

### Rules of `*.config`
* Comments must start with `#`.
* All variables should be uppercase.
* Variables and values should be separated by spaces.
