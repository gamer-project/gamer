- [Developer Guide](#developer-guide)
  - [Adding new source files](#adding-new-source-files)
  - [Adding new library paths](#adding-new-library-paths)
  - [Adding new compiler flag types](#adding-new-compiler-flag-types)
  - [Rules of Makefile_base](#rules-of-makefile_base)
# Developer Guide
This script consists of five parts: `Packages`, `Global variables`, `Classes`, `Functions`, and `Main execution`.

## Adding new source files
Edit `Makefile_base` to add new source files.

## Adding new simulation options
1. Add a Python argument reader for the new simulation option under `load_arguments()`. Here is a simple example of the argument reader:

   ```python
   parser.add_argument( "--new_argument", type=int, metavar="INTEGER", gamer_name="NAME_IN_GAMER",
                        default=0,
                        help="Your help message.\n"
                      )
   ```

   * Please check out the available options at [argparse document](https://docs.python.org/3/library/argparse.html#quick-links-for-add-argument).
   * `gamer_name` is the simulation option name in GAMER.
   * `default` is the default value of the argument. If the argument default depends on other arguments,
     set `default=None` and assign the default value under `set_conditional_defaults()`.
     For example, the default of `bitwise_reproducibility` is `True` when enabling `--debug` but otherwise is `False`.

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

## Adding new library paths
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

## Adding new compiler flag types
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

## Rules of `*.config`
* Comments must start with `#`.
* All variables should be uppercase.
* Variables and values should be separated by spaces.

## Rules of `Makefile_base`
* The strings to be replaced by `configure.py` must be sandwiched by `@@@`.

### Format
All compile-time simulation options in the `Makefile` are in the following two formats:

```Makefile
SIMU_OPTION += -DOPTION1
SIMU_OPTION += -DOPTION2=OPTION2_ADOPTED
```
which will enable `OPTION1` and assign `OPTION2_ADOPTED` to `OPTION2`.
For example, to (i) enable gravity and (ii) adopt the CTU fluid scheme, set

```Makefile
SIMU_OPTION += -DGRAVITY
SIMU_OPTION += -DFLU_SCHEME=CTU
```

> [!CAUTION]
> * Option values (if any) must be set explicitly since there are no default values.
> For example, `SIMU_OPTION += -DFLU_SCHEME` without assigning any value to the option `FLU_SCHEME` is invalid.
> * Do not insert any space before and after the equal sign `=`.
> For example, use `-DFLU_SCHEME=CTU` instead of `-DFLU_SCHEME = CTU`.


### Compilers and flags
To choose a compiler and compilation flags, set the following variables in the `Makefile`:

```Makefile
CXX        = @@@CXX@@@              # C++ compiler
CXXFLAG    = @@@CXXFLAG@@@          # compilation flags
LIB        = @@@LIBFLAG@@@          # openmp flag
OPENMPFLAG = @@@OPENMPFLAG@@@       # libraries and linker flags
NVCC       = $(CUDA_PATH)/bin/nvcc  # CUDA compiler
```

The prefixes `$(MPI_PATH)/bin/` and `$(CUDA_PATH)/bin/` in the
above examples are optional to force the Makefile to use the
correct compiler, where `MPI_PATH` and `CUDA_PATH` are the library
paths described in [[External Libraries | Installation: External Libraries]].

### Library Paths
Set the following library paths in the `Makefile` to help the compiler locate them (if necessary):

``` Makefile
CUDA_PATH    := @@@CUDA_PATH@@@
FFTW2_PATH   := @@@FFTW2_PATH@@@
FFTW3_PATH   := @@@FFTW3_PATH@@@
MPI_PATH     := @@@MPI_PATH@@@
HDF5_PATH    := @@@HDF5_PATH@@@
GRACKLE_PATH := @@@GRACKLE_PATH@@@
GSL_PATH     := @@@GSL_PATH@@@
LIBYT_PATH   := @@@LIBYT_PATH@@@
```

Only the paths of libraries being used need to be set. In addition,
it is usually unnecessary to set the paths that have been embedded
into the compiling command (e.g., when using `CC` and `module load`
in a Cray computer system).