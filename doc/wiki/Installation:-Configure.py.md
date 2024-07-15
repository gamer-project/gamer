This page shows how to use the Python script `configure.py` to tailor the `Makefile`
for your simulation and machine. The script supports both Python2 and Python3.
- [User Guide](#user-guide)
  - [Simulation options](#simulation-options)
  - [Library paths and compilation flags](#library-paths-and-compilation-flags)
  - [Running the script](#running-the-script)
- [Developer Guide](#developer-guide)
  - [Adding new source files](#adding-new-source-files)
  - [Adding new library paths](#adding-new-library-paths)
  - [Adding new compiler flag types](#adding-new-compiler-flag-types)
  - [Rules of Makefile_base](#rules-of-makefile_base)

## User Guide
### Simulation options
The following commands list all available [[simulation options|Installation:-Simulation-Options]].
```bash
python configure.py -h    # show a short help message
python configure.py -lh   # show a detailed help message

```

### Library paths and compilation flags
Edit the machine configuration file in `gamer/configs` to specify the library paths and compilation flags.
To set up your own machine configuration file, make a copy of `template.config` and modify it.

### Running the script
Run the following command to generate a new `Makefile`.
```bash
python configure.py --machine=your_configuration_file [--your_arguments]
```

For example, the following command will set the compiler, flags, and library paths
based on `gamer/configs/pleiades.config`, use `FFTW2`, and enable both gravity and GPU.

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
Here is a simple example of the argument reader:
```python
parser.add_argument( "--new_argument", type=int, metavar="INTEGER", gamer_name="NAME_IN_GAMER",
                     default=0,
                     help="Your help message.\n"
                   )
```
1. Add a Python argument reader for the new simulation option under `load_arguments()`.
2. If the argument default depends on other arguments
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

3. [Optional] If the argument depends on other arguments, add `depend={"depend_arg1":depend_value1, "depend_arg2":depend_value2}` so the argument will be loaded only if
`depend_arg1==depend_value1` and `depend_arg2==depend_value2`.
```python
parser.add_argument( "--new_argument", type=int, metavar="INTEGER", gamer_name="NEW_SIMUALTION_OPTION",
                     default=0,
                     depend={"depend_arg1":depend_value1, "depend_arg2":depend_value2},
                     help="Your help message.\n"
                   )
```

4. [Optional] To validate the input values, add `constraint={ val1:{"arg1":["a", "b"], val2:{"arg2":"c"} }`,
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

5. [Optional] Add additional checks in `validation()` and warning messages in `warning()` under `Functions`.
* `validation()`
```python
def validation( paths, depends, constraints, **kwargs ):
    success = True
    ...

    if kwargs["new_argument"] < -1:
        color_print("ERROR: Your error message.", BCOLOR.FAIL)
        success = False
    ...

    if not success: raise BaseException(BCOLOR.FAIL+"The above validation failed."+BCOLOR.ENDC)
    return
```
* `warning()`
```python
def warning( paths, **kwargs ):
    ...
    if kwargs["new_argument"] == 0:
        color_print("Warning: Your warning message.", BCOLOR.WARNING)
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
    print("Using %s as the config."%(config))
    paths, compilers, flags = {}, {"CXX":"", "CXX_MPI":""}, {"CXXFLAG":"", "OPENMPFLAG":"", "LIBFLAG":"", "CUDAFLAG":"", "NEW_FLAG":""}
    ...
    return paths, compilers, flags
```

3. Add `NEW_FLAG -new_flag` in your machine configuration file `configs/YOUR.config`.
```
# 2. Compiler flags
... other flags ...
NEW_FLAG    -new_flag
...
```

### Rules of `Makefile_base`
The strings to be replaced by `configure.py` must be sandwiched by `@@@`.