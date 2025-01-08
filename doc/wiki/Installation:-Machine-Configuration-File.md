The machine configuration file is under `configs`. The configuration file contains the library paths, the compiler types, the compilation flags, and the GPU compute capability.

If there is a configuration file that matches the name of your machine, you should set it as the default by

```bash
sh tool/config/set_settings.sh --local --machine=your_machine
```

For example, setting `--machine=pleiades` with the above command will use the `configs/pleiades.config` machine configuration when compiling the code.

If the machine configuration file is not available for your machine or the existing one is not appropriate, you will need to create a new one with the following instructions.

> [!NOTE]
> If you want to set the default machine configuration file for all of the GAMER copies under your user account, use the `--global` option instead of `--local`.
Still, you can override the global setting for the individual GAMER copies with the `--local` option.
Furthermore, you can override the default setting with the `--machine=` option in `configure.py`.

## Machine Configuration File

To set up your machine configuration file, go to `configs` and make a copy of `template.config` to modify it.
    
```bash
cd configs
cp template.config your_machine.config
```

Please refer to the following sections to set up your machine configuration file. And don't forget to set the machine configuration file as the default by the `tool/config/set_settings.sh` command above.


### 0. Rules of the configuration file
* Comments must start with `#`.
* The variable name and its value must be separated by space(s).

### 1. Library paths

For example, `MPI_PATH` can be set by:
```
MPI_PATH /usr/local/mpich-3.2
```
Please also check out the installation of [[External Libraries | Installation:-External-Libraries]].

### 2. Compilers

There are two compilers in this section: C++ (`CXX`) and MPI (`CXX_MPI`).
> [!NOTE]
> `MPI_PATH/bin/` will be combined with `CXX_MPI` automatically

### 3. Compilation flags

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

### 4. GPU compute capability

The GPU compute capability can be calculated by `major_verison*100 + minor_version*10`. For example, for `GeForce RTX 4090`, set `GPU_COMPUTE_CAPABILITY 890` (8\*100 + 9\*10).

> [!TIP]
> * You can also set `GPU_COMPUTE_CAPABILITY` to `-1` to determine the value automatically using `get_gpu_compute_capability()` in `configure.py`.
> * Check your GPU compute capability:
>   1. https://developer.nvidia.com/cuda-gpus
>   1. https://en.wikipedia.org/wiki/CUDA#GPUs_supported

<br>

## Links
* [[ Option List | Installation:-Option-List ]]
* [[External Libraries | Installation:-External-Libraries ]]
* [[Back to the main page of Installation | Installation]]
