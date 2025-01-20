The machine configuration file is located under `configs`. It specifies the library paths, compiler types, compilation flags, and GPU compute capability.

> [!TIP]
> Check the `configs` directory to see if a configuration file is already available for your machine.

## Set Up the Machine Configuration File

To set up your machine configuration file, go to `configs` and make a copy of `template.config` to modify:
    
```bash
cd configs
cp template.config your_machine.config
```

Please refer to the following sections to set up your machine configuration file.

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
