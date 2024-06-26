This page includes three demos:
1. [CPU-only without OpenMP](#cpu-only-without-openmp)
2. [CPU-only with OpenMP](#cpu-only-with-openmp)
3. [Hybrid OpenMP/GPU](#hybrid-openmpgpu)

## CPU-only without OpenMP

1\. Move to the source directory.
``` bash
> cd src
```

2\. **[Deprecated; use `configure.py` described below instead]**: Edit the `Makefile` to validate the following settings.
Note that we have enabled
[[SERIAL | Installation:-Simulation-Options#SERIAL]]
and disabled
[[GPU | Installation:-Simulation-Options#GPU]],
[[LOAD_BALANCE | Installation:-Simulation-Options#LOAD_BALANCE]], and
[[OPENMP | Installation:-Simulation-Options#OPENMP]]
to run in a CPU-only mode
without OpenMP and MPI.
See [[Installation: Simulation Options|Installation:-Simulation-Options]]
for a complete list of all compile-time simulation options.
``` makefile
SIMU_OPTION += -DMODEL=HYDRO
#SIMU_OPTION += -DGRAVITY
#SIMU_OPTION += -DPARTICLE
#SIMU_OPTION += -DSUPPORT_GRACKLE

#SIMU_OPTION += -DGPU
SIMU_OPTION += -DSERIAL
#SIMU_OPTION += -DLOAD_BALANCE=HILBERT
#SIMU_OPTION += -DOPENMP
```

**Update**: You can use the Python script
[configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py)
to tailor `Makefile` rather than editing it manually. An example can be found at
`../example/test_problem/Hydro/Riemann/generate_make.sh`, for which you can disable
OpenMP by adding `--openmp=false`.

3\. Compile the code.
``` bash
> make clean
> make -j
   ...
   ...
Compiling GAMER --> Successful!
```

4\. Create a working directory.
``` bash
> cd ../bin
> mkdir shocktube
> cd shocktube
```

5\. Copy the GAMER executable and example files of the
test problem to the working directory.
```bash
> cp -r ../../example/test_problem/Hydro/Riemann/* .
> cp ../gamer .
> ls
clean.sh  gamer  Input__Flag_Lohner  Input__Parameter  Input__TestProb  plot__hydro_dens.gpt  plot__mhd.gpt  README  ReferenceSolution
```

6\. Run the code. It will display `~ GAMER OVER ~` if succeeds.
```bash
> ./gamer
   ...
   ...
Time: 9.0000000e-02 -> 9.3057198e-02,   Step:      29 ->      30,   dt_base: 3.0571979e-03
Time: 9.3057198e-02 -> 9.6446514e-02,   Step:      30 ->      31,   dt_base: 3.3893166e-03
Time: 9.6446514e-02 -> 9.9834923e-02,   Step:      31 ->      32,   dt_base: 3.3884083e-03
Time: 9.9834923e-02 -> 1.0000000e-01,   Step:      32 ->      33,   dt_base: 1.6507717e-04
Output_DumpData_Part (DumpID = 10) ...
Output_DumpData_Part (DumpID = 10) ... done
End_GAMER ...
End_MemFree ... done
End_GAMER ... done


~ GAME OVER ~
```

7\. The code should have generated several [[log files|Outputs#log-files]] `Record__*`
and a series of 1D text data files `Xline_y0.000_z0.000_*`.

``` bash
> ls
clean.sh              plot__mhd.gpt       Record__PatchCount          Xline_y0.000_z0.000_000001  Xline_y0.000_z0.000_000007
gamer                 README              Record__Performance         Xline_y0.000_z0.000_000002  Xline_y0.000_z0.000_000008
Input__Flag_Lohner    Record__Dump        Record__TimeStep            Xline_y0.000_z0.000_000003  Xline_y0.000_z0.000_000009
Input__Parameter      Record__MemInfo     Record__Timing              Xline_y0.000_z0.000_000004  Xline_y0.000_z0.000_000010
Input__TestProb       Record__NCorrUnphy  ReferenceSolution           Xline_y0.000_z0.000_000005
plot__hydro_dens.gpt  Record__Note        Xline_y0.000_z0.000_000000  Xline_y0.000_z0.000_000006
```

8\. Plot a 1D text data file. You can use the sample [gnuplot](http://www.gnuplot.info)
script `plot__hydro_dens.gpt` (replace `display` by other image viewers if necessary).
``` bash
> gnuplot plot__hydro_dens.gpt
> display Fig__Riemann_Density_000010.png
```
[[images/shocktube.png | alt=shocktube]]

9\. Check the performance. See [[Log Files | Outputs#log-files]] for more detailed
timing analysis.
``` bash
> tail -n 3 Record__Note

Total Processing Time : 75.954923 s

```

## CPU-only with OpenMP

Next, we enable OpenMP for the same test problem.
Repeat the steps above with the following modifications.

1\. **[Deprecated]**: Enable
[[OPENMP | Installation:-Simulation-Options#OPENMP]]
in the `Makefile` and recompile the code.
``` makefile
SIMU_OPTION += -DOPENMP
```
**Update**: If you are using
[configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py),
you can enable OpenMP by adding `--openmp=true`.

**Caution: remember to clean the previous configurations by `make clean`
before `make` and copy the new executable to the working directory.**

2\. Set the number of OpenMP threads by editing the runtime parameter
[[OMP_NTHREAD | MPI-and-OpenMP#OMP_NTHREAD]]
in the input file
[[Input__Parameter | Runtime-Parameters#input__parameter]].
The following example uses 4 threads.
```
OMP_NTHREAD      4      # number of OpenMP threads (<=0=auto) [-1]
```
**Caution: all input files `Input__*` must be put in the
same directory as the executable `gamer`.**

3\. Remove all old log and data files, if necessary.
You can use the helper script `clean.sh`.
```bash
> sh clean.sh
> ls
clean.sh  gamer  Input__Flag_Lohner  Input__Parameter  Input__TestProb  plot__hydro_dens.gpt  plot__mhd.gpt  README  ReferenceSolution
```

4\. Run the code in the same way as without OpenMP.
```bash
> ./gamer
```

5\. Validate the OpenMP settings by searching for the keyword "OpenMP"
in the log file
[[Record__Note | Simulation-Logs:-Record__Note]].
You should see something like
```
OpenMP Diagnosis
***********************************************************************************
OMP__SCHEDULE                   DYNAMIC
OMP__SCHEDULE_CHUNK_SIZE        1
OMP__NESTED                     OFF

CPU core IDs of all OpenMP threads (tid == thread ID):
------------------------------------------------------------------------
 Rank        Host  NThread  tid-00  tid-01  tid-02  tid-03
    0    golub123        4       2       5       4       7
***********************************************************************************
```
Check the following things:
* **The number under `NThread` is the same as the runtime
parameter
[[OMP_NTHREAD | MPI-and-OpenMP#OMP_NTHREAD]]
you just set**
* **Different threads use different CPU cores**

6\. Check the performance. It should be about
[[OMP_NTHREAD | MPI-and-OpenMP#OMP_NTHREAD]]
times faster
than the case without OpenMP.
``` bash
> tail -n 3 Record__Note

Total Processing Time : 20.460586 s

```

## Hybrid OpenMP/GPU

To enable both GPU and OpenMP, repeat the steps in
[CPU-only with OpenMP](#cpu-only-with-openmp) with the
following modifications.

1\.  **[Deprecated]**: Edit the `Makefile` as described below and recompile the code.
* Enable [[GPU | Installation:-Simulation-Options#GPU]].
* Set [[GPU_COMPUTE_CAPABILITY | Installation:-Simulation-Options#GPU_COMPUTE_CAPABILITY]]
to match the GPU on your system (e.g., `GPU_COMPUTE_CAPABILITY=890` for GeForce RTX 4090).
* Set `CUDA_PATH` to the location of your CUDA installation.

The following example assumes that the `TURING` architecture
is adopted and CUDA is installed at `/usr/local/cuda`.
``` makefile
SIMU_OPTION += -DGPU
SIMU_OPTION += -DDGPU_COMPUTE_CAPABILITY=860

CUDA_PATH   := /usr/local/cuda
```

**Update**: If you are using
[configure.py](https://github.com/gamer-project/gamer/wiki/Installation%3A-Configure.py),
you can first set `CUDA_PATH` and `GPU_COMPUTE_CAPABILITY` in your machine configuration file in `gamer/config/YOUR_MACHINE.config`
and then add `--machine=YOUR_MACHINE --openmp=true --gpu=true` when running `configure.py`.

2\. Remove all old log and data files, if necessary.
Run the code with the new executable.
```bash
> sh clean.sh
> ./gamer
```

3\. Validate the GPU settings by searching for the keyword "Device Diagnosis"
in the log file
[[Record__Note | Simulation-Logs:-Record__Note]].
You should see something like
```
Device Diagnosis
***********************************************************************************
MPI_Rank =   0, hostname =   golub123, PID = 47842

CPU Info :
CPU Type        : Intel(R) Xeon(R) CPU E5-2670 v2 @ 2.50GHz
CPU MHz         : 2499.982
Cache Size      : 25600 KB
CPU Cores       : 10
Total Memory    : 63.0 GB

GPU Info :
Number of GPUs                    : 2
GPU ID                            : 0
GPU Name                          : Tesla K40m
CUDA Driver Version               : 8.0
CUDA Runtime Version              : 7.0
CUDA Major Revision Number        : 3
CUDA Minor Revision Number        : 5
Clock Rate                        : 0.745000 GHz
Global Memory Size                : 11439 MB
Constant Memory Size              : 64 KB
Shared Memory Size per Block      : 48 KB
Number of Registers per Block     : 65536
Warp Size                         : 32
Number of Multiprocessors:        : 15
Number of Cores per Multiprocessor: 192
Total Number of Cores:            : 2880
Max Number of Threads per Block   : 1024
Max Size of the Block X-Dimension : 1024
Max Size of the Grid X-Dimension  : 2147483647
Concurrent Copy and Execution     : Yes
Concurrent Up/Downstream Copies   : Yes
Concurrent Kernel Execution       : Yes
GPU has ECC Support Enabled       : Yes
***********************************************************************************
```
This example shows that we are running on the computing node "golub123"
with 2 GPUs, and we are using the one with "GPU ID = 0",
a "Tesla K40m" GPU.

4\. Check the performance. It should be noticeably faster
than the case without GPU.
``` bash
> tail -n 3 Record__Note

Total Processing Time : 9.417532 s

```


<br>

## Links
* [[Next demo -- 3D Blast Wave: hybrid MPI/OpenMP/GPU + yt analysis | Quick-Start:-3D-Blast-Wave]]
* [[Back to the main page of Quick Start | Quick-Start]]

