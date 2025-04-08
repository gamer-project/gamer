This page includes three demos:
1. [CPU-only without OpenMP](#cpu-only-without-openmp)
2. [CPU-only with OpenMP](#cpu-only-with-openmp)
3. [Hybrid OpenMP/GPU](#hybrid-openmpgpu)

## CPU-only without OpenMP

1\. Move to the source directory.
``` bash
cd src
```

2\. Generate `Makefile` by [[configure.py | Installation]].
``` bash
cp ../example/test_problem/Hydro/Riemann/generate_make.sh ./
sh generate_make.sh --openmp=false
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
   ...
   ...
========================================
Makefile is created.
========================================
</pre>
</details>

> [!NOTE]
> We have set
[[--mpi | Installation:-Option-List#--mpi]]=false,
[[--gpu | Installation:-Option-List#--gpu]]=false, and
[[--openmp | Installation:-Option-List#--openmp]]=false
to run in a CPU-only mode
without OpenMP and MPI.
See [[Option List | Installation:-Option-List#Option-List]]
for a complete list of all available options of `configure.py`.

3\. Compile the code.
``` bash
make clean
make
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
   ...
   ...
Compiling GAMER --> Successful!
</pre>
</details>

4\. Create a working directory.
``` bash
cd ../bin
mkdir shocktube
cd shocktube
```

5\. Copy the GAMER executable and example files of the
test problem to the working directory.
```bash
cp -r ../../example/test_problem/Hydro/Riemann/* .
cp ../gamer .
ls
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
clean.sh  gamer  Input__Flag_Lohner  Input__Parameter  Input__TestProb  plot__hydro_dens.gpt  plot__mhd.gpt  README  ReferenceSolution
</pre>
</details>

6\. Run the code. It will display `~ GAMER OVER ~` if it succeeds.
```bash
./gamer
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
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
</pre>
</details>

7\. The code will generate several [[log files|Outputs#log-files]] `Record__*`
and a series of 1D text data files `Xline_y0.000_z0.000_*`.

``` bash
ls
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
clean.sh              plot__mhd.gpt       Record__PatchCount          Xline_y0.000_z0.000_000001  Xline_y0.000_z0.000_000007
gamer                 README              Record__Performance         Xline_y0.000_z0.000_000002  Xline_y0.000_z0.000_000008
Input__Flag_Lohner    Record__Dump        Record__TimeStep            Xline_y0.000_z0.000_000003  Xline_y0.000_z0.000_000009
Input__Parameter      Record__MemInfo     Record__Timing              Xline_y0.000_z0.000_000004  Xline_y0.000_z0.000_000010
Input__TestProb       Record__NCorrUnphy  ReferenceSolution           Xline_y0.000_z0.000_000005
plot__hydro_dens.gpt  Record__Note        Xline_y0.000_z0.000_000000  Xline_y0.000_z0.000_000006
</pre>
</details>

8\. Plot a 1D text data file. You can use the sample [gnuplot](http://www.gnuplot.info)
script `plot__hydro_dens.gpt` (replace `display` by other image viewers if necessary).
``` bash
gnuplot plot__hydro_dens.gpt
display Fig__Riemann_Density_000010.png
```
<details>
<summary><u><i>Execution results</i></u></summary>

[[/images/shocktube.png | alt=shocktube]]
</details>


9\. Check the performance. See [[Log Files | Outputs#log-files]] for more detailed
timing analysis.
``` bash
tail -n 3 Record__Note
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Total Processing Time : 75.954923 s
</pre>
</details>

## CPU-only with OpenMP

Next, we enable OpenMP for the same test problem.
Repeat the steps above with the following modifications.

1\. Re-generate `Makefile` by [[configure.py | Installation]] and recompile `gamer`.
``` bash
sh generate_make.sh --openmp=true
make clean
make -j4
```

> [!CAUTION]
> Remember to copy the new executable to `bin/shocktube`.

2\. Set the number of OpenMP threads by editing the runtime parameter
[[OMP_NTHREAD | Runtime-Parameters:-MPI-and-OpenMP#OMP_NTHREAD]]
in the input file
[[Input__Parameter | Runtime-Parameters:-Input__Parameter]].
The following example uses 4 threads.
```
OMP_NTHREAD      4      # number of OpenMP threads (<=0=auto) [-1]
```
> [!CAUTION]
> All input files `Input__*` must be put in the same directory as the executable `gamer`.

3\. Remove all old log and data files, if necessary.
You can use the helper script `clean.sh`.
```bash
sh clean.sh
ls
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
clean.sh  gamer  Input__Flag_Lohner  Input__Parameter  Input__TestProb  plot__hydro_dens.gpt  plot__mhd.gpt  README  ReferenceSolution
</pre>
</details>

4\. Run the code in the same way as without OpenMP.
```bash
./gamer
```

5\. Validate the OpenMP settings by searching for the keyword "OpenMP"
in the log file
[[Record__Note | Simulation-Logs:-Record__Note]].
You should see something like
<pre>
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
</pre>
Check the following things:
* **The number under `NThread` is the same as the runtime
parameter
[[OMP_NTHREAD | Runtime-Parameters:-MPI-and-OpenMP#OMP_NTHREAD]]
you just set**
* **Different threads use different CPU cores**

6\. Check the performance. It should be about
[[OMP_NTHREAD | Runtime-Parameters:-MPI-and-OpenMP#OMP_NTHREAD]]
times faster
than the case without OpenMP.
``` bash
tail -n 3 Record__Note
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Total Processing Time : 20.460586 s
</pre>
</details>

## Hybrid OpenMP/GPU

To enable both GPU and OpenMP, repeat the steps in
[CPU-only with OpenMP](#cpu-only-with-openmp) with the
following modifications.

1\. Re-generate `Makefile` by [[configure.py | Installation]] and recompile `gamer`.
``` bash
sh generate_make.sh --openmp=true --gpu=true
make clean
make -j4
```

> [!CAUTION]
> * Please make sure that the `GPU_COMPUTE_CAPABILITY` is set properly in your [[machine configuration file | Installation:-Machine-Configuration-File]]
> * Remember to copy the new executable to `bin/shocktube`.

2\. Remove all old log and data files, if necessary.
Run the code with the new executable.
```bash
sh clean.sh
./gamer
```

3\. Validate the GPU settings by searching for the keyword "Device Diagnosis"
in the log file
[[Record__Note | Simulation-Logs:-Record__Note]].
You should see something like
<pre>
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
</pre>
This example shows that we are running on the computing node `golub123`
with 2 GPUs, and we are using the one with `GPU ID = 0`,
a `Tesla K40m` GPU.

4\. Check the performance. It should be noticeably faster
than the case without GPU.
``` bash
tail -n 3 Record__Note
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Total Processing Time : 9.417532 s
</pre>
</details>


<br>

## Links
* [[Next demo -- 3D Blast Wave: hybrid MPI/OpenMP/GPU + yt analysis | Quick-Start:-3D-Blast-Wave]]
* [[Back to the main page of Quick Start | Quick-Start]]

