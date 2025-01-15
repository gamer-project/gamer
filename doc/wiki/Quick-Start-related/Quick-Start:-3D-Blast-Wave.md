This demo will illustrate the following features:
* Hybrid MPI/OpenMP/GPU
* HDF5 output
* yt analysis

It is recommended to first check [[Quick Start: 1D Shock Tube]].

***

1\. Install the external packages required for this demo.
* [yt](http://yt-project.org)
* [HDF5](https://support.hdfgroup.org/HDF5)

2\. Generate `Makefile` by [[configure.py | Installation]].
``` bash
cp ../example/test_problem/Hydro/BlastWave/generate_make.sh ./
sh generate_make.sh --openmp=true --gpu=true --mpi=true --hdf5=true
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

3\. Compile the code.
``` bash
make clean
make -j
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
   ...
   ...
Compiling GAMER --> Successful!
</pre>
</details>

4\. Create a working directory and copy the GAMER executable and
example files of the test problem.
``` bash
cd ../bin
mkdir blastwave
cd blastwave
cp -r ../../example/test_problem/Hydro/BlastWave/* .
cp ../gamer .
ls
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Input__Flag_Lohner  Input__Parameter  Input__TestProb  README  clean.sh  gamer  plot_profile.gpt  plot_slice.py
</pre>
</details>

5\. Edit the input file
[[Input__Parameter | Runtime-Parameters:-Input__Parameter]]
to set the total number of OpenMP threads per MPI process.
The following example will launch 10 threads per process.
```
OMP_NTHREAD      10      # number of OpenMP threads (<=0=auto) [-1]
```

> [!IMPORTANT]
> See [[Hybrid MPI/OpenMP/GPU | MPI-and-OpenMP#hybrid-mpiopenmpgpu]]
for the recommended configuration of the number of MPI processes and OpenMP threads.

6\. Run the code with MPI (e.g., `mpirun, mpiexec, aprun`).
> [!IMPORTANT]
> Please consult the documentation of your system, especially about how to compile and launch "hybrid MPI/OpenMP" jobs.

The following example uses OpenMPI 1.8.4 to launch 2 MPI processes
and 10 threads per process, assuming that there are two 10-core CPUs.
```bash
mpirun -np 2 -map-by ppr:1:socket:pe=10 ./gamer
```
> [!TIP]
> See more options of MPI:
>   * `-map-by`: [Mapping options and modifiers](https://www.ibm.com/docs/en/smpi/10.2?topic=affinity-mapping-options-modifiers)
>   * `-np`: [mpirun command options](https://www.ibm.com/docs/da/smpi/10.2?topic=command-mpirun-options)
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
   ...
   ...
Time: 4.8438672e-03 -> 4.8932853e-03,   Step:     114 ->     115,   dt_base: 4.9418079e-05
Time: 4.8932853e-03 -> 4.9427910e-03,   Step:     115 ->     116,   dt_base: 4.9505779e-05
Time: 4.9427910e-03 -> 4.9923843e-03,   Step:     116 ->     117,   dt_base: 4.9593262e-05
Time: 4.9923843e-03 -> 5.0000000e-03,   Step:     117 ->     118,   dt_base: 7.6156993e-06
Output_DumpData_Total_HDF5 (DumpID = 10) ...
Output_DumpData_Total_HDF5 (DumpID = 10) ... done
Output_DumpData_Part (DumpID = 10) ...
Output_DumpData_Part (DumpID = 10) ... done
End_GAMER ...
End_MemFree ... done
End_GAMER ... done


~ GAME OVER ~
</pre>
</details>

7\. The code should have generated several [[log files|Outputs#log-files]]
`Record__*`, a series of 1D text data files `Diag_*`, and a series of
HDF5 snapshots `Data_*`.
```bash
ls
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Data_000000  Data_000005  Data_000010  Diag_000004  Diag_000009         README               Record__Note         clean.sh
Data_000001  Data_000006  Diag_000000  Diag_000005  Diag_000010         Record__Dump         Record__PatchCount   gamer
Data_000002  Data_000007  Diag_000001  Diag_000006  Input__Flag_Lohner  Record__LoadBalance  Record__Performance  plot_profile.gpt
Data_000003  Data_000008  Diag_000002  Diag_000007  Input__Parameter    Record__MemInfo      Record__TimeStep     plot_slice.py
Data_000004  Data_000009  Diag_000003  Diag_000008  Input__TestProb     Record__NCorrUnphy   Record__Timing
</pre>
</details>

```bash
h5ls -r Data_000010
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
/                        Group
/GridData                Group
/GridData/Dens           Dataset {6016, 8, 8, 8}
/GridData/Engy           Dataset {6016, 8, 8, 8}
/GridData/MomX           Dataset {6016, 8, 8, 8}
/GridData/MomY           Dataset {6016, 8, 8, 8}
/GridData/MomZ           Dataset {6016, 8, 8, 8}
/Info                    Group
/Info/InputPara          Dataset {SCALAR}
/Info/KeyInfo            Dataset {SCALAR}
/Info/Makefile           Dataset {SCALAR}
/Info/SymConst           Dataset {SCALAR}
/Tree                    Group
/Tree/Corner             Dataset {6016, 3}
/Tree/Father             Dataset {6016}
/Tree/LBIdx              Dataset {6016}
/Tree/Sibling            Dataset {6016, 26}
/Tree/Son                Dataset {6016}
</pre>
</details>

8\. Validate the OpenMP configuration by following step 5 in
[[Quick Start: 1D Shock Tube -- CPU-only with OpenMP | Quick-Start:-1D-Shock-Tube#cpu-only-with-openmp]].
The example setup given above leads to something like
<pre>
OpenMP Diagnosis
***********************************************************************************
OMP__SCHEDULE                   DYNAMIC
OMP__SCHEDULE_CHUNK_SIZE        1
OMP__NESTED                     OFF

CPU core IDs of all OpenMP threads (tid == thread ID):
------------------------------------------------------------------------
 Rank        Host  NThread  tid-00  tid-01  tid-02  tid-03  tid-04  tid-05  tid-06  tid-07  tid-08  tid-09
    0    golub121       10       2       0       4       8      10      12      14      16      18       6
    1    golub121       10       1       3       5       7       9      11      13      15      17      19
***********************************************************************************
</pre>
> [!NOTE]
> Here the two MPI processes run on the same node `golub121` and all OpenMP threads use different CPU cores.

9\. Validate the GPU configuration by following step 3 in
[[Quick Start: 1D Shock Tube -- Hybrid OpenMP/GPU | Quick-Start:-1D-Shock-Tube#hybrid-openmpgpu]].
Especially, make sure that the MPI processes running on the
same node access different `GPU ID` (unless that is what you want).

10\. Plot a HDF5 snapshot with yt. You can use the sample script `plot_script/plot_slice.py`.
```
cd plot_script
python plot_slice_gas.py -h
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
usage: plot_slice_gas.py [-h] -s IDX_START -e IDX_END [-d DIDX] [-i PREFIX]

Plot gas density slices for the blast wave test

optional arguments:
  -h, --help    show this help message and exit
  -s IDX_START  first data index
  -e IDX_END    last data index
  -d DIDX       delta data index [1]
  -i PREFIX     data path prefix [../]
</pre>
</details>

Let's plot `Data_000010`.
``` bash
python plot_slice_gas.py -s 10 -e 10
```
<details>
<summary><u><i>Execution results</i></u></summary>

<pre>
Command-line arguments:
-------------------------------------------------------------------
plot_slice_gas.py -s 10 -e 10
-------------------------------------------------------------------

yt : [WARNING  ] 2017-12-11 20:34:09,787 Cannot determine code units ==> Use units_override to specify the units
yt : [WARNING  ] 2017-12-11 20:34:09,788 Assuming length unit = 1.0 cm
yt : [WARNING  ] 2017-12-11 20:34:09,789 Assuming time unit = 1.0 s
yt : [WARNING  ] 2017-12-11 20:34:09,789 Assuming mass unit = 1.0 g
yt : [INFO     ] 2017-12-11 20:34:09,842 Parameters: current_time              = 0.005
yt : [INFO     ] 2017-12-11 20:34:09,842 Parameters: domain_dimensions         = [32 32 32]
yt : [INFO     ] 2017-12-11 20:34:09,843 Parameters: domain_left_edge          = [ 0.  0.  0.]
yt : [INFO     ] 2017-12-11 20:34:09,844 Parameters: domain_right_edge         = [ 1.  1.  1.]
yt : [INFO     ] 2017-12-11 20:34:09,845 Parameters: cosmological_simulation   = 0
yt : [INFO     ] 2017-12-11 20:34:14,783 xlim = 0.000000 1.000000
yt : [INFO     ] 2017-12-11 20:34:14,784 ylim = 0.000000 1.000000
yt : [INFO     ] 2017-12-11 20:34:14,786 xlim = 0.000000 1.000000
yt : [INFO     ] 2017-12-11 20:34:14,787 ylim = 0.000000 1.000000
yt : [INFO     ] 2017-12-11 20:34:14,810 Making a fixed resolution buffer of (('gas', 'density')) 800 by 800
yt : [INFO     ] 2017-12-11 20:34:18,923 Saving plot Data_000010_Slice_z_density.png
</pre>
</details>

```bash
display Data_000010_Slice_z_density.png
```
<details>
<summary><u><i>Execution results</i></u></summary>

[[/images/blastwave.png | alt=blastwave]]
</details>

<br>

## Links
* [[Previous demo -- 1D Shock Tube: OpenMP with/without GPU acceleration | Quick Start: 1D Shock Tube ]]
* [[Back to the main page of Quick Start | Quick Start]]
