## Enabling GPU

To enable GPU:
1. Generate the `Makefile` and recompile the code (see [[Installation]] for details)
    1. Set
[[GPU_COMPUTE_CAPABILITY | Installation:-Machine-Configuration-File#4-GPU-compute-capability]]
according to the GPU architecture on your system in the [[configuration file | Installation:-Machine-Configuration-File]]
    2. Generate `Makefile` with
[[--gpu | Installation:-Option-List#--gpu]]=`true`
    3. Recompile the code with `make clean; make`

2. [Query the GPUs on your system](#query-gpus) [optional]

3. [Set and Validate GPU IDs](#set-and-validate-gpu-ids)


## Compilation Options

Related options:
[[--gpu | Installation:-Option-List#--gpu]], &nbsp;


## Runtime Parameters
[[Runtime parameters: GPU | Runtime-Parameters:-GPU]]

Other related parameters: none


## Remarks

### Query GPUs

To query all GPUs on a node, use the command
``` bash
nvidia-smi
```
Here is an example on a node with 2 Tesla K40m GPUs:
<pre>
+-----------------------------------------------------------------------------+
| NVIDIA-SMI 375.66                 Driver Version: 375.66                    |
|-------------------------------+----------------------+----------------------+
| GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
| Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
|===============================+======================+======================|
|   0  Tesla K40m          Off  | 0000:05:00.0     Off |                    0 |
| N/A   28C    P0    72W / 235W |   1071MiB / 11439MiB |     30%      Default |
+-------------------------------+----------------------+----------------------+
|   1  Tesla K40m          Off  | 0000:42:00.0     Off |                    0 |
| N/A   26C    P0    75W / 235W |   1071MiB / 11439MiB |     36%      Default |
+-------------------------------+----------------------+----------------------+

+-----------------------------------------------------------------------------+
| Processes:                                                       GPU Memory |
|  GPU       PID  Type  Process name                               Usage      |
|=============================================================================|
|    0     35286    C   ./gamer                                       1067MiB |
|    1     35287    C   ./gamer                                       1067MiB |
+-----------------------------------------------------------------------------+
</pre>

It shows that the
[CUDA device compute mode](http://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__TYPES.html#group__CUDART__TYPES_1g7eb25f5413a962faad0956d92bae10d0)
of both GPUs are set to `Default` (corresponding to `cudaComputeModeDefault`),
and there are currently two running jobs using GPU ID 0 and 1, respectively.


### Set and Validate GPU IDs

On a node with <var>N</var><sub>GPU</sub>, each GPU has a unique
ID in the range 0 to <var>N</var><sub>GPU</sub>-1. GAMER uses the runtime
parameter [[OPT__GPUID_SELECT | Runtime-Parameters:-GPU#OPT__GPUID_SELECT]] to set the GPU ID
associated with each MPI process.

* `OPT__GPUID_SELECT = -2`: set by CUDA runtime. Typically, this
option should work together with the `cudaComputeModeExclusive`
[CUDA device compute mode](http://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__TYPES.html#group__CUDART__TYPES_1g7eb25f5413a962faad0956d92bae10d0),
by which different MPI ranks in the same node will be assigned with
different GPUs automatically. Otherwise, all MPI ranks will use GPU 0,
which is likely undesirable. The `cudaComputeModeExclusive` compute
mode can be set by `nvidia-smi -c 1`, which requires root privileges.

* `OPT__GPUID_SELECT = -1`: set by MPI ranks. Specifically, it will set GPU ID
to **MPI_Rank % <var>N</var><sub>GPU</sub>**, where % is the integer modulus
operator. **This is the recommended method when running on a system with
multiple GPUs on each node.** However, one must be careful about the order
of MPI ranks among different nodes to ensure full utilization of all GPUs.
For example, if you have two MPI ranks with MPI_Rank=0 and 2 running a node
with <var>N</var><sub>GPU</sub>=2, both ranks will access GPU 0
(since both 0%2 and 2%2 are equal to 0) and GPU 1 will become idle,
which is undesirable. One straightforward approach is to adopt a
"SMP-style" rank ordering, by which ranks are placed consecutively until the
node is filled up, then on to the next node. More detailed illustration
can be found in the
[Blue Waters User Guide](https://bluewaters.ncsa.illinois.edu/topology-considerations).
Please also consult your system documentation.

* `OPT__GPUID_SELECT >= 0`: simply set GPU ID to `OPT__GPUID_SELECT`.
Valid inputs are 0 to <var>N</var><sub>GPU</sub>-1.

See also [[Hybrid MPI/OpenMP/GPU | MPI-and-OpenMP#hybrid-mpiopenmpgpu]].

To validate the ID and configuration of the GPU adopted by each
MPI process, search for the keyword "Device Diagnosis" in the log file
`Record__Note` generated during the initialization of GAMER. You should
see something like
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
This example shows that the MPI rank 0 is using GPU 0
on the node `golub123`, which has 2 GPUs in total.


<br>

## Links
* [[MPI and OpenMP | MPI and OpenMP]]
* [[Main page of Runtime Parameters | Runtime Parameters]]
