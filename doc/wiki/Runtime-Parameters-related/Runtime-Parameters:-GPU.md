Parameters described on this page:
[OPT__GPUID_SELECT](#OPT__GPUID_SELECT), &nbsp;
[FLU_GPU_NPGROUP](#FLU_GPU_NPGROUP), &nbsp;
[POT_GPU_NPGROUP](#POT_GPU_NPGROUP), &nbsp;
[CHE_GPU_NPGROUP](#CHE_GPU_NPGROUP), &nbsp;
[GPU_NSTREAM](#GPU_NSTREAM) &nbsp;


Parameters below are shown in the format: &ensp; **`Name` &ensp; (Valid Values) &ensp; [Default Value]**

<a name="OPT__GPUID_SELECT"></a>
* #### `OPT__GPUID_SELECT` &ensp; (-2=CUDA, -1=MPI rank, &#8805;0=input) &ensp; [-1]
    * **Description:**
See [[Set and Validate GPU IDs | GPU#set-and-validate-gpu-ids]].
    * **Restriction:**
Must be smaller than the total number of GPUs in a node.

<a name="FLU_GPU_NPGROUP"></a>
* #### `FLU_GPU_NPGROUP` &ensp; (>0; &#8804;0 &#8594; set to default) &ensp; [depend on the GPU spec]
    * **Description:**
Number of patch groups updated by the GPU/CPU fluid solvers at a single time.
See also [[Performance Optimizations: GPU | Performance Optimizations:-GPU]].
    * **Restriction:**
Must be a multiple of [GPU_NSTREAM](#GPU_NSTREAM).

<a name="POT_GPU_NPGROUP"></a>
* #### `POT_GPU_NPGROUP` &ensp; (>0; &#8804;0 &#8594; set to default) &ensp; [depend on the GPU spec]
    * **Description:**
Number of patch groups updated by the GPU/CPU Poisson solvers at a single time.
See also [[Performance Optimizations: GPU | Performance Optimizations:-GPU]].
    * **Restriction:**
Must be a multiple of [GPU_NSTREAM](#GPU_NSTREAM).

<a name="CHE_GPU_NPGROUP"></a>
* #### `CHE_GPU_NPGROUP` &ensp; (>0; &#8804;0 &#8594; set to default) &ensp; [depend on the GPU spec]
    * **Description:**
Number of patch groups updated by the GPU/CPU GRACKLE solvers at a single time.
See also [[Performance Optimizations: GPU | Performance Optimizations:-GPU]].
The GPU version is currently not supported.
    * **Restriction:**

<a name="GPU_NSTREAM"></a>
* #### `GPU_NSTREAM` &ensp; (>0; &#8804;0 &#8594; set to default) &ensp; [depend on the GPU spec]
    * **Description:**
Number of CUDA streams for the asynchronous memory copy between CPU and GPU.
See also [[Performance Optimizations: GPU | Performance Optimizations:-GPU]].
    * **Restriction:**
See the restrictions on [FLU_GPU_NPGROUP](#FLU_GPU_NPGROUP) and
[POT_GPU_NPGROUP](#POT_GPU_NPGROUP).


## Remarks


<br>

## Links
* [[Main page of Runtime Parameters | Runtime Parameters]]
* [[Main page of GPU | GPU]]
