## Format

All compile-time simulation options in the `Makefile` are in
the following two formats:

    SIMU_OPTION += -DOPTION1
    SIMU_OPTION += -DOPTION2=OPTION2_ADOPTED

which will enable `OPTION1` and assign `OPTION2_ADOPTED` to
`OPTION2`. For example, to (i) enable gravity and (ii) adopt the
CTU fluid scheme, set

    SIMU_OPTION += -DGRAVITY
    SIMU_OPTION += -DFLU_SCHEME=CTU

To disable an option, just comment it out with `#`. For example,
to disable gravity, use

    #SIMU_OPTION += -DGRAVITY

> [!CAUTION]
> * Option values (if any) must be set explicitly since there are no default values.
> For example, `SIMU_OPTION += -DFLU_SCHEME` without assigning any value to the option `FLU_SCHEME` is invalid.
> * Do not insert any space before and after the equal sign `=`.
> For example, use `-DFLU_SCHEME=CTU` instead of `-DFLU_SCHEME = CTU`.

## Option List

All compile-time simulation options are listed below. They are
classified into the following categories:
* [Physical Modules](#physical-modules)
* [Hydro](#hydro-options)
* [Gravity](#gravity-options)
* [Particles](#particle-options)
* [In Situ Python Analysis](#in-situ-python-analysis-options)
* [Miscellaneous](#miscellaneous-options)

> [!CAUTION]
> Some combinations are mandatory (e.g., `RSOLVER` must be set when
`FLU_SCHEME=CTU`), while some combinations are prohibited
(e.g., `PARTICLE` is not supported when both `GRAVITY` and `TRACER` are
disabled). See the "Restriction" of each option carefully.

### Physical Modules
| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="MODEL"></a> `MODEL`  | `HYDRO`<br>`ELBDM` | Physical models, where `ELBDM` is for &psi;DM | Must be set in any cases; `ELBDM` is not released yet |
| <a name="GRAVITY"></a> `GRAVITY` | | Enable [[gravity\|Gravity]] | Must enable `SUPPORT_FFTW`; may need to set [[FFTW2/3_PATH\|Installation: External Libraries]] |
| <a name="PARTICLE"></a> `PARTICLE` | | Enable [[particles\|Particles]] | Must enable `GRAVITY` or `TRACER` |
| <a name="SUPPORT_GRACKLE"></a> `SUPPORT_GRACKLE` | | Enable [[GRACKLE\|Chemistry and Radiation]] | May need to set [[GRACKLE_PATH\|Installation: External Libraries]]; only support `EOS=EOS_GAMMA/EOS_COSMIC_RAY`; doesn't support `COMOVING` |

### Hydro Options
-- see [[Hydro]] for the related runtime parameters and other settings

| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="FLU_SCHEME"></a> FLU_SCHEME | RTVD<br>MHM<br>MHM_RP<br>CTU | Hydro schemes. RTVD: relaxing TVD; MHM: MUSCL-Hancock; MHM_RP: VL scheme; CTU: corner transport upwind | `MHD` only supports `MHM`, `MHM_RP`, and `CTU`; `SRHD` only supports `MHM` and `MHM_RP`; `COSMIC_RAY` only supports `MHM_RP` |
| <a name="LR_SCHEME"></a> LR_SCHEME | PLM<br>PPM | Spatial reconstruction. PLM: piecewise linear; PPM: piecewise parabolic | Useless for FLU_SCHEME=RTVD |
| <a name="RSOLVER"></a> RSOLVER | EXACT<br>ROE<br>HLLE<br>HLLC<br>HLLD | Riemann solvers | Useless for `FLU_SCHEME=RTVD`; `EXACT` is experimental; pure hydrodynamics supports `EXACT/ROE/HLLE/HLLC`; `MHD` supports `ROE/HLLE/HLLD`; `SRHD` and `COSMIC_RAY` support `HLLE/HLLC` |
| <a name="DUAL_ENERGY"></a> DUAL_ENERGY | DE_ENPY | Enable dual energy formalism | Not supported for FLU_SCHEME=RTVD |
| <a name="NCOMP_PASSIVE_USER"></a> NCOMP_PASSIVE_USER | &#8805; 0 | Number of user-defined passive scalars | See [[here\|Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes]] for details; not supported for FLU_SCHEME=RTVD |
| <a name="MHD"></a> MHD | | Magnetohydrodynamics | |
| <a name="SRHD"></a> SRHD | | Special relativistic hydrodynamics | Must adopt `EOS=EOS_TAUBMATHEWS` |
| <a name="COSMIC_RAY"></a> COSMIC_RAY | | Cosmic rays | Must adopt `EOS=EOS_COSMIC_RAY` |
| <a name="CR_DIFFUSION"></a> CR_DIFFUSION | | Cosmic-ray diffusion | Must enable both `COSMIC_RAY` and `MHD` |
| <a name="EOS"></a> EOS | EOS_GAMMA<br>EOS_ISOTHERMAL<br>EOS_COSMIC_RAY<br>EOS_TAUBMATHEWS<br>EOS_USER | [[Equation of state \|equation-of-state]] | The following options only support `EOS_GAMMA`: `FLU_SCHEME=RTVD/CTU`, `RSOLVER=EXACT/ROE`, `COMOVING`, `DUAL_ENERGY`; see also [BAROTROPIC_EOS](#BAROTROPIC_EOS) |
| <a name="BAROTROPIC_EOS"></a> BAROTROPIC_EOS || Is [EOS](#EOS) barotropic? | Must be disabled for `EOS_GAMMA` and enabled for `EOS_ISOTHERMAL` |

### Gravity Options
-- see [[Gravity]] for the related runtime parameters and other settings

| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="POT_SCHEME"></a> POT_SCHEME | SOR<br>MG | Poisson solver. SOR: successive-overrelaxation (recommended); MG: multigrid | Must be set when GRAVITY is enabled |
| <a name="STORE_POT_GHOST"></a> STORE_POT_GHOST | | Store the ghost-zone potential (recommended when PARTICLE is enabled) | Must be enabled when both STAR_FORMATION and STORE_PAR_ACC are adopted |
| <a name="UNSPLIT_GRAVITY"></a> UNSPLIT_GRAVITY | | Use operator-unsplit method to couple gravity to the adopted physical model (recommended) | Not supported for MODEL=ELBDM |
| <a name="COMOVING"></a> COMOVING | | Cosmological simulations | |

### Particle Options
-- see [[Particles]] for the related runtime parameters and other settings

| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="TRACER"></a> TRACER | | Enable tracer particles | |
| <a name="STORE_PAR_ACC"></a> STORE_PAR_ACC | | Store particle acceleration (recommended) | |
| <a name="STAR_FORMATION"></a> STAR_FORMATION | | Enable star formation | Must enable STORE_POT_GHOST when using STORE_PAR_ACC |
| <a name="FEEDBACK"></a> FEEDBACK | | Enable feedback from particles to grids (and vice versa) | see [[here\|Feedback]] for details |
| <a name="PAR_NATT_USER"></a> PAR_NATT_USER | &#8805; 0 | Number of user-defined particle attributes | See [[here\|Adding-New-Simulations#particle-attributes]] for details |

### In Situ Python Analysis Options
-- see [[In Situ Python Analysis | In-Situ-Python-Analysis]] for the related runtime parameters and other settings

| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="SUPPORT_LIBYT"></a> SUPPORT_LIBYT | | Enable libyt for in situ Python analysis | May need to set [[LIBYT_PATH\|Installation: External Libraries#libyt]] |
| <a name="LIBYT_USE_PATCH_GROUP"></a> LIBYT_USE_PATCH_GROUP | | Use patch groups instead of patches as the grid unit for better performance (recommended) | Must enable `SUPPORT_LIBYT` |
| <a name="LIBYT_INTERACTIVE"></a> LIBYT_INTERACTIVE | | Activate interactive Python prompt in in situ analysis | Must enable `SUPPORT_LIBYT` and compile libyt in interactive mode |

### Miscellaneous Options
-- AMR, GPU, parallelization, optimizations, ...

| Option | Value | Description | Restriction |
|:---:|:---:|---|---|
| <a name="NLEVEL"></a> NLEVEL | &#8805; 1 | Maximum number of AMR levels including the root level. Do not confuse with the [[MAX_LEVEL \| Runtime Parameters:-Refinement#MAX_LEVEL]] runtime parameter. | |
| <a name="MAX_PATCH"></a> MAX_PATCH | &#8805; 8 | Maximum number of patches that can be allocated on each AMR level (recommended value: 1000000 or even larger since it is not memory-consuming) | |
| <a name="PATCH_SIZE"></a> PATCH_SIZE | &#8805; 8 | Number of cells along each direction in a single patch | Must be an even number |
| <a name="GPU"></a> GPU | | Enable GPU acceleration | Must specify `GPU_COMPUTE_CAPABILITY` as well; may need to set [[CUDA_PATH\|Installation: External Libraries]] |
| <a name="GPU_COMPUTE_CAPABILITY"></a> GPU_COMPUTE_CAPABILITY | Three digits | [GPU Compute Capability](https://developer.nvidia.com/cuda-gpus) (e.g., `890` for GeForce RTX 4090) | Must enable GPU |
| <a name="GAMER_DEBUG"></a> GAMER_DEBUG | | Run GAMER in a debug mode | |
| <a name="BITWISE_REPRODUCIBILITY"></a> BITWISE_REPRODUCIBILITY | | Enable [[bitwise reproducibility\|Bitwise Reproducibility]]. It may deteriorate performance, especially for runs with a large number of particles. | |
| <a name="TIMING"></a> TIMING | | Record the wall time of various GAMER routines in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]] (recommended) |
| <a name="TIMING_SOLVER"></a> TIMING_SOLVER | | Record the wall time of individual GPU solvers in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]]. It will disable the CPU/GPU overlapping and thus deteriorate performance notably. | Must enable TIMING |
| <a name="FLOAT8"></a> FLOAT8 | | Enable double precision floating-point accuracy for grid fields. Note that it could have a serious impact on GPU performance. | |
| <a name="FLOAT8_PAR"></a> FLOAT8_PAR | | Enable double precision floating-point accuracy for particles. It will be set to `FLOAT8` by default. | |
| <a name="SERIAL"></a> SERIAL | | Run GAMER in a serial mode (i.e., no MPI; but OpenMP is still supported) | Must disable LOAD_BALANCE |
| <a name="LOAD_BALANCE"></a> LOAD_BALANCE | HILBERT | Enable load balancing using a space-filling curve (see [[MPI and OpenMP]]) | Must disable SERIAL; may need to set [[MPI_PATH\|Installation: External Libraries]] |
| <a name="OPENMP"></a> OPENMP | | Enable OpenMP (see [[MPI and OpenMP]]) | Must set the compilation flag [[OPENMPFLAG\|Installation: Compiler and Flags]] |
| <a name="SUPPORT_HDF5"></a> SUPPORT_HDF5 | | Enable HDF5 output (see [[Outputs]]) | May need to set [[HDF5_PATH\|Installation: External Libraries]] |
| <a name="SUPPORT_GSL"></a> SUPPORT_GSL | | Enable GNU scientific library | May need to set [[GSL_PATH\|Installation: External Libraries]] |
| <a name="SUPPORT_FFTW"></a> SUPPORT_FFTW | FFTW2<br>FFTW3 | Enable FFTW | May need to set [[FFTW2/3_PATH\|Installation: External Libraries]] |
| <a name="RANDOM_NUMBER"></a> RANDOM_NUMBER | RNG_GNU_EXT<br>RNG_CPP11 | Random number generators. RNG_GNU_EXT: GNU extension `drand48_r`; RNG_CPP11: c++11 `<random>` | Use RNG_GNU_EXT for compilers supporting GNU extensions (**they may not be supported on macOS**); use RNG_CPP11 for compilers supporting c++11 (**one may need to add `-std=c++11` to `CXXFLAG`** &#8594; see [[Compiler and Flags\|Installation:-Compiler-and-Flags]]) |

<br>

## Links
* [[Makefile configuration -- Compiler and Flags | Installation: Compiler and Flags]]
* [[Makefile configuration -- External Libraries | Installation: External Libraries]]
* [[Back to the main page of Installation | Installation]]
