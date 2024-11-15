# Generate Makefile
To get the `Makefile`, please execute the following command:

```bash
python configure.py --machine=your_configuration_file [--your_arguments]
```

`your_configuration_file` is the configuration file you got from [[Machine Configuration File | Installation:-Machine-Configuration-File]], and `[--your_arguments]` should match your simulation requirements. Please check out [Option List](#option-list) for all the available options.

For example, the following command uses `configs/pleiades.config` machine configuration, sets the FFT method to `FFTW2`, and enables gravity and GPU.

``` bash
python configure.py --machine=pleiades --fftw=FFTW2 --gravity=true --gpu=true
```

> [!TIP]
> An example script `generate_make.sh` to generate Makefile can be found in each test problem folder,
e.g., `example/test_problem/Hydro/AcousticWave/generate_make.sh`.

# Option List
All available options in `configure.py` including compile-time simulation options and `configure.py` only options are listed below.
* [Configure.py Only](#configurepy-only)

Compile-time simulation options are classified into the following categories:
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

## `configure.py` Only

| Option | Value | Description |
|:---:|:---:|---|
| `-h`        | -               | Show a short help message. |
| `-lh`       | -               | Show a detailed help message. |
| `--machine` | Filename string | Select the `*.config` file under `configs` directory. |

## Physical Modules

<table>
   <tr>
      <td align=center rowspan=4 width=130px><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--model</code></td>
      <td align=center><a name="MODEL"></a> <code>MODEL</code></td>
   </tr>
   <tr><td><code>HYDRO</code>, <code>ELBDM</code></td></tr>
   <tr><td>Physical models, where <code>ELBDM</code> is for &psi;DM</td></tr>
   <tr><td>Must be set in any cases; <code>ELBDM</code> is not released yet</td></tr>
   <tr>
      <td align=center rowspan=4><code>--gravity</code></td>
      <td align=center><a name="GRAVITY"></a> <code>GRAVITY</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Enable [[gravity | Gravity]]</td></tr>
   <tr><td>Must enable <code>SUPPORT_FFTW</code>; may need to set [[FFTW2/3_PATH | Installation: External Libraries]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--particle</code></td>
      <td align=center><a name="PARTICLE"></a> <code>PARTICLE</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Enable [[particles | Particles]]</td></tr>
   <tr><td>Must enable <code>GRAVITY</code> or <code>TRACER</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--grackle</code></td>
      <td align=center><a name="SUPPORT_GRACKLE"></a> <code>SUPPORT_GRACKLE</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Enable [[GRACKLE | Chemistry and Radiation]]</td></tr>
   <tr><td>May need to set [[GRACKLE_PATH | Installation: External Libraries]]; only support <code>EOS=EOS_GAMMA/EOS_COSMIC_RAY</code>; does NOT support <code>COMOVING</code></td></tr>
</table>

### Hydro Options
-- see [[Hydro]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td align=center><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--flu_scheme</code></td>
      <td align=center><a name="FLU_SCHEME"></a> <code>FLU_SCHEME</code></td>
   </tr>
   <tr><td><code>RTVD</code>, <code>MHM</code>, <code>MHM_RP</code>, <code>CTU</code></td></tr>
   <tr><td>Hydro schemes. <code>RTVD</code>: relaxing TVD; <code>MHM</code>: MUSCL-Hancock; <code>MHM_RP</code>: VL scheme; <code>CTU</code>: corner transport upwind</td></tr>
   <tr><td><code>MHD</code> only supports <code>MHM</code>, <code>MHM_RP</code>, and <code>CTU</code>; <code>SRHD</code> only supports <code>MHM</code> and <code>MHM_RP</code>; <code>COSMIC_RAY</code> only supports <code>MHM_RP</code></td></tr>
   <tr>
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
</table>

| Option | GAMER Name | Value | Description | Restriction |
|:---:|:---:|:---:|---|---|
| `--slope` | <a name="LR_SCHEME"></a> LR_SCHEME | PLM<br>PPM | Spatial reconstruction. PLM: piecewise linear; PPM: piecewise parabolic | Useless for FLU_SCHEME=RTVD |
| `--flux` | <a name="RSOLVER"></a> RSOLVER | EXACT<br>ROE<br>HLLE<br>HLLC<br>HLLD | Riemann solvers | Useless for `FLU_SCHEME=RTVD`; `EXACT` is experimental; pure hydrodynamics supports `EXACT/ROE/HLLE/HLLC`; `MHD` supports `ROE/HLLE/HLLD`; `SRHD` and `COSMIC_RAY` support `HLLE/HLLC` |
| `--dual` | <a name="DUAL_ENERGY"></a> DUAL_ENERGY | DE_ENPY | Enable dual energy formalism | Not supported for FLU_SCHEME=RTVD |
| `--passive` | <a name="NCOMP_PASSIVE_USER"></a> NCOMP_PASSIVE_USER | &#8805; 0 | Number of user-defined passive scalars | See [[here\|Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes]] for details; not supported for FLU_SCHEME=RTVD |
| `--mhd` | <a name="MHD"></a> MHD | | Magnetohydrodynamics | |
| `--srhd` | <a name="SRHD"></a> SRHD | | Special relativistic hydrodynamics | Must adopt `EOS=EOS_TAUBMATHEWS` |
| `--cosmic_ray` | <a name="COSMIC_RAY"></a> COSMIC_RAY | | Cosmic rays | Must adopt `EOS=EOS_COSMIC_RAY` |
| `--cr_diffusion` | <a name="CR_DIFFUSION"></a> CR_DIFFUSION | | Cosmic-ray diffusion | Must enable both `COSMIC_RAY` and `MHD` |
| `--eos` | <a name="EOS"></a> EOS | EOS_GAMMA<br>EOS_ISOTHERMAL<br>EOS_COSMIC_RAY<br>EOS_TAUBMATHEWS<br>EOS_USER | [[Equation of state \|equation-of-state]] | The following options only support `EOS_GAMMA`: `FLU_SCHEME=RTVD/CTU`, `RSOLVER=EXACT/ROE`, `COMOVING`, `DUAL_ENERGY`; see also [BAROTROPIC_EOS](#BAROTROPIC_EOS) |
| `--barotropic` | <a name="BAROTROPIC_EOS"></a> BAROTROPIC_EOS || Is [EOS](#EOS) barotropic? | Must be disabled for `EOS_GAMMA` and enabled for `EOS_ISOTHERMAL` |

## Gravity Options
-- see [[Gravity]] for the related runtime parameters and other settings

| Option | GAMER Name | Value | Description | Restriction |
|:---:|:---:|:---:|---|---|
| `--pot_scheme` | <a name="POT_SCHEME"></a> POT_SCHEME | SOR<br>MG | Poisson solver. SOR: successive-overrelaxation (recommended); MG: multigrid | Must be set when GRAVITY is enabled |
| `--store_pot_ghost` | <a name="STORE_POT_GHOST"></a> STORE_POT_GHOST | | Store the ghost-zone potential (recommended when PARTICLE is enabled) | Must be enabled when both STAR_FORMATION and STORE_PAR_ACC are adopted |
| `--unsplit_gravity` | <a name="UNSPLIT_GRAVITY"></a> UNSPLIT_GRAVITY | | Use operator-unsplit method to couple gravity to the adopted physical model (recommended) | Not supported for MODEL=ELBDM |
| `--comoving` | <a name="COMOVING"></a> COMOVING | | Cosmological simulations | |

## Particle Options
-- see [[Particles]] for the related runtime parameters and other settings

| Option | GAMER Name | Value | Description | Restriction |
|:---:|:---:|:---:|---|---|
| `--tracer` | <a name="TRACER"></a> TRACER | | Enable tracer particles | |
| `--store_par_acc` | <a name="STORE_PAR_ACC"></a> STORE_PAR_ACC | | Store particle acceleration (recommended) | |
| `--star_formation` | <a name="STAR_FORMATION"></a> STAR_FORMATION | | Enable star formation | Must enable STORE_POT_GHOST when using STORE_PAR_ACC |
| `--feedback` | <a name="FEEDBACK"></a> FEEDBACK | | Enable feedback from particles to grids (and vice versa) | see [[here\|Feedback]] for details |
| `--pat_attribute` | <a name="PAR_NATT_USER"></a> PAR_NATT_USER | &#8805; 0 | Number of user-defined particle attributes | See [[here\|Adding-New-Simulations#particle-attributes]] for details |

## In Situ Python Analysis Options
-- see [[In Situ Python Analysis | In-Situ-Python-Analysis]] for the related runtime parameters and other settings

| Option | GAMER Name | Value | Description | Restriction |
|:---:|:---:|:---:|---|---|
| `--libyt` | <a name="SUPPORT_LIBYT"></a> SUPPORT_LIBYT | | Enable libyt for in situ Python analysis | May need to set [[LIBYT_PATH\|Installation: External Libraries#libyt]] |
| `--libyt_use_patch_group` | <a name="LIBYT_USE_PATCH_GROUP"></a> LIBYT_USE_PATCH_GROUP | | Use patch groups instead of patches as the grid unit for better performance (recommended) | Must enable `SUPPORT_LIBYT` |
| `--libyt_interactive` | <a name="LIBYT_INTERACTIVE"></a> LIBYT_INTERACTIVE | | Activate interactive Python prompt in in situ analysis | Must enable `SUPPORT_LIBYT` and compile libyt in interactive mode |

## Miscellaneous Options
-- AMR, GPU, parallelization, optimizations, ...

| Option | GAMER Name | Value | Description | Restriction |
|:---:|:---:|:---:|---|---|
| `--nlevel` | <a name="NLEVEL"></a> NLEVEL | &#8805; 1 | Maximum number of AMR levels including the root level. Do not confuse with the [[MAX_LEVEL \| Runtime Parameters:-Refinement#MAX_LEVEL]] runtime parameter. | |
| `--max_patch` | <a name="MAX_PATCH"></a> MAX_PATCH | &#8805; 8 | Maximum number of patches that can be allocated on each AMR level (recommended value: 1000000 or even larger since it is not memory-consuming) | |
| `--patch_size` | <a name="PATCH_SIZE"></a> PATCH_SIZE | &#8805; 8 | Number of cells along each direction in a single patch | Must be an even number |
| `--gpu` | <a name="GPU"></a> GPU | | Enable GPU acceleration | Must specify `GPU_COMPUTE_CAPABILITY` as well; may need to set [[CUDA_PATH\|Installation: External Libraries]] |
| | <a name="GPU_COMPUTE_CAPABILITY"></a> GPU_COMPUTE_CAPABILITY | Three digits | [GPU Compute Capability](https://developer.nvidia.com/cuda-gpus) (e.g., `890` for GeForce RTX 4090) | Must enable GPU |
| `--debug` | <a name="GAMER_DEBUG"></a> GAMER_DEBUG | | Run GAMER in a debug mode | |
| `--bitwise_reproducibility` | <a name="BITWISE_REPRODUCIBILITY"></a> BITWISE_REPRODUCIBILITY | | Enable [[bitwise reproducibility\|Bitwise Reproducibility]]. It may deteriorate performance, especially for runs with a large number of particles. | |
| `--timing` | <a name="TIMING"></a> TIMING | | Record the wall time of various GAMER routines in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]] (recommended) |
| `--timing_solver` | <a name="TIMING_SOLVER"></a> TIMING_SOLVER | | Record the wall time of individual GPU solvers in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]]. It will disable the CPU/GPU overlapping and thus deteriorate performance notably. | Must enable TIMING |
| `--double` | <a name="FLOAT8"></a> FLOAT8 | | Enable double precision floating-point accuracy for grid fields. Note that it could have a serious impact on GPU performance. | |
| `--double_par` | <a name="FLOAT8_PAR"></a> FLOAT8_PAR | | Enable double precision floating-point accuracy for particles. It will be set to `FLOAT8` by default. | |
| | <a name="SERIAL"></a> SERIAL | | Run GAMER in a serial mode (i.e., no MPI; but OpenMP is still supported) | Must disable LOAD_BALANCE |
| `--mpi` | <a name="LOAD_BALANCE"></a> LOAD_BALANCE | HILBERT | Enable load balancing using a space-filling curve (see [[MPI and OpenMP]]) | Must disable SERIAL; may need to set [[MPI_PATH\|Installation: External Libraries]] |
| `--openmp` | <a name="OPENMP"></a> OPENMP | | Enable OpenMP (see [[MPI and OpenMP]]) | Must set the compilation flag [[OPENMPFLAG\|Installation: Compiler and Flags]] |
| `--hdf5` | <a name="SUPPORT_HDF5"></a> SUPPORT_HDF5 | | Enable HDF5 output (see [[Outputs]]) | May need to set [[HDF5_PATH\|Installation: External Libraries]] |
| `--gsl` | <a name="SUPPORT_GSL"></a> SUPPORT_GSL | | Enable GNU scientific library | May need to set [[GSL_PATH\|Installation: External Libraries]] |
| `--fftw` | <a name="SUPPORT_FFTW"></a> SUPPORT_FFTW | FFTW2<br>FFTW3 | Enable FFTW | May need to set [[FFTW2/3_PATH\|Installation: External Libraries]] |
| `--rng` | <a name="RANDOM_NUMBER"></a> RANDOM_NUMBER | RNG_GNU_EXT<br>RNG_CPP11 | Random number generators. RNG_GNU_EXT: GNU extension `drand48_r`; RNG_CPP11: c++11 `<random>` | Use RNG_GNU_EXT for compilers supporting GNU extensions (**they may not be supported on macOS**); use RNG_CPP11 for compilers supporting c++11 (**one may need to add `-std=c++11` to `CXXFLAG`** &#8594; see [[Compiler and Flags\|Installation:-Compiler-and-Flags]]) |

> [!CAUTION]
> On macOS, we recommend using the GNU compiler and set
[[RANDOM_NUMBER | Installation:-Simulation-Options#RANDOM_NUMBER]] to `RNG_CPP11`
in the `Makefile` (or via `--rng=RNG_CPP11` in `configure.py`).

<br>

# Links
* [[Makefile configuration -- Compiler and Flags | Installation: Compiler and Flags]]
* [[Makefile configuration -- External Libraries | Installation: External Libraries]]
* [[Back to the main page of Installation | Installation]]
