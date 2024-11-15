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
      <td align=center rowspan=4><code>--slope</code></td>
      <td align=center><a name="LR_SCHEME"></a> <code>LR_SCHEME</code></td>
   </tr>
   <tr><td><code>PLM</code>, <code>PPM</code></td></tr>
   <tr><td>Spatial reconstruction. <code>PLM</code>: piecewise linear; <code>PPM</code>: piecewise parabolic</td></tr>
   <tr><td>Useless for <code>FLU_SCHEME=RTVD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--flux</code></td>
      <td align=center><a name="RSOLVER"></a> <code>RSOLVER</code></td>
   </tr>
   <tr><td><code>EXACT</code>, <code>ROE</code>, <code>HLLE</code>, <code>HLLC</code>, <code>HLLD</code></td></tr>
   <tr><td>Riemann solvers</td></tr>
   <tr><td>Useless for <code>FLU_SCHEME=RTVD</code>; <code>EXACT</code> is experimental; pure hydrodynamics supports <code>EXACT/ROE/HLLE/HLLC</code>; <code>MHD</code> supports <code>ROE/HLLE/HLLD</code>; <code>SRHD</code> and <code>COSMIC_RAY</code> support <code>HLLE/HLLC</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--dual</code></td>
      <td align=center><a name="DUAL_ENERGY"></a> <code>DUAL_ENERGY</code></td>
   </tr>
   <tr><td><code>DE_ENPY</code></td></tr>
   <tr><td>Enable dual energy formalism</td></tr>
   <tr><td>Not supported for <code>FLU_SCHEME=RTVD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--passive</code></td>
      <td align=center><a name="NCOMP_PASSIVE_USER"></a> <code>NCOMP_PASSIVE_USER</code></td>
   </tr>
   <tr><td>&#8805; 0</td></tr>
   <tr><td>Number of user-defined passive scalars</td></tr>
   <tr><td>See [[here | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes]] for details; not supported for <code>FLU_SCHEME=RTVD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--mhd</code></td>
      <td align=center><a name="MHD"></a> <code>MHD</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Magnetohydrodynamics</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--srhd</code></td>
      <td align=center><a name="SRHD"></a> <code>SRHD</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Special relativistic hydrodynamics</td></tr>
   <tr><td>Must adopt <code>EOS=EOS_TAUBMATHEWS</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--cosmic_ray</code></td>
      <td align=center><a name="COSMIC_RAY"></a> <code>COSMIC_RAY</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Cosmic rays</td></tr>
   <tr><td>Must adopt <code>EOS=EOS_COSMIC_RAY</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--cr_diffusion</code></td>
      <td align=center><a name="CR_DIFFUSION"></a> <code>CR_DIFFUSION</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Cosmic-ray diffusion</td></tr>
   <tr><td>Must enable both <code>COSMIC_RAY</code> and <code>MHD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--eos</code></td>
      <td align=center><a name="EOS"></a> <code>EOS</code></td>
   </tr>
   <tr><td><code>EOS_GAMMA</code>, <code>EOS_ISOTHERMAL</code>, <code>EOS_COSMIC_RAY</code>, <code>EOS_TAUBMATHEWS</code>, <code>EOS_USER</code></td></tr>
   <tr><td>[[Equation of state |equation-of-state]]</td></tr>
   <tr><td>The following options only support <code>EOS_GAMMA</code>: <code>FLU_SCHEME=RTVD/CTU</code>, <code>RSOLVER=EXACT/ROE</code>, <code>COMOVING</code>, <code>DUAL_ENERGY</code>; see also [BAROTROPIC_EOS](#BAROTROPIC_EOS)</td></tr>
   <tr>
      <td align=center rowspan=4><code>--barotropic</code></td>
      <td align=center><a name="BAROTROPIC_EOS"></a> <code>BAROTROPIC_EOS</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Is [EOS](#EOS) barotropic?</td></tr>
   <tr><td>Must be disabled for <code>EOS_GAMMA</code> and enabled for <code>EOS_ISOTHERMAL</code></td></tr>
</table>

## Gravity Options
-- see [[Gravity]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td align=center><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--pot_scheme</code></td>
      <td align=center><a name="POT_SCHEME"></a> <code>POT_SCHEME</code></td>
   </tr>
   <tr><td><code>SOR</code>, <code>MG</code></td></tr>
   <tr><td>Poisson solver. <code>SOR</code>: successive-overrelaxation (recommended); <code>MG</code>: multigrid</td></tr>
   <tr><td>Must be set when <code>GRAVITY</code> is enabled</td></tr>
   <tr>
      <td align=center rowspan=4><code>--store_pot_ghost</code></td>
      <td align=center><a name="STORE_POT_GHOST"></a> <code>STORE_POT_GHOST</code></td>
   </tr>
   <tr><td>-</td></tr>
   <tr><td>Store the ghost-zone potential (recommended when <code>PARTICLE</code> is enabled)</td></tr>
   <tr><td>Must be enabled when both <code>STAR_FORMATION</code> and <code>STORE_PAR_ACC</code> are adopted</td></tr>
   <tr>
| `--unsplit_gravity` | <a name="UNSPLIT_GRAVITY"></a> UNSPLIT_GRAVITY | | Use operator-unsplit method to couple gravity to the adopted physical model (recommended) | Not supported for MODEL=ELBDM |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--comoving` | <a name="COMOVING"></a> COMOVING | | Cosmological simulations | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
</table>

## Particle Options
-- see [[Particles]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td align=center><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
| `--tracer` | <a name="TRACER"></a> TRACER | | Enable tracer particles | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--store_par_acc` | <a name="STORE_PAR_ACC"></a> STORE_PAR_ACC | | Store particle acceleration (recommended) | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--star_formation` | <a name="STAR_FORMATION"></a> STAR_FORMATION | | Enable star formation | Must enable STORE_POT_GHOST when using STORE_PAR_ACC |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--feedback` | <a name="FEEDBACK"></a> FEEDBACK | | Enable feedback from particles to grids (and vice versa) | see [[here\|Feedback]] for details |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--pat_attribute` | <a name="PAR_NATT_USER"></a> PAR_NATT_USER | &#8805; 0 | Number of user-defined particle attributes | See [[here\|Adding-New-Simulations#particle-attributes]] for details |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
</table>

## In Situ Python Analysis Options
-- see [[In Situ Python Analysis | In-Situ-Python-Analysis]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td align=center><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
| `--libyt` | <a name="SUPPORT_LIBYT"></a> SUPPORT_LIBYT | | Enable libyt for in situ Python analysis | May need to set [[LIBYT_PATH\|Installation: External Libraries#libyt]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--libyt_use_patch_group` | <a name="LIBYT_USE_PATCH_GROUP"></a> LIBYT_USE_PATCH_GROUP | | Use patch groups instead of patches as the grid unit for better performance (recommended) | Must enable `SUPPORT_LIBYT` |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--libyt_interactive` | <a name="LIBYT_INTERACTIVE"></a> LIBYT_INTERACTIVE | | Activate interactive Python prompt in in situ analysis | Must enable `SUPPORT_LIBYT` and compile libyt in interactive mode |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
</table>

## Miscellaneous Options
-- AMR, GPU, parallelization, optimizations, ...

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td align=center><b>GAMER Name</b></td>
   </tr>
   <tr><td align=center><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
| `--nlevel` | <a name="NLEVEL"></a> NLEVEL | &#8805; 1 | Maximum number of AMR levels including the root level. Do not confuse with the [[MAX_LEVEL \| Runtime Parameters:-Refinement#MAX_LEVEL]] runtime parameter. | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--max_patch` | <a name="MAX_PATCH"></a> MAX_PATCH | &#8805; 8 | Maximum number of patches that can be allocated on each AMR level (recommended value: 1000000 or even larger since it is not memory-consuming) | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--patch_size` | <a name="PATCH_SIZE"></a> PATCH_SIZE | &#8805; 8 | Number of cells along each direction in a single patch | Must be an even number |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--gpu` | <a name="GPU"></a> GPU | | Enable GPU acceleration | Must specify `GPU_COMPUTE_CAPABILITY` as well; may need to set [[CUDA_PATH\|Installation: External Libraries]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| | <a name="GPU_COMPUTE_CAPABILITY"></a> GPU_COMPUTE_CAPABILITY | Three digits | [GPU Compute Capability](https://developer.nvidia.com/cuda-gpus) (e.g., `890` for GeForce RTX 4090) | Must enable GPU |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--debug` | <a name="GAMER_DEBUG"></a> GAMER_DEBUG | | Run GAMER in a debug mode | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--bitwise_reproducibility` | <a name="BITWISE_REPRODUCIBILITY"></a> BITWISE_REPRODUCIBILITY | | Enable [[bitwise reproducibility\|Bitwise Reproducibility]]. It may deteriorate performance, especially for runs with a large number of particles. | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--timing` | <a name="TIMING"></a> TIMING | | Record the wall time of various GAMER routines in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]] (recommended) |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--timing_solver` | <a name="TIMING_SOLVER"></a> TIMING_SOLVER | | Record the wall time of individual GPU solvers in the file [[Record__Timing \| Simulation-Logs:-Record__Timing]]. It will disable the CPU/GPU overlapping and thus deteriorate performance notably. | Must enable TIMING |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--double` | <a name="FLOAT8"></a> FLOAT8 | | Enable double precision floating-point accuracy for grid fields. Note that it could have a serious impact on GPU performance. | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--double_par` | <a name="FLOAT8_PAR"></a> FLOAT8_PAR | | Enable double precision floating-point accuracy for particles. It will be set to `FLOAT8` by default. | |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| | <a name="SERIAL"></a> SERIAL | | Run GAMER in a serial mode (i.e., no MPI; but OpenMP is still supported) | Must disable LOAD_BALANCE |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--mpi` | <a name="LOAD_BALANCE"></a> LOAD_BALANCE | HILBERT | Enable load balancing using a space-filling curve (see [[MPI and OpenMP]]) | Must disable SERIAL; may need to set [[MPI_PATH\|Installation: External Libraries]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--openmp` | <a name="OPENMP"></a> OPENMP | | Enable OpenMP (see [[MPI and OpenMP]]) | Must set the compilation flag [[OPENMPFLAG\|Installation: Compiler and Flags]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--hdf5` | <a name="SUPPORT_HDF5"></a> SUPPORT_HDF5 | | Enable HDF5 output (see [[Outputs]]) | May need to set [[HDF5_PATH\|Installation: External Libraries]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--gsl` | <a name="SUPPORT_GSL"></a> SUPPORT_GSL | | Enable GNU scientific library | May need to set [[GSL_PATH\|Installation: External Libraries]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--fftw` | <a name="SUPPORT_FFTW"></a> SUPPORT_FFTW | FFTW2<br>FFTW3 | Enable FFTW | May need to set [[FFTW2/3_PATH\|Installation: External Libraries]] |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr>
| `--rng` | <a name="RANDOM_NUMBER"></a> RANDOM_NUMBER | RNG_GNU_EXT<br>RNG_CPP11 | Random number generators. RNG_GNU_EXT: GNU extension `drand48_r`; RNG_CPP11: c++11 `<random>` | Use RNG_GNU_EXT for compilers supporting GNU extensions (**they may not be supported on macOS**); use RNG_CPP11 for compilers supporting c++11 (**one may need to add `-std=c++11` to `CXXFLAG`** &#8594; see [[Compiler and Flags\|Installation:-Compiler-and-Flags]]) |
      <td align=center rowspan=4><code></code></td>
      <td align=center><a name=""></a> <code></code></td>
   </tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
   <tr><td></td></tr>
</table>

> [!CAUTION]
> On macOS, we recommend using the GNU compiler and set
[[RANDOM_NUMBER | Installation:-Simulation-Options#RANDOM_NUMBER]] to `RNG_CPP11`
in the `Makefile` (or via `--rng=RNG_CPP11` in `configure.py`).

<br>

# Links
* [[Makefile configuration -- Compiler and Flags | Installation: Compiler and Flags]]
* [[Makefile configuration -- External Libraries | Installation: External Libraries]]
* [[Back to the main page of Installation | Installation]]
