# Generate Makefile
To get the `Makefile`, please execute the following command:

```bash
python configure.py --machine=your_configuration_file [--your_arguments]
```

`your_configuration_file` is the configuration file you got from [[Machine Configuration File | Installation:-Machine-Configuration-File]], and `[--your_arguments]` should match your simulation requirements. Please check out [Option List](#option-list) for all the available options.

For example, the following command uses `configs/pleiades.config` machine configuration, sets the FFTW method to `FFTW2`, and enables gravity and GPU.

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
* [ELBDM](#elbdm-options)
* [Gravity](#gravity-options)
* [Particles](#particle-options)
* [Microphysics](#microphysics-options)
* [In Situ Python Analysis](#in-situ-python-analysis-options)
* [Parallelization](#parallelization-options)
* [Miscellaneous](#miscellaneous-options)

> [!CAUTION]
> Some combinations are mandatory (e.g., `RSOLVER` must be set when
`--flu_scheme=CTU`), while some combinations are prohibited
(e.g., `--particle` is not supported when both `--gravity` and `--tracer` are
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
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--model</code></td>
      <td><a name="--model"></a> <code>MODEL</code></td>
   </tr>
   <tr><td><code>HYDRO</code>, <code>ELBDM</code></td></tr>
   <tr><td>Physical models, where <code>ELBDM</code> is for &psi;DM</td></tr>
   <tr><td>Must be set in any cases; <code>ELBDM</code> is not released yet</td></tr>
   <tr>
      <td align=center rowspan=4><code>--gravity</code></td>
      <td><a name="--gravity"></a> <code>GRAVITY</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable [[gravity | Gravity]]</td></tr>
   <tr><td>Must enable <code>--fftw</code>; may need to set <code>FFTW2/3_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--particle</code></td>
      <td><a name="--particle"></a> <code>PARTICLE</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable [[particles | Particles]]</td></tr>
   <tr><td>Must enable <code>--gravity</code> or <code>--tracer</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--grackle</code></td>
      <td><a name="--grackle"></a> <code>SUPPORT_GRACKLE</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable [[GRACKLE | Chemistry and Radiation]]</td></tr>
   <tr><td>May need to set <code>GRACKLE_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]; only support <code>--eos=GAMMA/COSMIC_RAY</code>; does not support <code>--comoving</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--passive</code></td>
      <td><a name="--passive"></a> <code>NCOMP_PASSIVE_USER</code></td>
   </tr>
   <tr><td>&#8805; 0</td></tr>
   <tr><td>Number of user-defined passive scalars</td></tr>
   <tr><td>See [[here | Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes]] for details; not supported for <code>--flu_scheme=RTVD</code></td></tr>
</table>

## Hydro Options
-- see [[Hydro]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4 width=160px><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--flu_scheme</code></td>
      <td><a name="--flu_scheme"></a> <code>FLU_SCHEME</code></td>
   </tr>
   <tr><td><code>RTVD</code>, <code>MHM</code>, <code>MHM_RP</code>, <code>CTU</code></td></tr>
   <tr><td>Hydro schemes. <code>RTVD</code>: relaxing TVD; <code>MHM</code>: MUSCL-Hancock; <code>MHM_RP</code>: VL scheme; <code>CTU</code>: corner transport upwind</td></tr>
   <tr><td><code>--mhd</code> only supports <code>MHM</code>, <code>MHM_RP</code>, and <code>CTU</code>; <code>--srhd</code> only supports <code>MHM</code> and <code>MHM_RP</code>; <code>--cosmic_ray</code> only supports <code>MHM_RP</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--slope</code></td>
      <td><a name="--slope"></a> <code>LR_SCHEME</code></td>
   </tr>
   <tr><td><code>PLM</code>, <code>PPM</code></td></tr>
   <tr><td>Spatial reconstruction. <code>PLM</code>: piecewise linear; <code>PPM</code>: piecewise parabolic</td></tr>
   <tr><td>Useless for <code>--flu_scheme=RTVD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--flux</code></td>
      <td><a name="--flux"></a> <code>RSOLVER</code></td>
   </tr>
   <tr><td><code>EXACT</code>, <code>ROE</code>, <code>HLLE</code>, <code>HLLC</code>, <code>HLLD</code></td></tr>
   <tr><td>Riemann solvers</td></tr>
   <tr><td>Useless for <code>--flu_scheme=RTVD</code>; <code>EXACT</code> is experimental; pure hydrodynamics supports <code>EXACT/ROE/HLLE/HLLC</code>; <code>--mhd</code> supports <code>ROE/HLLE/HLLD</code>; <code>--srhd</code> and <code>--cosmic_ray</code> support <code>HLLE/HLLC</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--dual</code></td>
      <td><a name="--dual"></a> <code>DUAL_ENERGY</code></td>
   </tr>
   <tr><td><code>DE_ENPY</code></td></tr>
   <tr><td>Enable dual energy formalism</td></tr>
   <tr><td>Not supported for <code>--flu_scheme=RTVD</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--mhd</code></td>
      <td><a name="--mhd"></a> <code>MHD</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Magnetohydrodynamics</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--srhd</code></td>
      <td><a name="--srhd"></a> <code>SRHD</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Special relativistic hydrodynamics</td></tr>
   <tr><td>Must adopt <code>--eos=TAUBMATHEWS</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--cosmic_ray</code></td>
      <td><a name="--cosmic_ray"></a> <code>COSMIC_RAY</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Cosmic rays</td></tr>
   <tr><td>Must adopt <code>--eos=COSMIC_RAY</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--eos</code></td>
      <td><a name="--eos"></a> <code>EOS</code></td>
   </tr>
   <tr><td><code>GAMMA</code>, <code>ISOTHERMAL</code>, <code>COSMIC_RAY</code>, <code>TAUBMATHEWS</code>, <code>USER</code></td></tr>
   <tr><td>[[Equation of state | equation-of-state]]</td></tr>
   <tr><td>The following options only support <code>GAMMA</code>: <code>--flu_scheme=RTVD/CTU</code>, <code>--flux=EXACT/ROE</code>, <code>--comoving</code>, <code>--dual</code>; see also [[--barotropic | Installation:-Generate-Makefile#--barotropic]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--barotropic</code></td>
      <td><a name="--barotropic"></a> <code>BAROTROPIC_EOS</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Is [[--eos | Installation:-Generate-Makefile#--eos]] barotropic?</td></tr>
   <tr><td>Must be disabled for <code>--eos=GAMMA</code> and enabled for <code>--eos=ISOTHERMAL</code></td></tr>
</table>

## ELBDM Options
-- ELBDM is not supported yet!

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--conserve_mass</code></td>
      <td><a name="--conserve_mass"></a> <code>CONSERVE_MASS</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enforce the mass conservation</tr>
   <tr><td>Must set <code>--model=ELBDM</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--laplacian_four</code></td>
      <td><a name="--laplacian_four"></a> <code>LAPLACIAN_4TH</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable the fourth-order Laplacian</td></tr>
   <tr><td>Must set <code>--model=ELBDM</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--self_interaction</code></td>
      <td><a name="--self_interaction"></a> <code>QUARTIC_SELF_INTERACTION</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Include the quartic self-interaction potential</td></tr>
   <tr><td>Must set <code>--model=ELBDM</code> and enable <code>--gravity</code>. Does not support <code>--comoving</code></td></tr>
</table>

## Gravity Options
-- see [[Gravity]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--pot_scheme</code></td>
      <td><a name="--pot_scheme"></a> <code>POT_SCHEME</code></td>
   </tr>
   <tr><td><code>SOR</code>, <code>MG</code></td></tr>
   <tr><td>Poisson solver. <code>SOR</code>: successive-overrelaxation (recommended); <code>MG</code>: multigrid</td></tr>
   <tr><td>Must be set when <code>--gravity</code> is enabled</td></tr>
   <tr>
      <td align=center rowspan=4><code>--store_pot_ghost</code></td>
      <td><a name="--store_pot_ghost"></a> <code>STORE_POT_GHOST</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Store the ghost-zone potential (recommended when <code>--particle</code> is enabled)</td></tr>
   <tr><td>Must be enabled when both <code>--star_formation</code> and <code>--store_par_acc</code> are adopted</td></tr>
   <tr>
      <td align=center rowspan=4><code>--unsplit_gravity</code></td>
      <td><a name="--unsplit_gravity"></a> <code>UNSPLIT_GRAVITY</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Use operator-unsplit method to couple gravity to the adopted physical model (recommended)</td></tr>
   <tr><td>Not supported for <code>--model=ELBDM</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--comoving</code></td>
      <td><a name="--comoving"></a> <code>COMOVING</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Cosmological simulation</td></tr>
   <tr><td>-</td></tr>
</table>

## Particle Options
-- see [[Particles]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--tracer</code></td>
      <td><a name="--tracer"></a> <code>TRACER</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable tracer particles</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--store_par_acc</code></td>
      <td><a name="--store_par_acc"></a> <code>STORE_PAR_ACC</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Store particle acceleration (recommended)</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--star_formation</code></td>
      <td><a name="--star_formation"></a> <code>STAR_FORMATION</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable star formation</td></tr>
   <tr><td>Must enable <code>--store_pot_ghost</code> when using <code>--store_par_acc</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--feedback</code></td>
      <td><a name="--feedback"></a> <code>FEEDBACK</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable feedback from particles to grids (and vice versa)</td></tr>
   <tr><td>see [[here | Feedback]] for details</td></tr>
   <tr>
      <td align=center rowspan=4><code>--par_attribute</code></td>
      <td><a name="--par_attribute"></a> <code>PAR_NATT_USER</code></td>
   </tr>
   <tr><td>&#8805; 0</td></tr>
   <tr><td>Number of user-defined particle attributes</td></tr>
   <tr><td>See [[here | Adding-New-Simulations#particle-attributes]] for details</td></tr>
   <tr>
      <td align=center rowspan=4><code>--double_par</code></td>
      <td><a name="--double_par"></a> <code>FLOAT8_PAR</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable double-precision floating-point accuracy for particles. It will be set to <code>--double</code> by default.</td></tr>
   <tr><td>-</td></tr>
</table>

## Microphysics Options
<table>
   <tr>
      <td align=center rowspan=4><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--cr_diffusion</code></td>
      <td><a name="--cr_diffusion"></a> <code>CR_DIFFUSION</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Cosmic-ray diffusion</td></tr>
   <tr><td>Must enable both <code>--cosmic_ray</code> and <code>--mhd</code></td></tr>
</table>

## In Situ Python Analysis Options
-- see [[In Situ Python Analysis | In-Situ-Python-Analysis]] for the related runtime parameters and other settings

<table>
   <tr>
      <td align=center rowspan=4 width=180px><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--libyt</code></td>
      <td><a name="--libyt"></a> <code>SUPPORT_LIBYT</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable libyt for in situ Python analysis</td></tr>
   <tr><td>May need to set <code>LIBYT_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--libyt_patchgroup</code></td>
      <td><a name="--libyt_patchgroup"></a> <code>LIBYT_USE_PATCH_GROUP</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Use patch groups instead of patches as the grid unit for better performance (recommended)</td></tr>
   <tr><td>Must enable <code>--libyt</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--libyt_interactive</code></td>
      <td><a name="--libyt_interactive"></a> <code>LIBYT_INTERACTIVE</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Activate interactive Python prompt in in situ analysis</td></tr>
   <tr><td>Must enable <code>--libyt</code> and compile libyt with <code>INTERACTIVE_MODE</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--libyt_reload</code></td>
      <td><a name="--libyt_reload"></a> <code>LIBYT_RELOAD</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable the interactive mode of libyt. This activates Python prompt and does not shut down a simulation when there are errors in an inline Python script.</td></tr>
   <tr><td>Must enable <code>--libyt</code> and compile libyt with <code>INTERACTIVE_MODE</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--libyt_jupyter</code></td>
      <td><a name="--libyt_jupyter"></a> <code>LIBYT_JUPYTER</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Allow for in situ analysis using Jupyter Notebook / JupyterLab through libyt.</td></tr>
   <tr><td>Must enable <code>--libyt</code> and compile libyt with <code>JUPYTER_KERNEL</code>.</td></tr>
</table>

## Parallelization Options

<table>
   <tr>
      <td align=center rowspan=4 width=230px><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--openmp</code></td>
      <td><a name="--openmp"></a> <code>OPENMP</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable OpenMP (see [[MPI and OpenMP]])</td></tr>
   <tr><td>Must set the compilation flag <code>OPENMPFLAG</code> in [[configuration file | Installation:-Machine-Configuration-File#3-Compilation-flags]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--mpi</code></td>
      <td><a name="--mpi"></a> <code>LOAD_BALANCE=HILBERT</code>, <a name="SERIAL"></a> <code>SERIAL</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td><code>true</code>: Enable load balancing using a space-filling curve (see [[MPI and OpenMP]]); <code>false</code>: Run GAMER in a serial mode, but OpenMP is still supported</td></tr>
   <tr><td>May need to set <code>MPI_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--overlap_mpi</code></td>
      <td><a name="--overlap_mpi"></a> <code>OVERLAP_MPI</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Overlap MPI communication with computation.</td></tr>
   <tr><td>Does not support yet!!! Must enable `--mpi`.</td></tr>
   <tr>
      <td align=center rowspan=4><code>--gpu</code></td>
      <td><a name="--gpu"></a> <code>GPU</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable GPU acceleration</td></tr>
   <tr><td>Must specify <code>GPU_COMPUTE_CAPABILITY</code> in configuration file as well; may need to set <code>CUDA_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
</table>

## Miscellaneous Options
-- AMR, optimizations, ...

<table>
   <tr>
      <td align=center rowspan=4 width=230px><b>Option</b></td>
      <td><b>GAMER Name</b></td>
   </tr>
   <tr><td><b>Value</b></td></tr>
   <tr><td><b>Description</b></td></tr>
   <tr><td><b>Restriction</b></td></tr>
   <tr>
      <td align=center rowspan=4><code>--nlevel</code></td>
      <td><a name="--nlevel"></a> <code>NLEVEL</code></td>
   </tr>
   <tr><td>&#8805; 1</td></tr>
   <tr><td>Maximum number of AMR levels including the root level. Do not confuse with the [[MAX_LEVEL | Runtime Parameters:-Refinement#MAX_LEVEL]] runtime parameter.</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--max_patch</code></td>
      <td><a name="--max_patch"></a> <code>MAX_PATCH</code></td>
   </tr>
   <tr><td>&#8805; 8</td></tr>
   <tr><td>Maximum number of patches that can be allocated on each AMR level (recommended value: 1000000 or even larger since it is not memory-consuming)</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--patch_size</code></td>
      <td><a name="--patch_size"></a> <code>PATCH_SIZE</code></td>
   </tr>
   <tr><td>&#8805; 8</td></tr>
   <tr><td>Number of cells along each direction in a single patch</td></tr>
   <tr><td>Must be an even number</td></tr>
   <tr>
      <td align=center rowspan=4><code>--debug</code></td>
      <td><a name="--debug"></a> <code>GAMER_DEBUG</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Run GAMER in a debug mode</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--bitwise_reproducibility</code></td>
      <td><a name="--bitwise_reproducibility"></a> <code>BITWISE_REPRODUCIBILITY</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable [[bitwise reproducibility | Bitwise Reproducibility]]. It may deteriorate performance, especially for runs with a large number of particles.</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--timing</code></td>
      <td><a name="--timing"></a> <code>TIMING</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Record the wall time of various GAMER routines in the file [[Record__Timing | Simulation-Logs:-Record__Timing]] (recommended)</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--timing_solver</code></td>
      <td><a name="--timing_solver"></a> <code>TIMING_SOLVER</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Record the wall time of individual GPU solvers in the file [[Record__Timing | Simulation-Logs:-Record__Timing]]. It will disable the CPU/GPU overlapping and thus deteriorate performance notably.</td></tr>
   <tr><td>Must enable <code>TIMING</code></td></tr>
   <tr>
      <td align=center rowspan=4><code>--double</code></td>
      <td><a name="--double"></a> <code>FLOAT8</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable double-precision floating-point accuracy for grid fields. Note that it could have a serious impact on GPU performance.</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--laohu</code></td>
      <td><a name="--laohu"></a> <code>LAOHU</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Work on the NAOC Laohu GPU cluster.</td></tr>
   <tr><td>-</td></tr>
   <tr>
      <td align=center rowspan=4><code>--hdf5</code></td>
      <td><a name="--hdf5"></a> <code>SUPPORT_HDF5</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable HDF5 output (see [[Outputs]])</td></tr>
   <tr><td>May need to set <code>HDF5_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--gsl</code></td>
      <td><a name="--gsl"></a> <code>SUPPORT_GSL</code></td>
   </tr>
   <tr><td><code>true</code>, <code>false</code></td></tr>
   <tr><td>Enable GNU scientific library</td></tr>
   <tr><td>May need to set <code>GSL_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--fftw</code></td>
      <td><a name="--fftw"></a> <code>SUPPORT_FFTW</code></td>
   </tr>
   <tr><td><code>FFTW2</code>, <code>FFTW3</code></td></tr>
   <tr><td>Enable FFTW</td></tr>
   <tr><td>May need to set <code>FFTW2/3_PATH</code> in [[configuration file | Installation:-Machine-Configuration-File#1-Library-paths]]</td></tr>
   <tr>
      <td align=center rowspan=4><code>--rng</code></td>
      <td><a name="--rng"></a> <code>RANDOM_NUMBER</code></td>
   </tr>
   <tr><td><code>RNG_GNU_EXT</code>, <code>RNG_CPP11</code></td></tr>
   <tr><td>Random number generators. <code>RNG_GNU_EXT</code>: GNU extension <code>drand48_r</code>; <code>RNG_CPP11</code>: c++11 <code>&lt;random&gt;</code></td></tr>
   <tr><td>Use <code>RNG_GNU_EXT</code> for compilers supporting GNU extensions (<b>they may not be supported on macOS</b>); use <code>RNG_CPP11</code> for compilers supporting c++11 (<b>one may need to add <code>-std=c++11</code> to <code>CXXFLAG</code></b> &#8594; see [[compilation flags | Installation:-Machine-Configuration-File#3-Compilation-flags]])</td></tr>
</table>

> [!CAUTION]
> On macOS, we recommend using the GNU compiler and set
[[--rng | Installation:-Generate-Makefile#--rng]] to `RNG_CPP11`.

<br>

# Links
* [[Configuration file | Installation:-Machine-Configuration-File]]
* [[External Libraries | Installation: External Libraries]]
* [[Back to the main page of Installation | Installation]]
