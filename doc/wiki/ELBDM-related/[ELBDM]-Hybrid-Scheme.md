The hybrid scheme in GAMER is designed to efficiently simulate the behavior of Extremely Light Bosonic Dark Matter (ELBDM; also referred to as Fuzzy Dark Matter or Wave Dark Matter) in cosmological simulations.

# Quick Start
1. Switch to `psidm` branch
``` bash
git checkout psidm
```
2. Copy test problem
``` bash
cp -r example/test_problem/LSS_Hybrid bin/
```
3. Compile GAMER
``` bash
cd src
cp ../bin/LSS_Hybrid/generate_make.sh ./
sh generate_make.sh
```
4. Download initial condition
``` bash
cd ../bin/LSS_Hybrid/
sh download_light_halo.sh
```
5. **To use spectral interpolation**
   * Download table

   ``` bash
   sh download_spectral_interpolation_tables.sh
   ```
   * Set `OPT__FLU_INT_SCHEME` and `OPT__REF_FLU_INT_SCHEME` to `8` in `Input__Parameter`
6. Run GAMER
``` bash
mpirun -map-by ppr:2:socket:pe=8 --report-bindings ./gamer 1>>log 2>&1
```
7. Plot results at $z=0$:
   * Slice plots
   ``` bash
   python plot_script/plot_slice.py -s 69 -e 69
   ```
   <details>
   <summary><u><i>Execution results</i></u></summary>

   ![Data_000069_Lv_10_Slice_z_density_x1](https://github.com/gamer-project/gamer/assets/6187378/d51ed4b8-af8f-4df2-a89a-f4aef66793e9)
   ![Data_000069_Lv_10_Slice_z_density_x3](https://github.com/gamer-project/gamer/assets/6187378/14ad6aff-c4f1-4f55-88c5-d623e6265dac)
   ![Data_000069_Lv_10_Slice_z_density_x10](https://github.com/gamer-project/gamer/assets/6187378/a891a4b8-3814-4131-a42c-14ba297392af)
   </details>

   * Projection plots
   ``` bash
   python plot_script/plot_projection.py -s 69 -e 69
   ```
   <details>
   <summary><u><i>Execution results</i></u></summary>

   ![Data_000069_Proj_x_density](https://github.com/gamer-project/gamer/assets/6187378/b412cf7e-8de5-4817-a039-58fecec72d13)
   ![Data_000069_Proj_x_density_grid](https://github.com/gamer-project/gamer/assets/6187378/4f416be0-0be2-4468-9182-9a77879545be)
   </details>

   * Power spectrum
   ``` bash
   python plot_script/plot_power_spectrum.py -s 0 -e 69
   ```
   <details>
   <summary><u><i>Execution results</i></u></summary>

   ![Data_000069_PS](https://github.com/gamer-project/gamer/assets/6187378/260649c1-69e7-4a44-8dc9-06ccf6fe798d)
   </details>

The wave scheme (`MATMUL`) is used for the dark grey and black grids.
The relative mass conservation error at $z=0$ with spectral interpolation on is $6.72e-04$.


# Motivation

The wave-based solver in GAMER uses a finite-difference or spectral method to discretise the Schrödinger equation.  The drawback of solving the wave formulation is the need to resolve the de Broglie wavelength of the matter wave even in regions where the density is low and smooth. The reason is that the velocity is related to the gradient of the phase of the wave function, i.e. a given velocity translates into a phase that varies on the scale of the de Broglie wavelength $\lambda_{dB} = \frac{2\pi}{m v}$.  If the latter is not resolved, the velocity field is not represented correctly. This is in contrast to conventional CDM simulations with AMR that only require higher spatial resolution in regions with higher density. In an FDM simulation with $m = 10^{-22}$ eV and a velocity of $v = 100$ km/s, the de Broglie wavelength $\lambda_{dB} \sim 1.2$ kpc is much smaller than the box size required for a large cosmological simulation. Because of their high computational demands, existing wave simulations are therefore usually limited to small box sizes.

In contrast to wave simulations, simulations of the fluid formulation do not need to resolve the de Broglie wavelength to correctly capture large-scale dynamics. In fact, they even allow the adoption of a Lagrangian picture via smoothed-particle hydrodynamics methods incorporating the quantum pressure. The possible simulation volumes are much closer to those attainable in traditional $N$-body and smoothed-particle hydrodynamics approaches for CDM. Their biggest drawback, however, is that they fail to correctly resolve regions of interference because the quantum pressure is ill-defined in regions of vanishing density. Therefore, fluid simulations of the Schrödinger-Poisson equations are generally not trustworthy on small scales.

The ELBDM hybrid scheme solves the fluid formulation of the Schrödinger-Poisson equations on a coarse grid on large scales and the wave formulation on a refined grid on small scales thereby combining the advantages of both approaches: The coarse grid on large scales need not resolve the de Broglie wavelength while the wave formulation on small scales correctly describes interference effects.

## Use cases
The hybrid scheme is particularly suited for cosmological simulations where a significant part of the simulation domain features a smooth density field and high infall velocities.

- Simulating the formation and evolution of cosmic structures under the influence of ELBDM physics.
- Investigating the core-halo structure of dark matter-dominated galaxies.
- Exploring the role of quantum pressure and interference in dark matter dynamics.

# Setup of the Hybrid Scheme for ELBDM in GAMER

## Configuring the configure.py

The configure.py script is the starting point for setting up a simulation. For the hybrid scheme, ensure the following flags are set:

- `--elbdm_scheme = HYBRID`: Determines the scheme type for ELBDM simulations. There are two options:
    - `WAVE`: Wave-only scheme.
    - `HYBRID`: Fluid-wave hybrid scheme.
- `--hybrid_scheme`: This option is specific to the hybrid scheme, defining the fluid dynamics treatment:
    - `UPWIND`: First-order upwind scheme. It is diffusive but stable.
    - `FROMM`: Second-order Fromm scheme. It has no limiter but can be unstable.
    - `MUSCL`: Second-order MUSCL scheme. It includes a limiter, balancing accuracy and stability.

    For zoom-in and fluid-only simulations, use `--hybrid_scheme MUSCL` since unrefined regions may become unstable otherwise. For fully refined simulations use `--hybrid_scheme FROMM` since it is slightly faster and more accurate than `MUSCL`. In case you encounter instabilities, test the more diffusive, first-order scheme `--hybrid_scheme UPWIND`.

## Configuring the `Input__Parameter` file

- `ELBDM_FIRST_WAVE_LEVEL` determines on which level the hybrid scheme switches to the wave solver. This level requires careful consideration as it directly impacts the accuracy and performance of your simulation. It must be set based on trial and error. If it is set too low, the hybrid scheme becomes inefficient because a large part of the simulation domain will use the wave solver. If it is set too high, simulation regions with destructive interference tend to be overrefined. See the `LSS_Hybrid` test problem for a rough estimate of how to set this parameter.
- `OPT__FLAG_INTERFERENCE` must be set in conjunction with the `Input__Flag_Interference` file. It refines fluid levels based on the quantum pressure and should generally be used in all hybrid simulations to ensure that regions where the fluid scheme fails are refined.
- The `DT__HYBRID_*` parameters control the time-step size based on different aspects of the hybrid solver. Set these to negative for automatic determination, which is recommended unless you have specific requirements. The constants were determined empirically. `DT__HYBRID_VELOCITY` with a value of up to `3.5` has never shown any instability in tests, but is set to `1.0` by default.
- `OPT__LB_EXCHANGE_FATHER` is crucial for proper load balancing in MPI runs with `ELBDM_MATCH_PHASE` and should always be set to 1 for the hybrid scheme. Its purpose is to ensure that all MPI ranks exchange information about the father patches of physical patches. This is required for backward phase matching in the restriction operation.
- `ELBDM_MATCH_PHASE` ensures that the phase at the fluid-wave-level transitions is unwrapped during the restriction operation. Alternatively, one could unwrap the phase directly in the fluid solver. The latter approach has the advantage that one could turn off `ELBDM_MATCH_PHASE`  (and `OPT__LB_EXCHANGE_FATHER` in MPI runs), but has the disadvantage that one needs to ensure that the resolution on the fluid solver levels is high enough to resolve $2\pi$ phase jumps.

## Configuring the `Input__Flag_Interference` file

The `Input__Flag_Interference` file is used in GAMER to control when to refine regions with interference within the ELBDM hybrid scheme. This file allows users to specify parameters that determine when and where to apply interference conditions throughout different levels of the simulation.

The `Input__Flag_Interference` file consists of rows of parameters, each corresponding to a specific level in the simulation. The columns represent different conditions that must be met for interference to be flagged. The parameters are as follows:

    Level: The refinement level in the adaptive mesh at which the parameters will apply.
    QP: The quantum potential threshold for flagging interference.
    Density: The matter density condition. In the given file, this is set to 0, implying that density is not used as a criterion.
    PhaseLap: The minimum phase laplacian to consider for interference effects.
    OnlyAtExtrema: A flag to indicate whether to consider interference only at extrema of the wave function (0 = No, 1 = Yes).

An example configuration might look as follows

    Level    QP    Density  PhaseLap  OnlyAtExtrema

    0      0.03     0       1.0          0
    1      0.03     0       1.0          0
    2      0.03     0       1.0          0
    ...
    10     0.03     0       1.0          0

In the above configuration, for all levels from 0 to 10, the simulation will flag interference when the quantum potential exceeds 0.03 or the phase Laplacian is at least 1.0, regardless of density and not exclusively at wave function extrema. The thresholds were determined empirically. Raising the Density parameter prevents refinement in low-density regions, but can make the scheme unstable. Setting `OnlyAtExtrema = 1` has also been observed to be unstable in cosmological simulations.

# Test problems
The ELBDM hybrid scheme is shipped with a series of dedicated test problems highlighting different features of the hybrid scheme.

## VortexPairLinear_Hybrid

This test case simulates a pair of vortices moving linearly in a two-dimensional space. It evolves the wave function of the vortices over a fixed period under periodic boundary conditions. The initial setup involves a background wave function with an added periodic component along the x and y axes. If needed, a wave in the z-direction can be added as an option. The fluid scheme is used on level 0 and the wave scheme is used on level 1 close to the vortices. This test demonstrates refinement using `Input__Flag_Interference` and also demonstrates the $2\pi$ phase jump connecting the two vortices. By enabling a wave in the z-direction, this test can also be used to verify the properties of the hybrid scheme in a true 3D simulation.

Characteristics:

    2D simulation with a possible addition in the z-direction.
    Vortices move linearly with a specific wave function described by given parameters (A, kx, ky, Omega, Phase0).
    If mass conservation is important, it can be ensured by enabling CONSERVE_MASS due to the periodic boundaries.

Run the test case in the default configuration (purely 2D) and plot the density and phase field at the final time using `python plot_slice -s 20 -e 20 -i ./`. The L1 errors are as follows: error of density = $3.89e-02$ and error of phase $1.71e-03$. The result should look as follows:

![Data_000020_z_axis](https://github.com/gamer-project/gamer/assets/6187378/b52a00b3-944d-498b-b288-510aa1996cc6)

## VortexPairRotating_Hybrid

This test problem is aimed at simulating a rotating pair of vortices in two dimensions. The boundary conditions for this simulation are based on the analytical solution, and it involves a wave function that is a function of the radial distance and azimuthal angle. The fluid scheme is used on level 0 and the wave scheme is used on level 1 close to the vortices.  Like the linear vortex pair test, this test demonstrates refinement using `Input__Flag_Interference` and also demonstrates the $2\pi$ phase jump connecting the two vortices.

Characteristics:

    2D simulation of a rotating vortex pair.
    The wave function is described in polar coordinates and evolves over time.
    This simulation references the work of Tzihong Chiueh and others, which discusses vortex turbulence in the context of linear Schrödinger wave mechanics.

Run the test case in the default configuration (purely 2D) and plot the density and phase field at the final time using `python plot_slice -s 9 -e 9 -i ./`. The L1 errors are as follows: error of density = $9.57e-02$ and error of phase $8.42e-02$. The result should look as follows:

![Data_000009_z_axis](https://github.com/gamer-project/gamer/assets/6187378/2b54547b-215c-4b0c-ac1b-6884d801f7ca)



## Perturbation

This simulation investigates the gravitational collapse that can occur in one, two, or three dimensions. The focus is on the evolution of small-amplitude plane waves superimposed on a homogeneous background density, a setup that can help understand the early stages of structure formation in the universe.

Default setup and characteristics:

    Designed to study gravitational collapse in multiple dimensions.
    Periodic boundary conditions are used.
    The hybrid scheme's refinement capabilities are showcased using Input__Flag_Interference.
    Compared to a comoving reference simulation, it supports multi-dimensional collapse without assuming a comoving frame.

Simulation parameters include:

    Perturbation_Amplitude: Maximum amplitude of the perturbations.
    Perturbation_BgAmplitude: Amplitude of the background density.
    Perturbation_N: Number of Fourier modes in the perturbation, with a maximum of 4.
    Perturbation_NDim: Dimensionality of the problem, with unused dimensions set to zero amplitude.

Run the test case in the default configuration (purely 2D) and plot the density and phase field at the final time using `python plot_slice -s 20 -e 20 -i ./`.
The result should look as follows:
![Data_000020_z_axis](https://github.com/gamer-project/gamer/assets/6187378/53d87857-fdce-46ec-b86d-32f36cff5945)

Note that the initial resolution is lower and adaptive mesh refinement is used in this test problem.

A reference solution obtained with a $1024^2$ cells simulation using the base-level spectral solver produces the following plot:
![Data_000020_z_axis](https://github.com/gamer-project/gamer/assets/6187378/b3bac9c5-b9e9-4952-95fe-c49c9af1fe97)

It is obvious that the two disagree. This is confirmed by plotting the respective 2D power spectra with the script `example/yt/extract_power_spectrum.py`. Therefore, it is clear that the default refinement settings still need modifications to obtain reliable results.
![i=20](https://github.com/gamer-project/gamer/assets/6187378/a3a59d29-1aa1-401e-a7c4-2954ec573c4f)


## LSS_Hybrid
This simulation is aimed at modeling cosmological large-scale structures using the hybrid scheme. The test problem is comprehensive, including scripts for setting up initial conditions (IC), converting wave IC to hybrid IC, and providing information on resolution and computational details relevant to the cosmological context.

Default setup and characteristics:

    Initial conditions can be obtained through the provided scripts, which also facilitate the conversion to the hybrid IC format.
    The simulation uses a base-level resolution and transitions between a fluid solver at lower levels to a wave solver at higher refinement levels.
    The resolution is set to strike a balance between performance and accuracy.
    Parameters such as ELBDM_MASS, A_INIT, OMEGA_M0, HUBBLE0, and BOX_SIZE are specified to define the cosmological framework of the simulation.

For more detailed information see the `README` file of the test problem and for reference plots, see the quick start section.

# The algorithm
The crucial issue in developing a hybrid code combining a fluid solver with a wave solver lies in reconstructing the wave function from the Madelung formulation of the Schrödinger-Poisson equations.
The Madelung transform
$$\psi(x, t) \equiv \sqrt{\frac{\rho(x, t)}{m}} e^{i S(x, t)}$$
splits the wave function up into its norm and its phase. We then take the gradient of the phase to derive the velocity:
$$v = \frac{\hbar}{m a}\nabla S.$$

It is clear that the wave function can be locally reconstructed from the density and phase fields $\rho(x, t)$ and $S(x, t)$. Reconstructing it from the velocity would either require solving a Poisson equation $\Delta v = S$ or integrating the velocity field via a suitably chosen line integral. The ELBDM Hybrid Scheme instead treats the density and phase fields as fundamental and evolves them using the Hamilton-Jacobi-Madelung equations
$$\partial_t  \rho + 3 H \rho + \frac{1}{a} \nabla \cdot \left(\rho \frac{\hbar}{ma} \nabla S\right)  = 0$$
$$\partial_t S + H S + \frac{\hbar}{2m a^2}  \left(\nabla S\right)^2 + \frac{m}{\hbar} \phi + \frac{m}{\hbar a^2} Q_P = 0 $$
$$\Delta \phi(x, t) - 4\pi G a^2 \rho_b(t) \delta(x, t) = 0$$
This allows the unique reconstruction of the wave function $\psi(x, t)$ without additional computational overhead. At the same time, it leads to what we call the _reverse boundary matching problem_ in the following.
The velocity $v$ can be easily obtained from the wave function via differentiation. Yet, the phase is only determined by the wave function up to a multiple of $2\pi$:
$$S(x, t) = \arctan\left(\frac{\Im(\psi(x, t))}{\Re(\psi(x, t)}\right).$$
We must therefore determine the correct phase by requiring continuity at the matching boundary. On a discrete grid, this translates into the requirement that the phase field $S(x, t)$ must not change by more than $2 \pi$ between neighbouring grid points at the matching boundary. In other words, we must resolve the de Broglie wavelength at the matching boundary in order for the reverse boundary matching problem to admit a unique solution. This immediately shows that a useful hybrid scheme based on this approach necessarily requires an AMR algorithm with at least two refinement levels for the phase equation: An outer refinement level where the grids need not resolve the de Broglie wavelength and a second refinement level where the de Broglie wavelength is resolved and the reverse boundary matching problem has a unique solution.

## Implementation
We evolve the continuity equation using a first-order upwind scheme (`UPWIND`), the second-order Fromm scheme without limiter (`FROMM`) or a MUSCL scheme with a linear subgrid model together with the van Albada limiter (`MUSCL`). As for the HJ scheme, we use a finite difference scheme. The convection term is discretised via the Sethian-Osher flux. The velocities at the cell faces are then computed as regular finite differences. We treat the cell averages $\bar{\rho}$ in the finite volume scheme as point values $\rho$ in the discretisation of the quantum pressure term. Technically, we use the MUSCL scheme as a conservative finite-difference method and not as a finite volume scheme. This has the consequence that the maximum order of accuracy we can reach with the linear subgrid MUSCL scheme is second order. The quantum pressure term is discretised as
$$\frac{\Delta \sqrt{\rho}}{\sqrt{\rho}} = \left(\frac{1}{2} \Delta \log(\rho) + \frac{1}{4} \left(\nabla \log(\rho)\right)^2\right)$$
with second-order central finite differences for both the gradient operator and the Laplacian.
The resulting semi-discrete scheme is discretised with a third-order RK method. A second-order discretisation also works well, but requires a more stringent CFL condition.
The combined scheme requires the CFL conditions
$$\Delta \tau \leq  \min\left[C_D \frac{m}{\hbar} \Delta x^2, C_K \Delta x \cdot m \left(2 \hbar \sum_{i=1}^3{|\partial_i S|}\right)^2\right]$$
where the constants $C_D$ (`DT__HYBRID_CFL`) and $C_K$ (`DT__HYBRID_VELOCITY`) were determined empirically.

## Stability
While the fluid scheme has been observed to be stable in fluid-only runs, hybrid runs pose additional challenges. When vortices form on wave levels, the phase field develops discontinuities and the quantum pressure becomes ill-defined. To handle this case, the fluid scheme implements the following additional checks.
- switch to a first-order scheme in time and space in fully refined patches
- switch to a second-order centered-in-space finite difference scheme wherever the fluid equation produces unphysical data.
This modification to the fluid scheme is work-in-progress and must be modified if it proves to be unstable in production runs.

# Known Issues in the Hybrid Scheme

## Bitwise Reproducibility

- **Issue**: Bitwise reproducibility currently fails in the hybrid scheme due to the conversion from RE/IM to DENS/PHAS when storing fields in HDF5.
- **Potential Solution**: Implement conversion from RE/IM <-> DENS/PHAS using high-precision routines to ensure bitwise identity for significant digits.

## Sudden Overrefinement

- **Issue**: Non-reproducible overrefinement, similar to what occurs in the heavy halo case, still occurs without clear cause. See [this Slack thread](https://calab-ntu.slack.com/archives/C02JV33GF5Y/p1695378046587619) for more information.

## Spurious Halos

### Causes
- Insufficient resolution in the wave region.
- Too large a value for `QP` in the QP refinement criterion.
- Use of monotonic interpolation which may not capture the complexity of wave functions.
- Interpolation in low-resolution regions.
- Phase interpolation errors, particularly in regions with complex interference patterns.

## Detecting Fluid Scheme Failures

- **Issue**: The fluid scheme may fail in low-density regions, especially around vortices where the wave scheme is used. These failures can propagate through levels via restriction operations.
- **Handling**: If the fluid solver fails (indicated by negative densities or NaN values), and a refined wave counterpart exists, overwrite the values using second-order finite-difference discretization. But some users still report instabilities, especially together with spectral interpolation


## Potential Improvements
- Global memory accesses in the ELBDM solver are uncoalesced, suggesting potential for optimization by implementing an effective 3D transpose on the GPU.
- Currently `ELBDM_FIRST_WAVE_LEVEL` is fixed. Ideally, the code should be able to adaptively increase the first wave level when the resolution turns out to be insufficient for the wave solver. This would require implementing a detection for insufficient resolution and converting patches back from the wave to the fluid representation.

## Zoom-In Simulation

- Ensure that transitions at the zoom-in boundary are handled correctly. Currently, the mismatch between fluid and wave schemes requires using the wave scheme outside of the zoom-in region.
