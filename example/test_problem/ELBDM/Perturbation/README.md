# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--model]]
   - [[--elbdm_scheme=ELBDM_HYBRID | Installation:-Option-List#--elbdm_scheme]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Study gravitational collapse in one, two and three dimensions (2D by default)
2. Use periodic boundary conditions
3. Demonstrate refinement in hybrid scheme by using `Input__Flag_Interference`
4. Compared to JeansInstabilityComoving, this test also supports 1D, 2D and 3D collapse

# Note
1. Evolve small-amplitude plane waves on homogeneous background density
2. `Perturbation_Amplitude` is maximum amplitude of individual perturbations, `Perturbation_BgAmplitude` is amplitude of background density, `Perturbation_N` is number of Fourier modes in perturbation (maximum is 4)
3. Dimensionality of the problem can be set via `Perturbation_NDim` (set amplitude of plane waves in other dimensions to zero)
3. For reference simulation using base-level spectral method, use `MODEL = WAVE` and `Input__Parameter_BaseSpectral` and `generate_make_BaseSpectral.sh`
