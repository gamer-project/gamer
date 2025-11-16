# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--elbdm_scheme | Installation:-Option-List#--elbdm_scheme]]=`HYBRID`
  - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Study gravitational collapse in one, two and three dimensions (2D by default)
- Use periodic boundary conditions
- Demonstrate refinement in hybrid scheme by using `Input__Flag_Interference`
- Compared to JeansInstabilityComoving, this test also supports 1D, 2D and 3D collapse

# Note
- Evolve small-amplitude plane waves on homogeneous background density

- `Perturbation_Amplitude` is maximum amplitude of individual perturbations,
  `Perturbation_BgAmplitude` is amplitude of background density,
  `Perturbation_N` is number of Fourier modes in perturbation (maximum is 4)

- Dimensionality of the problem can be set via `Perturbation_NDim` (set amplitude of plane waves in other dimensions to zero)

- For reference simulation using base-level spectral method, use [[--elbdm_scheme | Installation:-Option-List#--elbdm_scheme]]=`WAVE` and `Input__Parameter_BaseSpectral` and `generate_make_BaseSpectral.sh`
