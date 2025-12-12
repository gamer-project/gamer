# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
- Must disable
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Evolve vortex pair for one period
- Use the periodic boundary conditions
- Wavelengths along x and y can be adjusted to multiples of the box size
- Wave in z-direction is turned off by default


# Note
- Evolve vortex pair with linear motion along x in a 2D simulation
  - Wave function `psi_vorpair(x,y) = background + A*cos(ky*y)*exp( i*(kx*x-Omega*t+Phase0) )`,
    where `A` is a constant on the order of background, `kx` and `ky` are wavenumbers,
    `Omega=0.5/ELBDM_ETA*(kx^2+ky^2)`, and `Phase0` is a phase constant
- Optionally: Add wave in `z` direction
  - `psi(x, y, z) = psi_vorpair(x,y) + background_z * exp( i*(kz*z-ZWaveOmega*t) )`,
    where `background_z` is a constant on the order of background, `kz` is a wavenumber and
    `Omega=0.5/ELBDM_ETA*(kz^2)`
