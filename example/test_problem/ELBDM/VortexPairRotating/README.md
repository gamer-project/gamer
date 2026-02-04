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
- Use the analytical solution as the boundary conditions


# Note
- Evolve a rotating vortex pair in 2D
  - Wave function `psi(R,phi) = background - A*J1( sqrt(2*Eta*Omega)*R )*exp( i*(phi-Omega*t+Phase0) )`,
    where `A` is a constant on the order of background, `Eta=ELBDM_MASS/PLANCK_CONSTANT`,
    `phi` is azimuthal angle, `R` is radius, `Omega` is angular frequency, and `Phase0` is a phase constant

- Ref: [Tzihong Chiueh et al 2011 J. Phys. B: At. Mol. Opt. Phys. 44 115101](https://doi.org/10.1088/0953-4075/44/11/115101),
  Vortex turbulence in linear Schrödinger wave mechanics
