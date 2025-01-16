> [!CAUTION]
> Please do not edit this file(page) manually since the workflow will overwrite your changes.
>
> This file(page) is automatically generated by the workflow `Update test problem wiki page` using the script `tool/wiki/sync_test_problem_pages.py`.
>
> The workflow is triggered by push changes to any of `example/test_problem/*/*/README.md` and `tool/wiki/sync_test_problem_pages.py`.


# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--model]]
   - [[--elbdm_scheme=ELBDM_HYBRID | Installation:-Option-List#--elbdm_scheme]]
- Must disable
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--particle | Installation:-Option-List#--particle]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Evolve vortex pair for one period
2. Use the periodic boundary conditions
3. Wavelengths along x and y can be adjusted to multiples of the box size
4. Wave in z-direction is turned off by default


# Note
1. Evolve vortex pair with linear motion along x in a 2D simulation
   --> Wave function `psi_vorpair(x,y) = background + A*cos(ky*y)*exp( i*(kx*x-Omega*t+Phase0) )`
       where `A` is a constant on the order of background, `kx` and `ky` are wavenumbers,
       `Omega=0.5/ELBDM_ETA*(kx^2+ky^2)`, and `Phase0` is a phase constant
2. Optionally: Add wave in `z` direction
    -->   `psi(x, y, z) = psi_vorpair(x,y) + background_z * exp( i*(kz*z-ZWaveOmega*t) )`
       where `background_z` is a constant on the order of background, `kz` is a wavenumber and
       `Omega=0.5/ELBDM_ETA*(kz^2)`
3. If [[--conserve_mass | Installation:-Option-List#--conserve_mass]] is enabled, mass ought to be conserved in this test because of periodic boundary conditions