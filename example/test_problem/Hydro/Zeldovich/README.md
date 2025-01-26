# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]=`HYDRO`
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gsl | Installation:-Option-List#--gsl]]
   - [[--fftw | Installation:-Option-List#--fftw]]
   - [[--dual | Installation:-Option-List#--dual]]
- Must disable
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]
- Configure.py
   ```
   python configure.py --model=HYDRO --gravity=true --comoving=true --particle=true --gsl=true --fftw=FFTW2 --dual=DE_ENPY
   ```


# Default setup:
1. Cosmological constants
   | Parameter name | value                  | Note |
   |---             |---                     |--- |
   | A_INIT         | 0.00990099009900990099 | z_i = 100 |
   | OMEGA_M0       | 1.0                    | matter-only Universe |
   | HUBBLE0        | 0.6732117              |  |

2. Simulation
   | Parameter name     | value       |
   |---                 |---          |
   | BOX_SIZE           |  80 (Mpc/h) |
   | NX0_TOT_X/Y/Z      |  128/32/32  |
   | OPT__BC_*          |  1          |
   | MAX_LEVEL          |  0          |
   | GPU_NSTREAM        | -1          |
   | OPT__FLAG_RHO      |  1          |
   | OPT__OUTPUT_PART   |  4          |
   | OPT__OUTPUT_BASEPS |  0          |
   | OUTPUT_PART_X/Y/Z  |  0.0        |

3. Particle
   | Parameter name  | value |
   |---              |---    |
   | PAR_NPAR        | -1    |
   | PAR_INIT        | 1     |
   | PAR_IC_FORMAT   | 1     |
   | PAR_IC_MASS     | -1.0  |
   | PAR_IC_TYPE     | -1    |

For "gas-only" setup:
4-I. `Input__Parameter`
   | Parameter name     | value        |
   |---                 |---           |
   | GAMMA              | 1.6666666667 |
   | MOLECULAR_WEIGHT   | 0.6          |
   | DUAL_ENERGY_SWITCH | 3.0e-2       |
   | OPT__FREEZE_FLUID  | 0            |
   | OPT__OUTPUT_USER   | 1            |
4-II. `Input__TestProb`
   | Parameter name  | value |
   |---              |---    |
   | Gas_Par_Setup   | 1     |
   | n_Pert_Wave_Len | 1     |
   | zc_Collapse     | 1.0   |
   | GasTemp_Init    | 100   |

For "particle-only" setup:
4-I. `Input__Parameter`
   | Parameter name    | value |
   |---                |---    |
   | OPT__FREEZE_FLUID | 1     |
   | OPT__OUTPUT_USER  | 0     |
4-II. `Input__TestProb`
   | Parameter name  | value |
   |---              |---    |
   | Gas_Par_Setup   | 2     |
   | n_Pert_Wave_Len | 1     |
   | zc_Collapse     | 1.0   |
   | NPar_X          | 64    |


# Note
1. Zeldovich pancake collapse cosmological simulations

2. It is recommended to enable [[--double | Installation:-Option-List#--double]] for gas-only setup
   (`Gas_Par_Setup = 1`) to properly resolve thepressure and temperature of the cold flow

3. Since GAMER does not support genuine 1D/2D simulations, one should set the root-level grid size to 32
   (i.e., patch size*2*2) in the transverse directions. The density distribution is perturbed along only
   the x-axis and remains uniform across the y-z plane.

4. The perturbation wavelength set by `n_Pert_Wave_Len`, initial simulation scale factor
   [[A_INIT | Runtime-Parameters:-Cosmology#A_INIT]], pancake collapse redshift `zc_Collapse` are free
   parameters in this test problem. Note that `zc_Collapse` has to be smaller than the initial redsfhit.
   Analytical solutions are applicable for `z > zc_Collapse`.

5. The parameter [[PAR_NPAR | Runtime-Parameters:-Particles#PAR_NPAR]] is irrelevant for the gas-only setup, and replaced by
   `NPar_X` in `Input__TestProb` for the particle-only setup. Freezing the fluid component
   [[OPT__FREEZE_FLUID | Runtime-Parameters:-Hydro#OPT__FREEZE_FLUID]] = 1 is necessary for the particle-only setup. For the
   particle-only setup, `NPar_X = 64` is the minimum spatial resolution such that the collapsed structure
   can be properly resolved.

6. The analysis script `./plots/analysis_plot_Zeldovich.py` works for both gas-/particle-only setups. Note
   that the parameters under "# user-specified parameters" in the script NEED to match those adopted in
   `Input__TestProb`, especially the value of `Gas_Par_Setup`.

7. For references, see Sec. 9.2 of [Springel, MNRAS 401 791 (2010)](https://doi.org/10.1111/j.1365-2966.2009.15715.x)
   and Sec. 3 of [Heitmann et al., ApJS 160 28 (2005)](https://dx.doi.org/10.1086/432646).