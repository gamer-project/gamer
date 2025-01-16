# `configure.py` options
- Must enable
   - [[--model=HYDRO/ELBDM | Installation:-Option-List#--model]]
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Cosmological constants
   | Parameter name | value                  |
   |---             |---                     |
   | A_INIT         | 0.00990099009900990099 |
   | OMEGA_M0       | 0.315823               |
   | HUBBLE0        | 0.6732117              |

2. Simulation
   | Parameter name       | value           |
   |---                   |---              |
   | BOX_SIZE             | 30.0 (Mpc/h)    |
   | NX0_TOT              | 128             |
   | MAX_LEVEL            | 5               |
   | NPAR                 | 2097152 (128^3) |
   | OPT__FLAG_NPAR_PATCH | 2               |
   | OPT__FREEZE_FLUID    | 1               |

3. Particle initial condition
   | Parameter name | value |
   |---             |---    |
   | PAR_IC_FORMAT  | 1     |
   | PAR_IC_MASS    | -1.0  |
   | PAR_IC_TYPE    | -1    |


# Note
1. CDM cosmological large-scale structure simulations
2. Fiducial `PAR_IC` file can be downloaded with the command.
   ```bash
   sh download_ic.sh
   ```
3. GAMER currently doesn't support particle-only simulations. So the following temporary solutions are adopted
   1. Set gas density/energy to arbitrarily small (for [[--model | Installation:-Option-List#--model]])
      or wave function to zero (for [[--model=ELBDM | Installation:-Option-List#--model]]) in `SetGridIC()`.
   2. Enable [[OPT__FREEZE_FLUID | Hydro#OPT__FREEZE_FLUID]]
4. Default maximum spatial resolution is 30.0/(128*2^5)~7.3 kpc/h (comoving).