# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--comoving | Installation:-Option-List#--comoving]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Physical Parameters
  | Parameter name | Value       | Note   |
  |----------------|-------------|--------|
  | BOX_SIZE       | 1.4         | Mpc/h  |
  | ELBDM_MASS     | 8.0e-23     | eV/c^2 |
  | A_INIT         | 3.124024e-4 |        |
  | OMEGA_M0       | 0.2835      |        |
  | HUBBLE0        | 0.6955      |        |

- AMR
  | Parameter name        | Value | Note  |
  |-----------------------|-------|-------|
  | MAX_LEVEL             | 3     |       |
  | OPT__FLAG_RHO         | 1     |       |
  | OPT__FLAG_LOHNER_DENS | 1     |       |

- Flag tables
  - `Input__Flag_Rho`:    `example/input/Input__Flag_Rho8`
  - `Input__Flag_Lohner`: `example/input/Input__Flag_Lohner__FLASH1_LSS_DH3.0`


# Note
- Cosmological large-scale structure simulations

- Download the IC file
  ```bash
  sh download_ic.sh
  ```

- Default maximum resolution is low
  - Only `1.4/(256*2^3) ~ 0.7` kpc/h (comoving)
  - Can only marginally resolve a central soliton

- Some yt visualization scripts are put in `plot_script`
