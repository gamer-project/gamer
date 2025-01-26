# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]=`ELBDM`
   - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Code units
   - `UNIT_L` = Mpc/h, where h=0.6955 is the present dimensionless Hubble parameter
   - `UNIT_V` = 100 km/s
   - `UNIT_D` = `rho_bg` (background matter density at z=0)
      - Mass density and wavefunction are normalized to `rho_bg`

2. ELBDM
   | Parameter name         | value   |
   |---                     |---      |
   | ELBDM_MASS             | 8.0e-23 |
   | ELBDM_REMOVE_MOTION_CM | 1       |
   | ELBDM_TAYLOR3_AUTO     | 0       |

3. Boundary conditions
   | Parameter name | value |
   |---             |---    |
   | OPT__BC_FLU_*  | 1     |
   | OPT__BC_POT    | 2     |

4. AMR
   | Parameter name | value |
   |---             |---    |
   | MAX_LEVEL      | 0     |


# Libyt covering_grid setup
1. Code units
   - `UNIT_L` = Mpc/h, where h=0.6955 is the present dimensionless Hubble parameter
   - `UNIT_V` = 100 km/s
   - `UNIT_D` = rho_bg (background matter density at z=0)
      - Mass density and wavefunction are normalized to `rho_bg`

2. ELBDM
   | Parameter name         | value   |
   |---                     |---      |
   | ELBDM_MASS             | 8.0e-23 |
   | ELBDM_REMOVE_MOTION_CM | 1       |
   | ELBDM_TAYLOR3_AUTO     | 0       |

3. Boundary conditions
   | Parameter name | value |
   |---             |---    |
   | OPT__BC_FLU_*  | 1     |
   | OPT__BC_POT    | 2     |

4. AMR
   | Parameter name | value |
   |---             |---    |
   | MAX_LEVEL      | 0     |
   | OPT__FLAG_RHO  | 1     |

5. End simulation condition
   | Parameter name | value         |
   |---             |---            |
   | END_T          | 5.7116620e-03 |
   | END_STEP       | 32            |

6. libyt
   | Parameter name | value                                    |
   |---             |---                                       |
   | YT_SCRIPT      | libyt_script/inline_script_covering_grid |
   | YT_VERBOSE     | 1                                        |


# Note
1. Download the IC file
   ```bash
   sh download_ic.sh
   ```

2. Some examples of yt visualization scripts are put in `yt_script`

3. Simulate a single isolated halo extracted from a cosmological simulation

4. About libyt `covering_grid` test
   1. Use submit script `./libyt_script/submit_gamer.job` for job submission

   2. Put `./libyt_script/inline_script_covering_grid.py` under the same folder as gamer

   3. For determining the `left_edge` and `dims` in function `ds.covering_grid`:
      - `left_edge`: LV1 resolution is 0.175/512/2 ; region covered by LV1 box (by `ds.covering_grid`) is 0.175/512/2*512; 0.04375 = (0.175 - 0.175/512/2*512)/2
      - `dims`:      Plan to cover region with half of the simulation box length, i.e. will have 256X256X256 base level cells -> refine to MAX_LEVEL=1 -> LV1 cells is 512X512X512

   4. Use `plot_slice-dens_covering_grid.py` to extract density slices from `.npz` files

   5. Use `make_movie.sh` to convert `.png` pictures into `.mp4`