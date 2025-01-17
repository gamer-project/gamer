# `configure.py` options
Enable  : MODEL=HYDRO, SRHD, NCOMP_PASSIVE_USER = 4
Disable : GRAVITY, MHD, COMOVING, PARTICLE
- Must enable
   - [[--model | Installation:-Option-List#--model]]
   - [[--srhd | Installation:-Option-List#--srhd]]
   - [[--passive=4 | Installation:-Option-List#--passive]]
- Must disable
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--mhd | Installation:-Option-List#--mhd]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Units:
   | [L]   | [M]        | [V] | [T] = [L]/[V] | [D] = [M]/[L]^3    |
   |---    |---         |---  |---            |---                 |
   | 5 kpc | 1.0e5 Msun | 1 c | ~ 16.308 kyr  | ~ 5.41e-29 g/cm**3 |

2. Default resolution ~ 39 pc ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]= 3)


# Note
1. [ZuHone et al. 2023, Galaxies, vol. 11, issue 2, p. 51](https://doi.org/10.3390/galaxies11020051)
