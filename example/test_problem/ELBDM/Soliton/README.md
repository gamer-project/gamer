# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- [[NX0_TOT_X/Y/Z | Runtime-Parameters:-General#NX0_TOT_X]]=`96`
- [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]=`2`


# Note
- `Soliton_N`=`1`: test a single soliton
  - Set `Soliton_RSeed`<`0`
  - Set `Soliton_EmptyRegion`=`0.0`

- `Soliton_N`>`1`: soliton merger
  - Set `Soliton_RSeed`>=`0`
  - Set `Soliton_EmptyRegion`>`0.0` (e.g., `1.0e2`)
