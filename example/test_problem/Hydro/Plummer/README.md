# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`HYDRO`
  - [[--gravity | Installation:-Option-List#--gravity]]
- Must disable
  - [[--comoving | Installation:-Option-List#--comoving]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
None


# Note
- Test the evolution of the Plummer model
  - `Plummer_Collision` = `0`: test the stability of a single Plummer cloud
  - `Plummer_Collision` = `1`: test the collision of two Plummer clouds

- One can set [[--passive | Installation:-Option-List#--passive]]=`2`
  and `Plummer_AddColor`=`1` for `Plummer_Collision`=`1`
  - Different clouds will be assigned with different passive scalars

- Plummer cloud(s) can have both gas and particles
  - `Plummer_GasMFrac`: gas mass fraction

- To enable feedback
  1. Enable [[FB_USER | Runtime-Parameters:-Feedback#FB_USER]]
  2. Enable `Plummer_FB_Exp` and/or `Plummer_FB_Acc`

- To test mass conservation to the machine precision with mass accretion feedback
  1. Set `Plummer_FB_ExpMMax`=`0.0` (no mass loss in explosion feedback)
  2. Set [[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]]=`1` (periodic BC)
  3. Use the gnuplot script `plot_script/plot_mass.gpt`
