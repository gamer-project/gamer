# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--double | Installation:-Option-List#--double]]
  - [[--passive | Installation:-Option-List#--passive]]=`1`
- Must disable
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--particle | Installation:-Option-List#--particle]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- Evolve the plane wave for six periods
- Apply the periodic BC
  - Set [[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]]=`1`


# Note
- Only support 1D
  - Use `PWave_XYZ` to control the propagation direction
- A passive field named `WrappedPhase` is output in the data
- There are plotting scripts provided
  - For [[--elbdm_scheme | Installation:-Option-List#--elbdm_scheme]]=`WAVE`: plot_Phase.gpt and plot_WaveFunc.gpt
  - For [[--elbdm_scheme | Installation:-Option-List#--elbdm_scheme]]=`HYBRID`: plot_Phase_hybrid.gpt and plot_WaveFunc_hybrid.gpt
- Support using the analytical-solution BC ([[OPT__BC_FLU_* | Runtime-Parameters:-Hydro#OPT__BC_FLU_XM]]=`4`) and it must be used when:
  - `PWave_NWavelength` is not a integer
  - Using the hybrid scheme to run the case of travelling wave and there are only fluid levels at the boundary
- For [[--elbdm_scheme | Installation:-Option-List#--elbdm_scheme]]=`HYBRID`
  - The default setup is without refinement and adopts the fluid scheme;
    When setting [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]>`0`,
    the refinement is based on [[OPT__FLAG_INTERFERENCE | Runtime-Parameters:-Refinement#OPT__FLAG_INTERFERENCE]]
    and the refined levels adopt the wave scheme while the base level adopts the fluid scheme
