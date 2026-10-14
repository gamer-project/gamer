This page describes various grid refinement criteria and interpolation schemes.


## Enabling AMR

It only takes three steps to enable AMR:

* Set [[MAX_LEVEL | [Runtime-Parameters]-Refinement#MAX_LEVEL]]
* Turn on at least one of the refinement criteria `OPT__FLAG_*`
* Edit the corresponding input file(s)
[[Input__Flag_{} | [Runtime-Parameters]-Input__Flag_{}]]
to specify the refinement thresholds

For details, see the descriptions of various refinement criteria `OPT__FLAG_*`
in [[ [Runtime parameters] Refinement | [Runtime-Parameters]-Refinement ]].


## Compilation Options

Related options:
[[--nlevel | [Installation]-Option-List#--nlevel]], &nbsp;
[[--max_patch | [Installation]-Option-List#--max_patch]] &nbsp;


## Runtime Parameters
[[ [Runtime parameters] Refinement | [Runtime-Parameters]-Refinement ]] \
[[ [Runtime parameters] Interpolation | [Runtime-Parameters]-Interpolation ]]

Other related parameters:
none


## Remarks

### Potential outside the isolated boundaries
When adopting the isolated boundary conditions for gravity (i.e.,
[[OPT__BC_POT | [Runtime-Parameters]-Gravity#OPT__BC_POT]]=2), the ghost zones of
gravitational potential outside the simulation domain are currently
filled out by extrapolation.


<br>

## Links
* [[Main page of Runtime Parameters | Runtime-Parameters]]
