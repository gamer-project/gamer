> [!CAUTION]
> Please do not edit this file(page) manually since the workflow will overwrite your changes.
>
> This file(page) is automatically generated by the workflow `Update test problem wiki page` using the script `tool/wiki/sync_test_problem_pages.py`.
>
> The workflow is triggered by push changes to any of `example/test_problem/*/*/README.md` and `tool/wiki/sync_test_problem_pages.py`.


# `configure.py` options
- Must enable
   - [[--model=ELBDM | Installation:-Option-List#--model]]
   - [[--double | Installation:-Option-List#--double]]
   - [[--passive=1 | Installation:-Option-List#--passive]]
- Must disable
   - [[--gravity | Installation:-Option-List#--gravity]]
   - [[--particle | Installation:-Option-List#--particle]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Evolve the plane wave for six periods
2. Apply the periodic BC
   --> Set [[OPT__BC_FLU_* | Hydro#OPT__BC_FLU_XM]] = 1


# Note
1. Only support 1D --> Use `PWave_XYZ` to control the propagation direction