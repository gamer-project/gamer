> [!CAUTION]
> Please do not edit this file(page) manually since the workflow will overwrite your changes.
>
> This file(page) is automatically generated by the workflow `Update test problem wiki page` using the script `tool/wiki/sync_test_problem_pages.py`.
>
> The workflow is triggered by push changes to any of `example/test_problem/*/*/README.md` and `tool/wiki/sync_test_problem_pages.py`.


# `configure.py` options
- Must enable
   - [[--model | Installation:-Option-List#--model]]
- Must disable
   - [[--comoving | Installation:-Option-List#--comoving]]
   - [[--particle | Installation:-Option-List#--particle]]
   - [[--gravity | Installation:-Option-List#--gravity]]
- Available options
   - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
1. Adopt the Lohner's error estimator of both mass density and pressure as the refinement criteria
2. Refinement threshold of the Lohner's error estimator = 0.80
3. Maximum refinement level ([[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]) = 2


# Note
1. Ref: [Ryu, D., Ostriker, J. P., Kang, H., & Cen, R. 1993, ApJ, 414, 1](https://doi.org/10.1086/173051)
2. This test is good for testing the dual-energy formalism
   --> Without it the preshock region will be over-heated
