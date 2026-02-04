# `configure.py` options
- Must enable
  - [[--model | Installation:-Option-List#--model]]=`ELBDM`
  - [[--gravity | Installation:-Option-List#--gravity]]
  - [[--comoving | Installation:-Option-List#--comoving]]
- Must disable
  - [[--particle | Installation:-Option-List#--particle]]
  - [[--unsplit_gravity | Installation:-Option-List#--unsplit_gravity]]
- Available options
  - [[Miscellaneous Options | Installation:-Option-List#miscellaneous-options]]


# Default setup
- 32^3 uniform resolution
  - No grid refinement is adopted in this test

- [[DT__MAX_DELTA_A | Runtime-Parameters:-Timestep#DT__MAX_DELTA_A]] is set to an extremely large number so that the time-step will
  be controlled by the fluid solver

- The default setup in `Input__Parameter` and `Input__TestProb` gives an UNSTABLE
  solution, in which the solution grows exponentially. To test a STABLE solution,
  one can try the following parameters:

  - `Input__Parameter`
    | Parameter name | Value  |
    |----------------|--------|
    | BOX_SIZE       | 0.5    |
    | A_INIT         | 1.0e-5 |
    | OUTPUT_DT      | 2.0e-6 |

  - `Input__TestProb`
    | Parameter name | Value  |
    |----------------|--------|
    | Jeans_RealAmp0 | 4.0e-8 |

> [!CAUTION]
> This setup requires further verification


# Note
- Analytical solution reference: [Woo, T. & Chiueh, T. 2009, ApJ, 697, 850](https://doi.org/10.1088/0004-637X/697/1/850)

- Note that the imaginary part `I` grows much faster than the real part `R`.
  Consequently, the accuracy of linear prediction will be deteriorated
  when the assumption (`2R >> I^2`) starts to break down, especially in the
  higher-resolution tests, and hence the 2nd-order accuracy may no longer
  hold.

> [!CAUTION]
> 2nd-order accuracy has NOT been verified for the UNSTABLE solution
