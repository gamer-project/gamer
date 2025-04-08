This page includes the following topics:
* [Parameter List](#parameter-list)
* [File Format](#file-format)
* [Adding New Parameters](#adding-new-parameters)


## Parameter List
* [[All | Runtime-Parameters:-All]]
-- List all parameters in alphabetical order

* [[General | Runtime-Parameters:-General]]
-- Parameters applicable to all simulations, such as simulation domain, root-level grid, simulation end time...

* [[MPI and OpenMP | Runtime-Parameters:-MPI-and-OpenMP]]
-- Number of MPI processes and OpenMP threads, load-balancing parameters

* [[GPU | Runtime-Parameters:-GPU]]
-- GPU IDs and optimization parameters

* [[Units | Runtime-Parameters:-Units]]
-- Unit system

* [[Initial Conditions | Runtime-Parameters:-Initial-Conditions]]
-- Initialization methods, restarting simulations

* [[Hydro | Runtime-Parameters:-Hydro]]
-- Hydro solvers, physical constants, boundary conditions

* [[Gravity | Runtime-Parameters:-Gravity]]
-- Gravity solvers, physical constants, boundary conditions

* [[Particles | Runtime-Parameters:-Particles]]
-- Particle parameters

* [[Cosmology | Runtime-Parameters:-Cosmology]]
-- Cosmological parameters

* [[Chemistry and Radiation | Runtime-Parameters:-Chemistry-and-Radiation]]
-- GRACKLE parameters

* [[Star Formation | Runtime-Parameters:-Star-Formation]]
-- Star formation parameters

* [[Feedback | Runtime-Parameters:-Feedback]]
-- Feedback parameters

* [[Timestep | Runtime-Parameters:-Timestep]]
-- Timestep criteria

* [[Refinement | Runtime-Parameters:-Refinement]]
-- Grid refinement criteria

* [[Interpolation | Runtime-Parameters:-Interpolation]]
-- Interpolation schemes

* [[Outputs | Runtime-Parameters:-Outputs]]
-- Data output parameters

* [[Miscellaneous | Runtime-Parameters:-Miscellaneous]]
-- Miscellaneous parameters such as timing options, log files, and self-checking options...


## File Format
All parameters in this file follow the syntax:

```
NAME     VALUE     # COMMENTS
```

* Comment symbol: #.

* Empty and comment lines (i.e., lines starting with #) are ignored.

* Parameter defaults are described as [DEFAULT] in the commentary.
For example,

    ```
    REGRID_COUNT     1     # refine every REGRID_COUNT sub-steps [4]
    ```

    indicates that the default value of the parameter `REGRID_COUNT` is 4.

* Parameters set to "auto" (usually by assigning a negative value) do
not have deterministic defaults and will be reset later according to
the adopted compilation options and/or other runtime parameters. For
example,

    ```
    DT__FLUID     -1.0     # hydro CFL factor (<0.0=auto) [-1.0]
    ```

    indicates that if the parameter `DT__FLUID` is assigned with a
negative value (which is also the default value as suggested by `[-1.0]`),
it will be reset later depending on other configurations
(in this case, it's value depends on the adopted hydro scheme).

* For boolean options, 0=off and 1=on.

* All dimensional variables must be set consistently with the adopted unit
system (set by [[OPT__UNIT | Runtime-Parameters:-Units#OPT__UNIT]]) unless otherwise specified.

* Parameter names are case sensitive.

* The program will display warning messages to stderr when detecting
duplicate or unrecognizable parameters. For example,

    ```
    WARNING : unrecognizable parameter [DUAL_ENERGY_SWITCH       ] at line  113 !!
    WARNING : duplicate      parameter [OPT__VERBOSE             ] at line  160 !!
    ```

    It will display stdout messages when setting any parameters
to defaults. For example

    ```
    NOTE : parameter [UNIT_L                   ] is set to the default value [-1.00000000000000e+00]
    NOTE : parameter [UNIT_M                   ] is set to the default value [-1.00000000000000e+00]
    ```

    It will also display warning messages to stderr after resetting
any parameters (usually because they are either useless or set to
"auto" as described above). For example,

    ```
    WARNING : parameter [DT__FLUID                   ] is reset to [ 5.00000000000000e-01]
    WARNING : parameter [GPU_NSTREAM                 ] is reset to [ 1                   ] since GPU is disabled
    WARNING : parameter [OPT__NORMALIZE_PASSIVE      ] is reset to [ 0                   ] since there are no passive scalars
    ```


## Adding New Parameters

See [[Adding Parameters | Adding-Parameters]].
