All runtime parameters are specified in the input files `Input__*`:

* [Input__Parameter](#input__parameter): **MOST of the runtime parameters are put here**
* [Input__TestProb](#input__testprob): problem-specific parameters
* [Input__Flag_*](#input__flag_): grid refinement criteria
* [Input__DumpTable](#input__dumptable): timetable for dumping data
* [Input__MemoryPool](#input__memorypool): numbers of patches to be preallocated for the memory pool
* [Input__Note](#input__note): simulation notes

Input file templates are put in `example/input/`.


## Input__Parameter

* [Parameter List](#parameter-list)
* [File Format](#file-format)
* [Adding New Parameters](#adding-new-parameters)


### Parameter List
* [[All]]
-- List all parameters in alphabetical order

* [[General | Runtime-Parameters:-General]]
-- Parameters applicable to all simulations, such as simulation domain, root-level grid, simulation end time...

* [[MPI and OpenMP | MPI-and-OpenMP#runtime-parameters]]
-- Number of MPI processes and OpenMP threads, load-balancing parameters

* [[GPU | GPU#runtime-parameters]]
-- GPU IDs and optimization parameters

* [[Units | Runtime-Parameters:-Units]]
-- Unit system

* [[Initial Conditions | Initial-Conditions]]
-- Initialization methods, restarting simulations

* [[Hydro | Hydro#runtime-parameters]]
-- Hydro solvers, physical constants, boundary conditions

* [[Gravity | Gravity#runtime-parameters]]
-- Gravity solvers, physical constants, boundary conditions

* [[Particles | Particles#runtime-parameters]]
-- Particle parameters

* [[Cosmology | Runtime-Parameters:-Cosmology]]
-- Cosmological parameters

* [[Chemistry and Radiation | Chemistry-and-Radiation#runtime-parameters]]
-- GRACKLE parameters

* [[Star Formation | Star-Formation#runtime-parameters]]
-- Star formation parameters

* [[Feedback | Feedback#runtime-parameters]]
-- Feedback parameters

* [[Timestep | Runtime-Parameters:-Timestep]]
-- Timestep criteria

* [[Refinement | Runtime-Parameters:-Refinement]]
-- Grid refinement criteria

* [[Interpolation | Runtime-Parameters:-Interpolation]]
-- Interpolation schemes

* [[Outputs | Outputs#runtime-parameters]]
-- Data output parameters

* [[Miscellaneous | Runtime-Parameters:-Miscellaneous]]
-- Miscellaneous parameters such as timing options, log files, and self-checking options...


### File Format
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


### Adding New Parameters

See [[ Adding Parameters|Adding-Parameters ]].


## Input__TestProb

Problem-specific runtime parameters. This file follows the same
syntax as [Input__Parameter](#file-format). See
`example/test_problem/*/*/Input__TestProb` for some examples.
The following shows
`example/test_problem/Hydro/Riemann/Input__TestProb`:

```
Riemann_Prob     0     # target Riemann problem (0 : Sod's shock tube
                       #                         1 : strong shock
                       #                         2 : two shocks
                       #                         3 : Einfeldt's 1-2-0-3
                       #                         4 : Einfeldt's 1-1-2-5
                       #                         5 : sonic rarefaction wave
Riemann_LR      +1     # wave propagation direction (>0/<0 --> positive/negative direction) [1]
Riemann_XYZ      0     # wave propagation axis (0/1/2 --> x/y/z) [0]
```


## Input__Flag_*

Refinement thresholds of various refinement criteria `OPT__FLAG_*`
described in [[Runtime Parameters -- Refinement | Runtime-Parameters:-Refinement]],
where each criterion has its own input file.
For instance, the following shows an example of the table `Input__Flag_Rho`
used by [[OPT__FLAG_RHO | Runtime-Parameters:-Refinement#OPT__FLAG_RHO]]:

```
# Level                         Density
      0                             8.0
      1                            64.0
      2                           512.0
```

In this example, cells on levels 0, 1, and 2 will be flagged for
refinement if their gas mass densities exceed 8.0, 64.0, and 512.0,
respectively.

To add new flag tables, edit `src/Init/Init_Load_FlagCriteria.cpp`.

Table format:

* Must have one and only one header line.

* The first column (i.e., levels 0, 1, 2 in the example above)
is actually useless and will not be loaded at all. The refinement
thresholds and other parameters (if any) must be put in column(s)
other than the first column.

* Empty and comment lines (i.e., lines starting with #) are NOT allowed
except in the first header line.

* Must contain at least [[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]
data lines excluding the header line.


## Input__DumpTable

Timetable for dumping data when adopting
[[OPT__OUTPUT_MODE | Outputs#OPT__OUTPUT_MODE]]=3.
The following example will dump data at t=1.0, 2.3, and 3.7.

```
#Dump ID             Dump Time
       0                   1.0
       1                   2.3
       2                   3.7
***************END LINE***************
```

Note that this text file currently has a quite inflexible format:

* Must have one and only one header line.

* The first column (i.e., the dump IDs 0, 1, 2 in the example above)
is actually useless and will not be loaded at all. The output timetable
must be put in the second column.

* Empty and comment lines (i.e., lines starting with #) are NOT allowed
except in the first header line.

* Timetable must be in ascending order.

* The loading routine will stop loading when detecting a line
starting with *. All remaining lines will be ignored.

* Simulation end time [[END_T | Runtime-Parameters:-General#END_T]]
loaded from `Input__Parameter` will be reset to the maximum output
time in the dump table if the latter is found to be smaller.


## Input__MemoryPool

Number of patches to be preallocated on each AMR level when enabling
[[OPT__MEMORY_POOL | Runtime-Parameters:-Refinement#OPT__MEMORY_POOL]].
The following example will preallocate 100, 800 and 6400 patches on
levels 0, 1, and 2, respectively, and will not preallocate any patch
above level 2.

```
# Level      Number of patches to be preallocated
      0                                       100
      1                                       800
      2                                      6400
```

Table format:

* Must have one and only one header line.

* The first column (i.e., levels 0, 1, 2 in the example above)
is actually useless and will not be loaded at all. The numbers of
patches to be preallocated must be put in the second column.

* Empty and comment lines (i.e., lines starting with #) are NOT allowed
except in the first header line.

* At most
[[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]+1
lines will be loaded. If there are fewer than
[[MAX_LEVEL | Runtime-Parameters:-Refinement#MAX_LEVEL]]+1
lines in the input file, no patches will be preallocated on the
unspecified levels.


## Input__Note

See [[Taking Notes | Running-the-Code#taking-notes]].
