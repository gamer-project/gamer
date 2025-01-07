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
