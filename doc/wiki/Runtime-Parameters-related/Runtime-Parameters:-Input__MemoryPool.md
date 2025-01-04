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
