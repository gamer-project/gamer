Timetable for dumping data when adopting
[[OPT__OUTPUT_MODE | Runtime-Parameters:-Outputs#OPT__OUTPUT_MODE]]=3.
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
