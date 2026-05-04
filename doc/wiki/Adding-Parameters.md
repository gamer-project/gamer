TBF.

Mandatory steps are marked by &#x1F4CC;.

1. &#x1F4CC; Add new global variables to both `src/GAMER/Main.cpp` and `include/Global.h`

> [!IMPORTANT]
> The length of the string variable must be `MAX_STRING`.

2. &#x1F4CC; Edit `src/Init/Init_Load_Parameter.cpp` to load new parameters

3. Add notes in `src/Auxiliary/Aux_TakeNote.cpp` [optional but recommended]

4. Add checks in `src/Auxiliary/Aux_Check.cpp` [optional but recommended]

5. Record new parameter in the HDF5 snapshots [optional but recommended]

    1. Add new parameters to `include/HDF5_Typedef.h`

    2. Add new parameters to `src/Output/Output_DumpData_Total_HDF5.cpp`

    3. Check new parameters during restart by editing `src/Init/Init_Restart_HDF5.cpp`
