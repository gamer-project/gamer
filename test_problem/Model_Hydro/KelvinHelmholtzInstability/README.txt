
   *********************************************
   ** HYDRO Kelvin-Helmholtz instability test **
   *********************************************

Note:
-------------------------------------------------------------------------------
1. Adopt the periodic boundary condition everywhere

2. Shear velocity is defined on the x-y plane

3. Random values are added to velocity in all three directions

4. When KH_AllRankSame == true, all MPI ranks will assign the same initial condition
   --> suitable for the weak scaling test

5. This test problem overwrite the original "Hydro_Init_StartOver_AssignData.cpp" file
   to remove the OpenMP pragma since rand() is not thread-safe
