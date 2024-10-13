
GAMER_ExtractUniform : construct the uniform-mesh data from the GAMER output file

==============================================================================

Version 1.4.0     09/24/2013
----------------------------
1. Output velocity field in ELBDM


Version 1.3.1     10/10/2012
----------------------------
#%!^%$!%!$#$

Version 1.3.0     10/04/2012
----------------------------
#%!^%$!%!$#$

Version 1.2.1     12/17/2011
----------------------------
#%!^%$!%!$#$

Version 1.2.0     12/02/2011
----------------------------
#%!^%$!%!$#$



******************
**     NOTE     **
******************
1. This tool can be used to cut a slice or 3D rectangle from the GAMER output
   file. User can specify the targeted refinement level by the option "-l Lv".
   Regions at levels lower than Lv will be refined to Lv, and regions at
   levels higher than Lv will be averaged to Lv. In other words, the output
   file of GAMER_ExtractUniform will contain a uniform-mesh data with the spatial
   resolution equal to the refinement level Lv.

   User can specify the range of the targeted slice or 3D rectangle by the
   options

      "-x/y/z starting coordinate" and "-X/Y/Z size".

   For example, if one specify

      "-n 7 -x 10 -y 20 -z 30 -X 50 -Y 60 -Z 70",

   the output data will contain the 3D data with

      x range = [start x : start x + size x] = [10:10+50] = [10: 60],
      y range = [start y : start y + size y] = [20:20+60] = [20: 80],
      z range = [start z : start z + size z] = [30:30+70] = [30:100].

   The default values are the entire simulation box.

   This tool can also be used to calculate the projected values along
   the x/y/z directions.

   It is a parallelized code, so the MPI libraries are required. User can
   specify the way of domain decomposition by the options

      "-p/q/r the number of MPI ranks in the x/y/z direction".

   Note that the way of domain decomposition specified here is not necessarily
   equal to that adopted in the parameter file of GAMER.

   The output format is as following:

      x   y   z   Density   Momentum x/y/z   Energy   Pressure


2. Command-line inputs:

   -b    output data in the binary form

   -h    display the synopsis

   -i    FILE
         name of the input file

   -I    INTERPOLATION_SCHEME
         0 : Min-Mod limiter
         1 : central difference

   -l    TARGETED_LEVEL
         targeted refinement level

   -n    OUTPUT_OPTION
         1 : output the x slice (yz plane)
         2 : output the y slice (xz plane)
         3 : output the z slice (xy plane)
         4 : output the x projection
         5 : output the y projection
         6 : output the z projection
         7 : output the 3-D data

   -o    SUFFIX
         suffix of the output file

   -O    load the old-format data

   -p    MPI_NRANK_X
         number of MPI ranks along the x direction for the domain decomposition

   -q    MPI_NRANK_Y
         number of MPI ranks along the y direction for the domain decomposition

   -r    MPI_NRANK_Z
         number of MPI ranks along the z direction for the domain decomposition

   -x    X_START
         starting x coordinate of the targeted region

   -y    Y_START
         starting y coordinate of the targeted region

   -z    Z_START
         starting z coordinate of the targeted region

   -X    X_SIZE
         size of the targeted region in the x direction

   -Y    Y_SIZE
         size of the targeted region in the y direction

   -Z    Z_SIZE
         size of the targeted region in the z direction

   -s    input grid scales instead of physical coordinates to specify the
         targeted range (default = off)

   -d    output the divergence of velocity

   -v    output the curl of velocity (vorticity)


3. Usage demo:

   (1) Output the xy plane with x=[100:500], y=[1000:1234], and z=5 at level 2.
       Use the MinMod limiter and 8 MPI processes (2 processes in each
       direction).

       mpirun -np 8 ./GAMER_ExtractUniform -i Input -l 2 -n 3 -x 100 -X 400
                                           -y 1000 -Y 234 -z 5 -p 2 -q 2 -r 2

   (2) Output the x projection with y=[0:BoxSize] and z=[0:BoxSize] at level 5.
       Use the central difference and 4 MPI processes in the z direction.
       Output the binary data.

       mpirun -np 4 ./GAMER_ExtractUniform -i Input -l 5 -n 4 -I 1 -r 4 -b

   (3) Output the 3D data with x=[x0:x0+dx], y=[y0:y0+dy], and z=[z0:z0+dz]
       at level 3. Use the MinMod limiter and a single MPI process. The
       suffix ".txt" is attached to the end of the output file name.

       mpirun -np 1 ./GAMER_ExtractUniform -i Input -o .txt -l 3 -n 7
                                           -x x0 -X dx -y y0 -Y dy -z z0 -Z dz


4. A simple MATLAB script "GAMER_PlotSlice.m" is put in the directory
   "Matlab_Plot". It can be used to plot the image of the output 2D data.



5. Caveat: For the fields "Potential (Pot)" and "ParticleDens (ParDens)", the
   extracted uniform-resolution results may be different when loading data
   from the C-binary and HDF5 snapshots. It is becaues HDF5 snapshots store
   the data of non-leaf patches as well, while C-binary snapshots only store
   the data of leaf patches and then use the restriction operation to obtain
   the data in the non-leaf patches.
