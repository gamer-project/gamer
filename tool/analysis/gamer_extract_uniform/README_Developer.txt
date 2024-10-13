
GAMER_ExtractUniform : construct the uniform-mesh data from a GAMER output file by spatial interpolation

==================================================================================================================


Version 1.7.2     04/09/2016
----------------------------
1. Add check the the boundary condition (only for HDF5 output)
   --> Terminate program if "-B 1" is not set for non-periodic simulation domain



1. Support the SERIAL mode
Version 1.7.1     03/14/2016
----------------------------
1. Support the SERIAL mode



Version 1.7.0     03/01/2016
----------------------------
1. Support file format >= 2000 (for both simple binary and HDF5)



Version 1.6.1     02/21/2015
----------------------------
1. Output the real physial coordinates even when "Shift2Center" is on



Version 1.6.0     07/15/2014
----------------------------
1. Support GAMER HDF5 output



Version 1.5.3     06/19/2014
----------------------------
1. Support converting both ELBDM and HYDRO peculiar velocity from comoving
   to physical coords (in km/s)
   --> Use the option "-k"
   --> For ELBDM, it must work with the option "-V"
       For HYDRO, it must work with either the option "-d" or "-v"



Version 1.5.2     06/12/2014
----------------------------
1. Rename "patch->ptr" as "amr->patch"



Version 1.5.1     05/15/2014
----------------------------
1. Support converting ELBDM peculiar velocity from comoving to physical coord.
   (in km/s)
   --> Use the option "-k"



Version 1.5.0     03/03/2014
----------------------------
1. Support loading the tree file (-T -e tree file)
   --> The tree file can be constructed by the tool "GAMER_ExtractTree.1.0.1"



1. Support different interpolation schemes (-w -I -m)
Version 1.4.2     12/06/2013
----------------------------
1. Support different interpolation schemes (-w -I -m)



Version 1.4.1     11/14/2013
----------------------------
1. Support periodic BC.
   --> Constraint: target size must be smaller than BoxSize - 2*BasePatchSize
2. New option -c --> Shift2Center, which can be used to
   (a) deal with the periodic BC.
   (b) make workload more balance
3. TakeNote now records the exact starting coord. and size instead of the
   input values. They can be different since the input values may not align
   with the underlying grids



Version 1.4.0     09/24/2013
----------------------------
1. Output velocity field in ELBDM
   --> can be plotted by Amira (see the introduction in this doc)



Version 1.3.1     10/10/2012
----------------------------
1. Support outputting the pressure with the option "-P"
   --> work for both text and binary outputs



Version 1.3.0     10/04/2012
----------------------------
1. Support non-periodic BC.
2. Record the command-line arguments



Version 1.2.2     12/25/2011
----------------------------
1. Properly declare "long"
2. Verify whether the input NGPU_X[0/1/2] are correct
3. Support OpenMP



Version 1.2.1     12/17/2011
----------------------------
1. Output both cell scales and physical coordinates
2. Record the range of physical coordinates instead of cell scales in the file
   name



Version 1.2.0     12/02/2011
----------------------------
1. Support new data format 1200
   --> Olde data format (1100 <= format < 1200) is still supported
2. Support both HYDRO and ELBDM models
3. New "MultiFiles__GAMER_ExtractUniform.sh" script



Version 1.1.0     09/20/2011
----------------------------

**********************************************************************
1. Fix the bug in "PreparePatch"
   --> In previous versions, buffer data obtained by interpolating on
       coarse-grid data are INCORRECT !!!!
**********************************************************************

2. Old data format is no longer supported
   --> the option "-O" is removed
3. Support outputting "potential, divergence(velocity), curl(velocity)"
   --> add the options "-d -v"
   --> potential will be outputted automatically if they are stored in the
       input file
4. Support using "physical coordinates" to specify the target
   --> It is the default setting
   --> To use "grid scale" instead, please turn on the option "-s"



Version 1.0.3     08/24/2010
----------------------------
1. Add a header line to the output file
2. A simple user manual "README_User.txt" is provided
3. The matlab script support working with the input files with headers



Version 1.0.2     08/05/2010
----------------------------
1. Work with the data format >= 1100 (version >= GAMER.1.0.beta3.2)
   --> reorder the fluid variables from "[PS][PS][PS][NCOMP]" to "[NCOMP][PS][PS][PS]"
2. Remove all TAB



Version 1.0.1     03/16/2010
----------------------------
1. Fix the bug in the function "LoadData"
   --> Declare the variable "Offset" as "long int" rather than "int"
   --> Work for larger input data



Version 1.0       02/01/2010
----------------------------
1. Inherit from GAMER_GetSlice and GAMER_GetUniformMesh
   --> can produce slice (-n 1/2/3) / projection (-n 4/5/6) /cube (-n 7) data

2. Specify the starting cooridnates --> command-line options "-x Start_X, -y Start_Y, -z Start_Z"
   Specify the targeted range       --> command-line options "-X Range_X, -Y Range_Y, -Z Range_Z"
   ** The coordinate system is consistent with the one adopted in GAMER **

3. Construct in parallel --> command-line options "-p NPROCESS_X, -q NPROCESS_Y, -r NPROCESS_Z"
   ** "-p, -q, -r" is NOT necessary to be the same values adopted in GAMER run **

4. Support binary output --> command-line option "-b"
   However, one must set "-p 1 -q 1 (slab decomposition)" for outputing the binary data

5. Work with the old data format (version < 1000) --> command-line option "-O"



Demo :
------
1. output x=[100,500], y=[1000,1234], z=5, level=2 plane (with MinMod limiter interpolation) using 8 processes
   mpirun -np 8 ./GAMER_ExtractUniform -i InputFile -o Comment -l 2 -n 3 -x 100 -X 400 -y 1000 -Y 234 -z 5 -p 2 -q 2 -r 2

2. output x projection, y=[0,BoxSize], z=[0,BoxSize], level=5 (with central interpolation), binary data
   using 4 processes
   mpirun -np 4 ./GAMER_ExtractUniform -i InputFile -o Comment -l 5 -n 4 -I 1 -r 4 -b

3. output cube, x=[x0,x0+dx], y=[y0,y0+dy], z=[z0,z0+dz], level=3, binary data using 1 process
   mpirun -np 1 ./GAMER_ExtractUniform -i InputFile -o Comment -l 3 -n 7 -b -x x0 -X dx -y y0 -Y dy -z z0 -Z dz



Binary outputs can be used to create image series by Amira in the following
steps:
---------------------------------------------------------------------------
1.  "View -> Background -> Adjust background color"
2.  "Load Time Series" --> select all binary files you want
3.  "TimeSeriesControl -> Display Time -> turn on [value text]"
4.  "Data -> BoundingBox -> make [Line width] thicker and probably
    change color"
5.  (OPTIONAL) "Data -> Compute -> Arithmetic -> ln(a)
    -> press the [auto-refresh] button"
6.  (3D)
    "Result -> Display -> Voltex -> set the colormap
    -> press the [auto-refresh] button"
    (2D)
    "Result -> Display -> OrthoSlice -> Mapping Type=Colarmap
    -> set the Colormap (can use the "Adjust range" option to set the range
       automatically)
    -> press the [auto-refresh] button"
    -> properly set the viewing angle (View YZ/XY/XZ)
7.  (OPTIONAL) Show colormap in the viewer window as legend by
    "Pool -> Show Objects -> choose the colarmap item
     -> right-clicking the colormap item -> DisplayColormap"
8.  Rotate the bounding box to a proper orientation
9.  "TimeSeriesControl -> set [Cached steps] = 0"
10. (for creating a single AVI movie)
    "TimeSeriesControl -> MovieMaker -> set [File format = AVI movie]
    -> set [Filename = Movie.avi] -> set [Frames] and [Frame rate]
    -> set [AVI encoder = ffdshow video encoder]
    -> set [Tiles X=2, Y=2] (if you want the image resolution to be higher
       than the screen)
    -> set [Resolution X=1280, Y=800]
    -> set [AntiAlias 1]
    -> Apply
    -> Ensure that the active window is not overlapped by anything during
       making movie"

    (for creating multiple PNG images)
    "TimeSeriesControl -> MovieMaker -> set [File format = PNG images]
    -> set [Filename = XXX]
    -> set [Frames = total number of data dumps]
    -> set [Tiles X=2, Y=2] (if you want the image resolution to be higher
       than the screen)
    -> set [Resolution X=1280, Y=800]
    -> set [AntiAlias 1]
    -> Apply
    -> Ensure that the active window is not overlapped by anything during
       making movie"


Vector plot in Amira:
---------------------------------------------------------------------------
1. Open Vx, Vy and Vz separately
2. Vx field -> Compute -> ChannelWorks
3. ChannelWorks -> Add Vy and Vz as Input2 and Input3, respectively (using
   right mouse button)
4. ChannelWorks -> Output = 3 channels, Channel 1/2/3 = Input 1/2/3 channel 1
5. Newly created object -> Display -> Vectors
6. The default plot is the "3D" vector. To plot the "2D" vector (e.g., ignoring the
   z component in a xy slice), one should (1) disable "perspective" view (the
   symbol with an "eye") and (2) enable "projection"
7. To have equal-size arrows, turn on "constant"



Add new output elements :
-------------------------
To add new output elements, one can follow what is done for the option
"OutputDivVel".



Unfinished work :
-----------------
1. Paste the sibling data after the first interpolation
2. Periodic boundary condition
3. Add the check for the size of the type "long int"
4. Support cubic decomposition for binary output
5. Support more accurate interpolation schemes
7. Set ExtBC according to the loaded data


