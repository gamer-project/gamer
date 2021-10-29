

******************************************************************************
GAMER_Data2Visit : transform 2D slice or 3D volume of GAMER binary outputs
                   to the silo format, which can be plotted by "VisIt"
******************************************************************************


NOTE :         12/01/2011
-------------------------
1. Work on both 2D slice and 3D volume
   --> User should provide the corner coordinates and the size of the targeted
       surface/volume
   --> corner coordinates  : -(x,y,z)
       volume size         : -(X,Y,Z)
       target surface      : -m (0,1,2,3) <--> (xy,yz,xz,3D)

   Example :

      ./GAMER_Data2Visit InputFilename -m 0 -x 0.0 -y 1.0 -z 0.3 \
                                       -X 1.0 -Y 0.5 -Z 1.5

   --> plot the z=0.3 plane, with xrange=[0.0 ... 1.0], yrange=[1.0 ... 1.5]
   --> the option "-Z 1.5" is useless when plotting a xy plane and hence is
       ignored in this case

2. VisIt working procedure:
   a. load file
      a1. open the file "Root.silo"
      a2. open the time-varying database "Group.visit" (movies can be
          generated in this mode)
      a3. File / Restore session :
          --> restore a previous image
      a4. File / Restore session with sources :
          --> use previous settings to generate new images

   b. Controls / Expressions
      --> generate new variables from the existing variables (e.g.,
          calculate vorticity vector from the input velocity vector)
      --> These expressions will appear as new scalar/vectors for plotting
      Example :
      b1. divergence(Vel)
      b2. curl(Vel)

   c. Plots
      c1. Plots / Add->Mesh
          --> plot the domain refinement result
          (1) Box   : outline of the targeted surface
          (2) Mesh  : outline of all grid cells
          (3) Patch : outline of all patches
      c2. Plots / Add->Pseudocolor
          --> plot the density, energy, and velocity magnitude
      c3. Plots / Add->Vector
          --> plot the velocity vectors
      c4. Plots / Add->Volume
          --> 3D volume rendering

   d. Apply operators
      Example :
      d1. Plots / Operators->Transforms->Transform
           --> apply rotate/scale/translate to the original data
      d2. Plots / Operators->Slicing->Slice
          --> cut a single slice from the 3D data
      % After adding operators, one can double click the operator items to
        edit their parameters

   e. Begin to draw : Plots / Draw

   f. Save results
      f1. File / Set Save options
          --> save a single image
      f2. File / Save movie
          --> save a series of images (to do this, one must load the
              time-varying file "Group.visit")
      f3. File / Save session as
          --> save all current settings to a sessions file, which can be used
              to restore the image latter
      f4. Options / Save Settings
          --> save all current "default" settings

3. Only patches within the targeted region will be allocated
   --> considerably save the total memory consumption

4. Size of the targeted region will always be truncated to fit the size of
   the root-level patch
   --> in other words, the root-level patch is the basic unit for specifying
       the targeted region (for both 2D and 3D plots)

5. The script "MultiFiles__GAMER_Data2Visit.sh" can be used to process
   multiple files

6. How to generate MULTIPLE figures automatically?
   a. Load the "Group.visit" file

   b. Properly set all settings for figure 1. For example,
      b1. enable log plot
      b2. set the color bar range
      b3. set the annotation properly
          b3.1: disable "User information"
          b3.2: 3D->X-Axis: disable "Title", set "Font scale" of
                "Labels" to 2
          b3.3: 3D->Y-Axis: disable "Title", set "Font scale" of
                "Labels" to 2
          b3.4: 3D->Z-Axis: disable "Title"
          b3.5: 3D->General 3D: set "Tick mark location" to "Outside",
                set "Axis type" to "Static edges", disable "Auto set
                tickes" and enable it again
      b4. draw figure 1 to ensure that all settings are ok

   c. File / Save movie
      c1. Format: png
      c2. Movie size: 1280:800
      c3. Press the bottom "->"
      c4. Set the file name
      --> Multiple PNG figures will be created after this step

   d. Use other softwares (e.g., the "convert" command) to convert all PNG
      figures to a movie file
      --> An example script "CreateMovie.sh" is put in the Run directory




Version 2.4.0  03/02/2016
-------------------------
1. Support file format >= 2000 (for both simple binary and HDF5)


Version 2.3.0  07/14/2014
-------------------------
1. Support GAMER HDF5 output
2. Update HDF5 version from 1.8.8 to 1.8.9


Version 2.2.3  06/12/2014
-------------------------
1. Rename "patch->ptr" as "amr->patch"


Version 2.2.2  08/20/2013
-------------------------
1. Output the gravitational potential


Version 2.2.1  04/09/2013
-------------------------
1. Output the effective total velocity (Ve) in ELBDM, which is defined as

      Ve = SQRT( Vt^2 + V^2 ) = SQRT( 2*K.E/DENS ) = 1.0/eta*( |grad(psi)|/|psi| ),

   where Vt is the thermal velocity and V is the bulk velocity


Version 2.2    02/13/2012
-------------------------
1. Output bulk velocity vector V and thermal velocity magnitude VT
   (internal energy = 0.5*rho*VT^2) in ELBDM


Version 2.1    12/01/2011
-------------------------
1. Support new data format 1200
   --> Olde data format (1100 <= format < 1200) is still supported
2. Support both HYDRO and ELBDM models
3. Support using both "cell scale" and "physical coordinates" to specify the
   targeted range
   --> Using "physical coordinates" is the default setting
   --> To use "cell scale" instead, please turn on the option "-s"
4. New "MultiFiles__GAMER_Data2Visit.sh" script



Version 2.0    04/11/2011
-------------------------
1. Support 3D volume rendering
   --> the option "-m/n" is replaced by "-X/Y/Z"
   --> the option "-s" is replaced by "-m"
2. Support pressure output
   --> use the option "-p"
3. Support vorticity vector output
   --> use the option "-w"
4. Record "Cycle (DumpID) and Time" in the silo files
5. A text file named "Group.visit" will be generated by the script
   "MultiFiles__GAMER_DataToVisit.sh" automatically. This file can be
   used to define a time-varying database in VisIt.



Version 1.1    02/09/2011
-------------------------
1. Work with "GAMER.1.0.beta4.1.2"
2. Integrate data at all levels
   --> No longer need to fix the color range manually
   --> Data can still be separated into different levels by the option "-d"
3. Use "physical coordinates" rather than "grid indices" to specify the
   targeted surface
4. New "MultiFiles__GAMER_DataToVisit_Slice.sh" script
   --> Use input filename as the first argument
   --> Input filename will be used to create the new directory
       (words in front of the last slash will be removed)
   --> For example:

       sh MultiFiles__GAMER_DataToVisit_Slice.sh ../../Data 0 2

       "../../Data_000000, ../../Data_000001, and ../../Data_000002" will be
       processed, and the results will be stored in the directory
       "Data_000000, Data_000001, and Data_000002" respectively



Version 1.0    07/09/2009
-------------------------
1. First version



Bug :
-------------------------
1. Vector plot seems to have a little problem
   --> occationally, different vectors will be connected by lines



Unfinished :   04/11/2011
-------------------------
1. Periodic boundary condition
2. Set the targeted maximum level "TargetLv"
   --> perform restrict operatin if TargetLv < NLEVEL-1
