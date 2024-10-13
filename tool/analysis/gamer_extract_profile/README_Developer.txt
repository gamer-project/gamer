
GAMER_ExtractProfile : Mode A: evaluate the shell averages of fluid variables
                       Mode B: get the maximum density

==================================================================================================================


Version 1.7.0     03/01/2016
----------------------------
1. Support file format >= 2000 (both simple binary and HDF5)



Version 1.6.1     12/17/2015
----------------------------
1. Support calculating the spherically symmetric gravitational potential
   --> -g -G
   --> also calculate potential energy when ELBDM_GetVir is on



Version 1.6.0     07/14/2014
----------------------------
1. Support GAMER HDF5 output



Version 1.5.1     06/12/2014
----------------------------
1. Rename "patch->ptr" as "amr->patch"



Version 1.5.0     02/28/2014
----------------------------
1. Support loading the tree file (-T -e tree file)
   --> The tree file can be constructed by the tool "GAMER_ExtractTree.1.0.1"
2. Tree structure will NOT be constructed if no ghost zones are required
   --> Improve I/O performance



Version 1.4.2     12/14/2013
----------------------------
1. Record Rho*(vr^2+wr^2) in "VirSurf", which does NOT equal to the sum of
   Rho*vr^2 + Rho*wr^2 because of different finite difference schemes
   --> avoid large error in wr when rho ~ 0
2. Record the sum of all surfaces terms in the virial condition, which should
   be equal to the negative gravitational potential energy. Both the results
   with Ek-L and Ek-G are included (Sum-L & Sum-G).



Version 1.4.1     12/12/2013
----------------------------
1. Support analyzing the surface terms in the ELBDM virial condition (-V)
   --> Create the file "VirSurf"
2. Remove the -K option, which will be enabled automatically with -V



Version 1.4.0     11/10/2013
----------------------------
1. Support preparing ghost zones for each patch
   --> Derived variables requiring stencils can be evaluated (e.g., kinematic
       energy in ELBDM)
2. Different interplation schemes are supported (-T -P -u)
3. Supporting evaluting the shell-average and accumulated kinematic energy in
   ELBDM (-K)
   --> Center-of-mass kinematic energy can be subtracted by -V



Version 1.3.0     10/07/2013
----------------------------
1. Support defining a maximum radius for the option -m to determine the sphere
   center
   --> By the command-line options "-R maximum radius for the option -m"
   --> Useful when there are two massive objects approaching each other and
       one would like to stick to the center of one of the objects even though
       it's peak density may be lower than the other one

2. Support outputting the center coords, which can be useful for tracing the
   trajectory of the target object and setting the sphere center
   automatically
   --> An example for setting the center coords. automatically has been put in
       the bash script MultiFiles__GAMER_ExtractProfile.sh
   --> Use the option "-c" to output the center coords



Version 1.2.3     05/28/2013
----------------------------
1. Support computing shell average at log bins
   --> By the command-line options "-L LogBinSize"



Version 1.2.2     04/18/2013
----------------------------
1. Output the average potential



Version 1.2.1     06/11/2012
----------------------------
1. Output the accumulated average density



Version 1.2.0     02/10/2012
----------------------------
1. Output the accumulated mass



Version 1.1.3     12/09/2011
----------------------------
1. The argument "-a" now takes an input parameter N
   --> "-a N"
   --> Set the shell width to N*(the size of cell closest to the targeted
       sphere center)
   --> Default value = 2
   --> Disable this option by inputting a negative number



Version 1.1.2     12/08/2011
----------------------------
1. Number of highest-density cells for determing the sphere center is
   adjustable
   --> use the option "-m ###" to input the number
   --> set -m to any negative number to disable this option



Version 1.1.1     12/07/2011
----------------------------
1. Support to determine the shell width automatically by searching the size of
   cell closest to the targeted sphere center
   --> enable it by the option "-a"

2. The "SetMaxRhoPos" function now uses the "center of mass" to get the sphere
   center



Version 1.1.0     12/02/2011
----------------------------
1. Support new data format 1200
   --> Olde data format (1100 <= format < 1200) is still supported
   --> format < 1100 is no longer supported, and the option "-O" is removed
2. Support both HYDRO and ELBDM models
3. Support using both "cell scale" and "physical coordinates" to specify the
   targeted range
   --> Using "physical coordinates" is the default setting
   --> To use "cell scale" instead, please turn on the option "-s"
4. New "MultiFiles__GAMER_GetCube.sh" script



Version 1.0.3     08/23/2010
----------------------------
1. Use different weightings (proportional to the cell volume) for cells in different refinement levels to
   evaluate the RMS.
2. A simple user manual "README_User.txt" is provided.



Version 1.0.2     08/05/2010
----------------------------
1. Work with the data format >= 1100 (version >= GAMER.1.0.beta3.2)
   --> reorder the fluid variables from "[PS][PS][PS][NCOMP]" to "[NCOMP][PS][PS][PS]"
2. Remove all TAB



Version 1.0.1     02/02/2010
----------------------------
1. Use different weightings (proportional to the cell volume) for cells in different refinement levels to
   evaluate the shell average.
2. The average velocity is re-defined as "Ave(momentum)/Ave(density)"



Version 1.0       01/22/2010
----------------------------
1. This new tool combines the two old tools "GAMMA_ShellAverage" and "GAMMA_GetMaxDensity".

   --> Use the command-line option "-S,-M" to determine the usage:
         -S : the "shell-average"  mode
         -M : the "maximum-densty" mode

   --> Two modes can be turned on simultaneously

2. Work with both the old (GAMER.1.0.beta1) and the new (GAMER.1.0.beta3.1 --> format version >= 1100) data format
   --> Use the command-line option "-O" to load the old-format data.



Demo :
------

   (1) Mode A "shell average":

      (a) Specify the sphere center, radius, and number of shells:

          ./GAMER_ExtractProfile -i Input -S -r 128 -n 128 -x 1 -y 2 -z 3 -p

      (b) Give the first guess of the sphere center, and then determine the
          sphere center by finding the maximum density:

          ./GAMER_ExtractProfile -i Input -S -r 128 -n 128 -x 1 -y 2 -z 3 -p -m

      (c) Determine everything automatically:

          ./GAMER_ExtractProfile -i Input -S -p -m


   (2) Mode B "maximum density":
      (1-a) - (1-c) can also be applied for this mode. For example:

       ./GAMER_ExtractProfile -i Input -M -p -m -t 100



Unfinished work :
-----------------
1. Remove the C.M. velocity for evaluating the shell-average of velocity.
2. Add the check for the size of the type "long int"
3. Add the warning message when the format version is not supported
   --> Either an old file or a unrecognizable file
4. OpenMP



