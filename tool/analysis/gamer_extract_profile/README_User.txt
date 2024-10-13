
GAMER_ExtractProfile : Mode A: evaluate the shell averages of fluid variables
                       Mode B: get the maximum density

==============================================================================

Version 1.0.3     08/23/2010
----------------------------
1. Two modes are supported in this tool.

   (1) Mode A "shell average":

         Evaluate the shell average and RMS of the targeted sphere. User
         can provide the center, radius, and number of shells for the
         targeted sphere.

         If the option "-m" is turned on, the user-provided sphere center
         "-x/y/z" will be used as the first guess of the sphere center. Then,
         the program will search for the eight most massive grids, and use
         their average position as the final sphere center.

         Five output files "AveRho, AveV_R, AveV_T, AvePre, AveEgy" will be
         created, which record the shell averages of density, radial
         velocity, tangential velocity, pressure, and energy, respectively.
         The output format is as following:

            Radius   NCount   AveXXX   RMSXXX   MaxXXX   MinXXX

            Radius : the radius of the shell
            NCount : the number of cells located within this shell
            AveXXX : the shell average of the variable XXX
            RMSXXX : the RMS of the variable XXX
            MaxXXX : the maximum value of the variable XXX in this shell
            MinXXX : the minimum value of the variable XXX in this shell

         This mode can be turned on by the option "-S".


   (2) Mode B "maximum density":

         Find the grids with the mass density exceeding the user-provided
         threshold (-t THRESHOLD).

         This mode can be turned on by the option "-M".


2. Command-line inputs:

   -h    display the synopsis

   -i    FILE
         name of the input file

   -m    use the position with the maximum density as the sphere center

   -M    turn on the mode B "maximum density"

   -n    NSHELL
         divide the targeted sphere into NSHELL shells.

   -o    SUFFIX
         suffix of the output file

   -O    load the old-format data

   -p    use the periodic boundary condition

   -r    RADIUS
         radius of the targeted sphere

   -S    turn on the mode A "shell average"

   -t    THRESHOLD
         density threshold for the mode B "maximum density"

   -x    X
         x coordinate of the center of the targeted sphere

   -y    Y
         y coordinate of the center of the targeted sphere

   -z    Z
         z coordinate of the center of the targeted sphere


3. Usage demo:

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



