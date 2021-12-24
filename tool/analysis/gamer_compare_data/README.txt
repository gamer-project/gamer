
GAMER_DataCompare : compare two binary outputs of GAMER

==================================================================================================================


Version 2.3    03/01/2016
-------------------------
1. Support new data format 2000. Stop supporting format < 2000



Version 2.2    12/03/2011
-------------------------
1. Support new data format 1200
   --> Olde data format (1100 <= format < 1200) is still supported
   --> format < 1100 is no longer supported, and the option "-O" is removed
2. Support both HYDRO and ELBDM models



Version 2.1    12/01/2010
-------------------------
1. Fix the bug when comparing the potential data
   --> The keyword "GRAVITY" is no longer necessary
   --> Only rely on the input parameter "opt__output_pot"
2. Support "xyzv" and "vxyz" data formats
3. Add the "isfinite" check



Version 2.0    01/21/2010
-------------------------
1. Suuport the new data format (version 1100) for the program GAMER.1.0.beta3.1.
2. Add the command line option "-F" to determine the data format of the input files.



Version 1.2    10/03/2009
-------------------------
1. add the option "-c"
   --> compare patches with the same corner coordinates
2. verify that all patches without son have been checked
   --> add the data member "check"



Version 1.1    09/07/2009
-------------------------
1. only compare the data stored in the patches without son
2. only compare the potential data if both two input files have potential data



Version 1.0    04/22/2009
-------------------------
1. output the absolute and relative error between two input files



Demo :
------
./GAMER_DataCompare -i InputFile1 -j InputFile2 -o OutputFile -e TolerantError

./GAMER_DataCompare -i InputFile1 -j InputFile2 -o OutputFile -e TolerantError -c

./GAMER_DataCompare -i InputFile1 -j InputFile2 -o OutputFile -e TolerantError -c -F 0



Unfinished work :
-----------------

