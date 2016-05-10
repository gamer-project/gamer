
   ************************************
   ** HYDRO Plummer hydrostatic test **
   ************************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories 

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Initial pressure is set to hydrostatic
   --> test whether the Plummer model is stable

2. Material will still propagate in from the boundary
   --> hydrostatics will not be able to hold for many free-fall time
       (it might be useful to use the analytical solution for the boundary
        ghost-zones )
