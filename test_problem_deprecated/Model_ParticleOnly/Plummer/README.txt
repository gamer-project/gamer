
   ************************
   ** Plummer model test **
   ************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories 

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Test the evolution of the Plummer model

2. The gnuplot script 'plot.gpt' can be used to compare different results

3. This test supports two modes
   (1) MODEL == ELBDM: grid density is set to zero. So essentially there are particles only
   (2) MODEL == HYDRO: gas and particles share the same density profile, and the mass ratio is
                       determined by the input parameter "Plummer_GasMFrac"
