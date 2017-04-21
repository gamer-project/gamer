
   *********************************************************
   ** HYDRO Jean's instability test in the physical frame **
   *********************************************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Only support the 3D case

2. Example parameters:
   Jeans_Cs = 1.0 -->   stable solution
   Jeans_Cs = 0.1 --> unstable solution

