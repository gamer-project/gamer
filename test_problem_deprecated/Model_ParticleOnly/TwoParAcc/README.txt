
   ******************************
   ** Two particles force test **
   ******************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. This test randomly places two particles and calculate the gravitational
acceleration and errors

2. This test rewrite the initialization procedure to improve the performance
   --> program is terminated in the function "Init_TestProb"
