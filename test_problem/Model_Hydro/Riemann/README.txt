
   ********************************
   ** HYDRO Riemann problem test ** 
   ********************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy.sh" to copy files to their corresponding
   directories 

      Init_TestProb.cpp      --> GAMER/src/Init
      Makefile               --> GAMER/src
      Input__*               --> GAMER/bin/Run

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. Test problem parameters cab be set in the file "Input__TestProb"
2. By default, the maximum refinement level = 2
3. By default, the refinement threshold of the Lohner error estimator = 0.50
