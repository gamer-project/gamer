
   ***************************
   ** HYDRO blast wave test **
   ***************************

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
1. Lohner's error estimator for both mass and energy density are adopted as the refinement criteria
2. By default, the maximum refinement level = 3
3. By default, the refinement threshold of the Lohner error estimator = 0.80
