
   ***************************
   ** ELBDM plane wave test **
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
1. No grid refinement is adopted in this test
2. This test will overwrite the files "Init_GAMER.cpp", "Flu_Close.cpp" and "Patch.h".
   Please backup these files in advance.
