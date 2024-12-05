
   *********************************************************
   ** ELBDM Jean's instability test in the physical frame **
   *********************************************************

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

2. The "Input__Parameter" example file gives the UNSTABLE solution, in which
   the solution grows exponentially. To test the STABLE solution, one can
   set the gravitational constant NEWTON_G to a smaller value (e.g., 1.0e2).

