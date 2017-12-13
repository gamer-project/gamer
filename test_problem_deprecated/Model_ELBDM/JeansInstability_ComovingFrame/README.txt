
   *********************************************************
   ** ELBDM Jean's instability test in the comoving frame **
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

2. DT__FLUID is set to 0.125, and DT__MAX_DELTA_A is set an extremely large
   number so that the time-step will be controlled by the time-step criterion
   of the fluid solver

3. Analytical solution reference: Woo, T., & Chiueh, T. 2009, ApJ, 697, 850

** Note that the imaginary part (I) grows much faster than the real part 
   (R), and the accuracy of the linear prediction will be deteriorated 
   when the assumption (2R >> I^2) is slightly broken down, especially in the
   higher resolution tests (and hence the 2nd-order accuracy may no longer
   hold).

4. The "Input__Parameter" example file gives the UNSTABLE solution, in which
   the solution grows with time. To test the STABLE solution, one can
   try the following parameters:

      BOX_SIZE  = 0.5
      A_INIT    = 1.0e-5
      OUTPUT_DT = 2.0e-6

   Besides, one should also set 

      Jeans_RealAmp0 = 4.e-8;

   in the file "Init_TestProb.cpp" to ensure that the double precision is
   enough in the 256^3 STABLE run.
