
   ******************************
   ** Two particles orbit test **
   ******************************

Procedure to run the test problem:
-------------------------------------------------------------------------------
1. Execute the script "Copy_ExtPot/Acc.sh" to copy files to their corresponding
   directories
   --> ExtPot: test external potential
       ExtAcc: test external acceleration

2. Compile GAMER by typing "make" in the directory "GAMER/src"
   ("make clean" first if the program has been compiled previously)

3. Run GAMER by executing the file "Dizzy" in the directory "GAMER/bin/Run"


Note:
-------------------------------------------------------------------------------
1. This test tests the binary orbit of two particles

2. Support both self-gravity, external potential and external acceleration
   --> External acceleration is supported only in MODEL=PAR_ONLY/HYDRO
       External potential    is supported only in MODEL=ELBDM

3. When using external potential, particles mass will be reset to extremely
   small values so that self-gravity can be ignore
   --> Because currently the ELBDM model must work with self-gravity
