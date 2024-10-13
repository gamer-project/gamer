

                        =================================================================
                        ||                                                             ||
                        ||    GAMER : GPU-accelerated Adaptive-MEsh-Refinement code    ||
                        ||                                                             ||
                        =================================================================

------------------------------------------------------------------------------------------------------------------

sibling index (SibID):

 0 (-,0,0)
 1 (+,0,0)
 2 (0,-,0)
 3 (0,+,0)
 4 (0,0,-)
 5 (0,0,+)
 6 (-,-,0)
 7 (+,-,0)
 8 (-,+,0)
 9 (+,+,0)
10 (0,-,-)
11 (0,+,-)
12 (0,-,+)
13 (0,+,+)
14 (-,0,-)
15 (-,0,+)
16 (+,0,-)
17 (+,0,+)
18 (-,-,-)
19 (+,-,-)
20 (-,+,-)
21 (+,+,-)
22 (-,-,+)
23 (+,-,+)
24 (-,+,+)
25 (+,+,+)


local patch index (LocalID) in one patch group :

0  (0,0,0)
1  (1,0,0)
2  (0,1,0)
3  (0,0,1)
4  (1,1,0)
5  (0,1,1)
6  (1,0,1)
7  (1,1,1)



Version 1.0     8/4/2008
------------------------
1. rename from the AMR.3.1 program

2. written by Yu-Chih Tsai

3. some bugs will be fixed in the version 1.1



Version 1.1     8/13/2008
-------------------------
1. add the UM_Init.cpp file for initialization from a density distribution stored in a uniform mesh

2. add the Auxiliary.cpp file for some auxiliary functions

3. correct all warning messages



Version 1.2     8/29/2008
-------------------------
1. add several options in the Makefile

2. move global variables into the Global.h file

3. redefine the Time_constant in the Define.h
   --> two times smaller than the origianl definition

4. modify the program structure to fit my style

5. output the data of a single patch
   --> Output_Patch

6. record the level occupancy percentage
   --> Aux_PatchCount

7. output the flux of a single patch stored for fix-up operation
   --> Output_FixUpFlux

8. output the data of a single patch plus its ghost zone prepared by the "Prepare_solve" function
   --> Output_PreparedPatch

9. output the flux at one surface of a single patch stored in the S_array after the "Prepare_solve" call
   --> Output_PatchFlux

10. verify all parameters at begining
   --> Aux_Check_Parameter

11. put the refinement criteria for different levels (for density only) in the "Zodiacal Constellation" file

12. put all input file in the "Input" directory

13. put all header files in the "Header" directory

14. include the "GAMMA.h" file will automatically include all header files

15. add the ABSOLUTE operation for the sound speed evaluation in the "Determine_dt" function

16. check the conservation law
   --> Aux_Check_Conservation

17. trigger some auxiliary functions by calling the "Auxiliary" function

18. check the negative density and negative pressure
   --> Aux_Check_Negative

19. move prototype declarations into the Prototype.h file

20. move the BLOCK_SIZE_Y macro to the Define.h file

21. rearragne all options in the Makefile

22. the result has been verified by comparing to GAMMA.1.1



Version 1.3     9/2/2008
------------------------
1. the basic timing functions is added

2. fix the bug in the "Determine_dt" function
   --> take absolute value for velocity
   --> do NOT modify the original value

3. fix the bug in the timing inside the "Integration" funtion

4. record the present time for each level
   --> store in the "Time" array

5. the Init_Reload has been verified
   --> it should give EXACTLY the same result as if the program has never been terminated

6. the output data now contains the "Time" of each level and also record the output "Step"

7. the result has been verified by comparing to GAMMA.1.2
   --> it will no longer give exactly the same result as GAMMA.1.2 due to the following changes in GAMMA.1.3

   a. Time_Step   : float -> double
   b. Time        : float -> doulbe
   c. output "Time" for each level
   d. output "Step"
   e. add ABSOLUTE operator in the Determine_dt function



Version 2.0     9/9/2008
------------------------
1. the Min_level and Max_level macros have been replaced by the NLEVEL macro
   --> now the level ranges from 0 to NLEVEL-1

2. the Prepare_solver function no longer needs to find the siblings of father for the base level
   --> by using the Table_26
   --> the program no longer needs the Min_level-1 level ( or the -1 level )

3. the physical domain is no longer necessary to be cubic

                        *********************************
                        **                             **
   --> it can deal with ** rectangular physical domain **
                        **                             **
                        *********************************

4. all macros are highlighted by using capital letters

5. the result of the case that NX0==NY0==NZ0 has been verified by comparing to GAMMA.1.3
   --> it will no longer give exactly the same result as GAMMA.1.3 due to the following changes in GAMMA.2.0

   a. the Init_AssignData function has been modified, the binary output will not be exactly the same
   b. the order of patch IDs in all levels are different (different sibling, father, and son)

6. the UM_INIT option has been verified by comparing to GAMMA.1.3

7. the RELOAD option has been certified (for NX0==NY0==NZ0)

8. the Flag_Init function has been removed
   --> the "Flag" function now works both during initialization and hereafter

9. the Refine_Init.cpp file has been removed
   --> the Refine_Init function has been moved into the Init.cpp file and renamed as Init_Refine
   --> the Sibling_Search function has been moved into the Miscellaneous.cpp file

10. correct the bug in the "pdelete" function
   --> do NOT set the fathers of patches in the level "level+1" as -1 ( it will result in error if level+1=NLEVEL)

11. the result of the case that NX0!=NY0!=NZ0 has been verified by comparing to the mhtvd3_3 program

12. load some parameters from the "Parameter" file, move other parameters into the "Makefile"
   --> the program can be triggered with different parameters without re-compiling

13. move some options from the "Makefile" into the "Parameter" file

14. for the case that NX0!=NY0!=NZ0, the OPT__RELOAD option will NOT give exactly the same result !!

15. record NX0, NY0, and NZ0 in the output file
   --> when doing Reload, it will check whether the N(X,Y,Z)0 recorded in the RESTART file is equal to the
       N(X,Y,Z)0 loaded from the "Parameter" file



Version 3.0     10/2/2008
-------------------------

************************
**                    **
**  PARALLEL VERSION  **
**                    **
************************

1. the result has been verified by comparing to GAMMA.2.0   (differenct by only round-off error)

2. the Reload function is complete and verified

3. currently the UM_Init function does NOT work



Version 3.1     10/6/2008
-------------------------
1. double check the total elapsing time of simulation by recording the CPU clock

2. reduce the useless work in the Interpolation function

3. rename some macros :

   PATCH_NX, PATCH_NY, PATCH_NZ     ->    PATCH_SIZE
   GHOST_WIDTH                      ->    GHOST_SIZE
   BUFFER                           ->    PARA_BUFFER_SIZE



Version 3.1.1   10/10/2008
--------------------------
1. simplify the function "GetExchangeFluxPatchID" by only scanning over the BoundP_IDMap in the level "lv"



Version 3.1.2   10/20/2008
--------------------------
1. remove the line

      if ( patch.num[lv] == 0 )  break;

   in the loop inside the function "Init_StartOver"

   --> it's a bug since that some initialization procedures will still have to be performed even though there is
       no patch in the current level
       e.g. initialize the BounP_IDMap as -1


2. correct the bug when " NX0_TOT[0,1,2] / NGPU_(X,Y,Z) == 16 "
   ~ by modify the function "Loop_width_5" and its usages in all calling functions



Version 3.1.3     11/6/2008
---------------------------
1. declare the variables "HeaderSize, DataSize[NGPU], offset" in the function Init_Reload as "long int"
   --> it's a bug in the version 3.1.2 which will cause the program failed if the RESTART file is too large



Version 4.0       3/10/2009
---------------------------

********************
**                **
**  SELF GRAVITY  **
**                **
********************

1. with self gravity --> too many changes ...

2. rearrange the order of directions : x,y,z -> z,y,x  (consistent with the solver)
   --> better performance

3. to make the result of re-start be exactly eqaul to the original result, three lines must be added into the
   program
   (1) AveDensity = 1.0;
      --> put it in the end of the function "Get_AverageDensity"
      --> otherwise the average density may be slightly different
   (2) Restrict( lv, 0, patch.evolve_counter[lv]%2, 0, 4 );
   (3) GetBufferData( lv, patch.evolve_counter[lv]%2, 1, 3 );
      --> put (2) and (3) in the end of the function "Refine"
      --> otherwise the coarse-patch value may NOT be exactly equal to the average of fine-patch value



Version 4.1       3/23/2009
---------------------------

**********************
**                  **
**  COMOVING FRAME  **
**                  **
**********************

1. fix the bug in the function "G_GetPot_AllLevel"
   --> we still need to transfer "1+(POT_GHOST_SIZE+1)/2" data even if "NPatchComma[lv+1][27] == 0", because
       the sibling ranks may have "NPatchComma[lv+1][27] != 0"

2. complete the parallel version of the function "UM_Init"

3. add the option "OPT__TUTELAR_DP"
   --> the file "Zodiacal_Constellation" can be used as the refinement criteria even if "OPT__INIT != 2"

4. remove all unnecessary lines in the function "Flag"

5. fix one bug in the functions "Integration_IndividualTimeStep" and "Integration_SharedTimeStep"
   --> add one more criterion "if ( NLEVEL == 1 )" to determine whether to call the function "GetBufferData"
       after the funtion "FixUp"
   --> otherwise the program will NOT get the correct buffer data if NLEVEL == 1

6. add the option "OPT__FIXUP_FLUX" to suspend the flux fix up operation
   --> since the density flux fix up operation will introduce the instability ...

7. add the option "OPT__DUMPTABLE" to load the dump table
   --> it will output data ONLY according to the dump table
   --> it will NOT output the initial and final data if they are not the output time listed in the dump table

8. add the parameter "FLAG_BUFFER_SIZE"
   --> F : the flagged grid due to over-density
       B : the flagged grid due to FLAG_BUFFER_SIZE

       e.g. (a). FLAG_BUFFER_SIZE == 1

                 B B B
                 B F B
                 B B B

            (b). FLAG_BUFFER_SIZE == 2

                 B B B B B
                 B B B B B
                 B B F B B
                 B B B B B
                 B B B B B



Version 4.2       3/30/2009
---------------------------
1. add two more criteria for the time-step determination
   a. DT__FREE_FALL
   b. DT__MAX_DELTA_A

2. remove the parameter "MAX_DELTA_A"

3. add the fucntion "Output_Part" for outputing slices or lines



Version 4.2.1     4/1/2009
--------------------------
1. the parameter "OPT__MONO_INTERPOLATION" is renamed as "OPT__INTERPO_SCHEME"
   --> add two more interpolation schemes
       a. MinMod limiter  --> use the minimum slope between the left and right slopes
       b. vanLeer limiter --> use the harmonic mean of the left and right slopes

   --> both these two interpolation schemes are monotonic and conservative

2. the density flux fixup operation seems to introduce instability only when working with the
   linear interpolation scheme (since this interpolation scheme is NOT conservative)



Version 4.2.2     4/6/2009
--------------------------
1. deallocate flux array in the function "Refine"
   --> otherwise the program may have memory leakage !!!

2. work with CUDA.2.1



Version 4.2.3     4/9/2009
--------------------------
1. even when the option "COMOVING" is turned on, the function "UM_Init" will NO LONGER add the
   background density (ONE)
   --> it is assumed that the file "UM_START" stores DENSITY rather than DELTA



Version 4.2.4     4/10/2009
---------------------------
1. disable the potential correction when NLEVEL == 1
   --> by resetting
       G_NCORR             = 0;
       G_NFINAL_SIB_RELAX  = 0;
       OPT__SIB_RELAX      = false;

2. add warning message for the vanLeer-limiter interpolation
   --> it may still introduce negative density !!
   --> currently the MinMod limiter is the most safe choise

3. set the sound speed automatically in the function "UM_INIT" for the cosmological simulation

4. move all symbolic constants from the Makefile to the file "Symbolic_Constant.h"

5. move most of the c++ headers into the file "GAMMA.h"

6. include more c++ headers
   --> now it can be compiled successfully in the NCHC GPU cluster (still need to remove some optimization flags)

7. slightly modify the output format of the function "Aux_PatchCount"



Version 4.3       4/13/2009
---------------------------
1. use GPU to advance fluid by gravity
   --> slightly faster than the CPU ( ~ 1.5x )
   --> most of the time spent on the data preparing and the data transferring between CPU and GPU

2. the free-fall time-step estimation is temporarily disabled
   ~~ because so far the function "CUCAL_AdvanceByGravity" cannot NOT return the maximum acceleration

3. use GPU to evaulte the hydrodynamical CFL condition
   --> add the option "OPT__GPU_HYDRO_CFL"

4. reset some variables in the end of the function "Init_Load_Parameter"



Version 4.3.1     4/22/2009
---------------------------
1. send ( &argc, &argv ) rather than ( argc, argv ) for the function "Init"

2. add the option "OPT__POSITIVE_PRES" to compensate the negative pressure
   --> both GPU and CPU solvers support this option



Version 4.4       4/29/2009
---------------------------
1. allocate two arrays "Flu_Array1_Out" and "Flux_Array" for storing the output of the GPU fluid solver
   --> reduce the amount of data transferred from GPU to CPU
   --> both CPU and GPU solvers have adopted this modification

2. type : flux_x_t, flux_y_t, flux_z_t --> flux_t

3. slightly modify the function "Init_Load_Parameter"
   --> to work with the GNU compiler

4. add the compile options "-ip -ipo0", remove the option "-axT, -fPIC"

5. use the symbolic constant M_PI

6. use the PoissoSolver.1.2.2



Version 4.5       5/1/2009
---------------------------
1. use the PoissonSolver.1.3
   --> eight patches are grouped into one patch group for the Poisson solver
   --> only work with PATCH_SIZE == 8  && POT_GHOST_SIZE == 1 --> POT_NXT == 18
   --> all corresponding self-gravity functions are modified to work with this new strategy
   --> CPU preparing time is reduced, GPU execution time is increased
   --> higher GPU/CPU execution time ratio

2. the files "CUCAL_PoissonSolver_10cube.cu" and "CUCAL_PoissonSolver_14cube.cu" are removed

3. the internal potential is initialized as zero in both CPU and GPU Poisson solvers

4. the conditional compilation "#if ( POT_GHOST_SIZE == 1 )" is removed



Version 4.5.1     6/5/2009
---------------------------
1. fix the bug associated with the incorrect record of physical time in the individual time step
   --> the line "Time[lv] += 0.5*Delta_A;" is replaced by "Time[lv] += 0.5*da;"
   --> add the input "da" for both the shared time step as well as the individual time step

2. remove the warning about turning on the TIMING option when using the individual time step with NLEVEL > 5

3. loose dumptable criterion : 1.e-10 --> 1.e-8



Version 4.5.2     7/6/2009
---------------------------
1. add the input "OUTPUT_DT" in the parameter file
   replace the option "OPT__DUMPTABLE" by "OUTPUT__OUTPUT_MODE"
   --> support the "constant time interval" output mode



Version 4.5.3     8/17/2009
---------------------------
1. do NOT invoke the GPU kernel if the number of blocks per grid is ZERO

2. simplify the Makefile, set the path for different libraries



Version 4.5.4     8/18/2009
---------------------------
1. add the function "Output_Part" in the function "Output"
   --> also add several output parameters for the function "Output_Part"

2. add serveral checks for loading the dump table

3. re-format the file "Record__PatchCount"

4. the option "OPT__OUTPUT_POT" will also apply to the option "OPT__OUTPUT_PART"




Version 4.6       09/??/2009
----------------------------
1. work with the PoissonSolver.2.1
   (a) support POT_GHOST_SIZE = 1 - 5
   (b) perform the spatial interpolation in GPU
   (c) support central and quadratic interpolations

2. work with the FluidSolver.2.4
   --> better performance

3. no more potential correction
   --> the accuracy of the Poisson solver now depends on the size of the potential ghost zone (POT_GHOST_SIZE)

4. redefine the file "Zodiacal_Constellation"
   --> now the table records the density threshold in each level

5. use the new function "InvokeSolver" to control the concurrent execution between CPU and GPU for all GPU solvers

6. integrate all interpolation schemes into the file "Interpolation.cpp"

7. flux array "fx[2], fy[2], fz[2]" --> "flux[6]"

8. replace redundant tables by TABLE_01 and TABLE_02

9. no memory copy in the function "Refine" if the patch already exists

10. correct the bug in the function "UM_Init"
    --> also add the check routines to verify the size of the file "UM_START"

11. the function "UM_Init" now can load "1" or "NCOMP" variables

12. rename: value -> fluid

13. remove the data structure "flux_t"

14. no memory allocation for outer buffer patches




BUG      10/21/2009
-------------------
1. the function <CUCAL_FluidSolver.cu.2.4> does NOT work in the new Geisha system

2. the function <Refine> occasionally fails in the new Geisha system

3. same executable file gives different results in geisha 1-4 and geisha 5-8

4. incorrect evaluation order of individual time-step : F->G->P --> F->P->G




Unfinished easy work
--------------------
1. the TIMING_SOLVER option does NOT work for multi GPUs (because of the MPI_Barrier)

2. add a small number to the freezing speed in the fluid solver

3. more CUDA check routines (e.g. cutilSafeCall ... )

4. add input parameter for future interpolation scheme
   --> NSide = ? 26:6
   --> size of coarse-grid ghost zones




Unfinished hard work
--------------------
2. use the fine-grid data for interpolation

4. re-order the IDs of a patch group and the IDs of 26 siblings
   --> first check the efficiency of table look up

16. write the gravity term as a conservative form

25. "Const*Time[lv]" or "Const*(Time[lv]-0.5*Delta_A)" in the function "G_Evolve_AllLevel"
    ... or call the function "G_Evolve_AllLevel" twice ??
    --> now this program seems to evolve a little bit faster than the program cmt2_4

26. higher-order monotonic interpolation for the functions "Refine" and "Interpolation"
   --> refer to the RAMESES or other cosmological AMR papers
   --> e.g. MinMod ( using 26 neighbors rather than 6 neighbors ?? )

27. dynamtically allocate the variable "NLEVEL"





possible optimization
---------------------
14. [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][5] --> [5][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE] ??

16. optimize the function "Refine" for FLAG_BUFFER_SIZE != PATCH_SIZE

19. subtract the average density inside GPU

22. do NOT advance fluid variables for patches with sons and not adjacent to the coarse-fine boundaries

23. inline functions for all TABLES





