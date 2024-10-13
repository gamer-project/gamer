

                        =================================================================
                        ||                                                             ||
                        ||    GAMER : GPU-accelerated Adaptive MEsh Refinement code    ||
                        ||                                                             ||
                        =================================================================


------------------------------------------------------------------------------------------------------------------

Version 1.0.beta6.0.0   ??/??/2012
----------------------------------
1. Support initializing FluSg and PotSg arbitrarily for both individual and
   shared time-step integrations

2. Support external force

3. Support non-periodic boundary conditions

4. Support monotonicity for all interpolation schemes

5. Support HDF5




===============
|| Bug fixes ||
===============
1. Fix the bug in "Init_Parallelization"
   --> The case "Flu_ParaBuf < InitGhostSize_Rho" will fail in the previous
       version (not enough ghost-zone data after refinement for performing
       interpolation on density)
   --> Flu_ParaBuf = MAX( Flu_ParaBuf, IntGhostSize_Rho );



Known issues         03/21/2012
-------------------------------
1. Out-of-core computation:
   (1) Failed for NDISK = 7
   (2) Failed for 2048^3 run (res != obj->u.c.nbytes)
   (3) OOC_UpdateBuffer takes too much time --> much longer than Buf_GetBufferData(send) !!
   (4) The performance improvement by using the TargetP list seems not to be good enough
2. WAF scheme does not work well in cosmological simulations
3. EXACT Riemann solver does not work well in cosmological simulations
4. Sometimes the program will become idle in Fermi GPUs
   (when the option "FERMI" is turned on)
5. OpenMP directive for grandson check is disabled in version beta5.0
   --> Please ensure that when multiple threads write to the same memory
       position, at least one of them will success.
6. OVERLAP_MPI will cause program to be terminated sometimes (in the 32-GPU
   runs) ...
   --> probably is due to bad multiple-thread support in OpenMPI 1.5.3
   --> segmentation fault when using mvapich2
7. High-resolution load-balance run (reload data N256_RefineRho4/Data_000025)
   using 8 nodes failed in Dirac (16- and 32-node run is ok).
   --> Error message :
       [[40961,1],7][btl_openib_component.c:3224:handle_wc]
       from dirac17 to: dirac20-ib error polling LP CQ with status LOCAL LENGTH ERROR
       status number 1 for wr_id 52572936 opcode 1  vendor error 105 qp_idx 3
   --> Probably due to sending too large arrays in MPI ??
8. Program hangs sometimes when using OpenMPI. Switching to MVAPICH2 on Dirac
   seems to be able to get rid of this issue !!
9. Disable OpenMP in DT__PHASE seems to make the code more stable !?
   --> check the DT__PHASE function again, is omp_nested important here??
   --> NONONO... still died !!!
   --> Something wrong in "Prepare_PatchGroupData" ??
10. Program died sometimes on Dirac during the "Gra_Advance(lv=0)" !!
   --> There seems to be something wrong in the base-leveel Poisson solver
       with FFTW !!
   --> It seems that the realloc function will fail in the file
       "LB_SFC_to_Slice.cpp"
       --> Setting "MemUnit" to a larger number to avoid calling the realloc
           function seems to make the code much stable (DON"T KNOW WHY...)
11. The memory is not properly released on the Dirac system sometimes
    (ok in Geisha and Laohu !!)



Future work (easy)   11/28/2011
-------------------------------
1.  Add a boolean parameter "Loaded" in the structure "AMR_t" to indicate whether or not the data of a targeted
    level has been loaded in the out-of-core computing.
2.  For the options "OPT__GPU_DT_CFL" and "OPT__GPU_DT_ACC", one must calculate the time-step by CPU for the
    newly created patches.
3.  Support all flag criteria in the function "Aux_Check_Refinement".
4.  Limit the maximum line width to be 80.
5.  Declare all non-negative variables as "uint".
6.  Move all global variables into the data structure (e.g., AMR_t)
7.  Use C function to grep the content in "Input__NoteScript"
    --> the system call "cat ..." will fail in the Dirac system
8.  Support PGI compiler
    (1) Add the following macros in the file Prototype.h
          #ifdef __cplusplus
          extern "C" {
          #endif
          ...
          #ifdef __cplusplus
          }
          #endif
    (2) maybe add the option "COMPILER=XXX" in the Makefile
    (3) Add the following macros for all files requiring the function "getline"
          #ifndef _GNU_SOURCE
          #define _GNU_SOURCE
          #endif
          #include <stdio.h>
9.  Define symbolic constants (or enum type) for all function input options
    --> e.g., pnew, TIMING_FUNC, Flu_Advance, ...
10. More flexible format of the file "Input__Parameter"
11. Add the error check of the return values for all function calls (e.g., MPI_Init, fftw, fopen, fclose ... ).
12. Move all #error checks in a isolated file, which should be compiled first.
14. Replace the "system" call by C library.
15. Replace all "getline" functions by the example "ReadParameter2".
17. Replace all macros by the inline functions.
18. Check the negative density and pressure during the initialization.
19. Unify the coefficient of the Poisson equation in the functions "Solver"
    and "Gra_AdvanceDt.cpp".
20. Input the name of the binary output file in Input__Parameter.
21. Limit the increase of time-step in each step by 2X.
24. Fully support "OPT__VERBOSE"
28. Dynamically allocate all XXX[MPI_Rank] arrays and all other automatic
    arrays with undetermined array size during compilation
31. Replace all "enum A" by "typedef int A" as suggested by Enzo
33. Move operations for "Stage == 0" in the function "Output_DumpData" to
    "Init_Output"
    --> Make "Output_DumpData" as simple as possible so that one can comment
        it out without causing any problem
34. Add the refinement criteria of "vorticity"
36. Give all arrays for CPU/GPU solvers more meaningful names ...
37. Add the Jeans length refinement criterion
    --> flag if Jeans length <= Constant*(cell size), (Constant = 4 for
        example)
    --> refer to "MODELING COLLAPSE AND ACCRETION IN TURBULENT GAS CLOUDS:
        IMPLEMENTATION AND COMPARISON OF SINK PARTICLES IN AMR AND SPH"
        by Christoph Federrath, 2010, ApJ, 713, 269
    --> refer to "Magnetic Fields in Population III Star Formation"
        by Matthew J. Turk et al.
39. Classify all global variables --> add the headers "HydroVar.h,
    MHDVar.h, ELBDMVar.h, Common.h" and their correspondiner structures
    --> no more global variables ??
    --> or just no more individual variables (all variables should be members of
        structures)
40. rename "patch" as "AMR", and "ptr" as "Patch"
41. Add ENFORCE_POSITIVE in Flu_FixUp and after performing interpolations
42. Check negative density and energy after interpolation in the debug mode
43. Input reference energy for OPT__CK_CONSERVATION during restart
44. Develop a tool to calculate dump table
45. Ensure that "long" is put properly everywhere. For example
    int  a = 10240;
    long b = a*a*a;        wrong
    long b = (long)a*a*a;  correct
46. Add the error checking routine after every memory allocation
47. Add the macro "GAMER_VERSION" and store it in the output file
48. Support using an arbitrary number of GPUs if LOAD_BALANCE is on
   (currently the code only supports this feature during RESTART)
49. Replace "MPI_DOUBLE" and "MPI_FLOAT" by "GAMER_MPI_REAL"
50. Data output
    (1) output different "partial" data at the same time
    (2) support projection
    (3) support file prefix and suffix
    (4)
51. Output reference solution recorded in the function "Aux_Check_Conservation"
    so that they can be reused after program restart



Future work (hard)   07/16/2011
-------------------------------
1.  The GPU/CPU fluid solver can save the output results in an arbitrary Sg.
2.  Let each "half step" contain both the forward and backward sweepings for the fluid solver.
3.  Use the fine-grid data for interpolation.
4.  Remove the constraint that eight nearby patches must reside in the same MPI process. Also let the patch
    indices in one patch group be arbitrary.
5.  Ensure the momentum and energy conservation even in a self-gravity system.
    --> Refer to Springel's AREPO paper, Ui-Li Pen's "A HIGH-RESOLUTION
        ADAPTIVE MOVING MESH HYDRODYNAMIC ALGORITHM" paper in 1998, Athena
        code, Collins, D. C.'s EnzoMHD paper ...
    -->Integrate the hydro and gravity solvers.
    --> Refer to Collins, D. C.'s EnzoMHD paper (as a start), Athena code
6.  Let the CFL condition be ajustable during one global step.
7.  No global variables
8.  Support both GPU_DT_CFL and GPU_DT_ACC
9.  Complete the GAMER_DEBUG mode
11. Adopt dual energy formalism (also advance entropy)
    --> Refer to Collins, D. C.'s EnzoMHD paper (as a start)
    --> Additional references: FLASH manual, AREPO paper, Cosmological PPM
        paper
12. Positively conservative hydro scheme
13.  Adopt higher-order monotonic interpolation schemes
    --> refer to the RAMESES or other cosmological AMR papers
    --> e.g. MinMod ( using 26 neighbors rather than 6 neighbors ?? )
14. Complete the pancake collapse test.
15. Add different boundary conditions (reflecting, flow-in, flow-out,
    user-defined ...)
    --> refer to Athena
    --> Athena also contains the 3D FFT solver for open boundaries !!!
    --> Also refer to Gadjet 2 for the Dirichlet B.C. using FFT
16. Support sending/receiving data large than 4 GB in MPI
17. Write our own Hilbert curve functions
18. Remove the constraint that dt(lv) = dt(lv-1)/2
    --> integrate individual time-step and shared time-step integrations
19. Remove (replace) all 1D coordinates which are not suitable for extremely
    large simulations
20. Try the visualization tools "FlashView" and "QuickFlash"
21. Make the refinement criteria more flexibie --> check magnitude or gradient
    for any components
22. Add particles
    --> Refer to all cosmological AMR codes
    --> Refer to Athena + particle paper
    --> Learn different CIC scheme (linear [ref ??], triangular-shaped [ref ??]), tri-cubic [ref Lekien & Marsden, 2005, Int.                                     J. Numer. Methods Egn., 63, 455 ~ already downloaded])



Optimizations        07/16/2011
-------------------------------
1.  Optimize the performance of the function "Flag" for "FLAG_BUFFER_SIZE != PATCH_SIZE".
2.  Do NOT advance fluid variables for patches with sons and not adjacent to the coarse-fine boundaries.
3.  Load fluid and potential separately with different TargetP_t list for the out-of-core computing.
    --> The OOC_Load must also support loading without the mode "_STRC"
4.  Re-dump the level "lv" data after refinement only if the TargetP_t list is changed (for the OOC computing).
    --> Even if the TargetP_t list is changed, we can just exchange the positions of the changed patches.
5.  Construct new sending and receiving data list for updating the buffer data after the "fix-up" operation.
    --> We don't need to load and dump the NoSonNoCFB data in the OOC computing.
6.  Optimize the data loading (restart procedure, HDF5 ?)
7.  Perform the spatial interpolation in GPU for the hydro solver
    (refer to GAMER.1.0.beta4.0.t30.bk__Interpolation_in_GPU)
8.  Reduce the MPI data transfer for the case "FLU_PARA_BUF < RHO_INT"
9.  Perform the MPI data transfer sequentially in x, y, and z directions, so
    that the extra MPI calls can be avoided (the total amount of transferred
    data is not changed).
    --> Particularly important for the out-of-core computation since that the
        amount of data stored and loaded at a time is increased.
10. Use the "__restrict__" qualifier.
11. OpenMP in all OOC functions
12. Support OpenMP in the functions "Refine", "FindFather" ...
13. Optimize the OpenMP performance in "Flu_FixUp"
    --> not all patches require the fix-up operation
    --> maybe "dynamic" scheduling will be better
14. Remove redundant table look-ups by precalculate results and store them in
    arrays like "LoopStart[26][3], LoopEnd[26][3]"
15. Support QuickSort with OpenMP
16. Support HeapSort with OpenMP
17. Support "OVERLAP_MPI" for the rectangular domain decomposition
18. Optimize the performance of "OVERLAP_MPI"
19. Optimize the load-balance efficiency for the cases that the box sizes are
    non-cubic and/or non-power-of-two.
20. Try to further optimize the patch distribution for LOAD_BALACNE in order to
    minimize the MPI time
21. Try "-prec-sqrt=false" in new CUDA version to check if it is stable now
22. Try "CUDA Thrust library" for sorting in the load-balance computing



Priority             01/17/2012
-------------------------------
1. Add the constant gravity for Kris
2. Add different kinds of boundary conditions



===================
||     ELBDM     ||
===================

Future work          12/22/2011
-------------------------------
1. Add the Jeans length refinement criterion
2. Check the Schrodinger+AMR related work
3. However about restricting "PHASE" instead of real/imag parts ?
4. Do not need to calculate time-step from phase rotation for patches
   surrounded by sons


Serious issue        11/28/2011
-------------------------------
1. Mass conservation
   --> There are three sources that will deteriorate the mass conservation
       (1) kinematic solver      --> not solved yet
       (2) coarse-fine boundary  --> not solved yet
       (3) grid refinement       --> solved by using conservative
           interpolation on density and (a) interpolation on phase or
           (b) rescaling real/imag parts by the refined density
2. Allowed time-step is extremely small ...
3. Is there any better solution than "interpolation on phase" to ensure the
   smoothness of momentum after grid refinement ?

