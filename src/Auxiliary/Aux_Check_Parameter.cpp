#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Parameter
// Description :  Check the initial parameter setting
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Parameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ... \n" );


// general errors
// =======================================================================================
#  if ( PATCH_SIZE%2 != 0 )
#     error : ERROR : PATCH_SIZE must be an even number !!
#  endif

#  if ( defined TIMING_SOLVER  &&  !defined TIMING )
#     error : ERROR : option TIMING_SOLVER must work with the option TIMING !!
#  endif

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP, the macro "_OPENMP" is NOT defined !!
#  endif

#  if ( defined OVERLAP_MPI  &&  !defined LOAD_BALANCE )
#     error : ERROR : option OVERLAP_MPI must work with the option LOAD_BALANCE !!
#  endif

#  if ( !defined GRAVITY  &&  defined UNSPLIT_GRAVITY )
#     error : ERROR : please turn off UNSPLIT_GRAVITY when GRAVITY is off !!
#  endif

#  if ( defined UNSPLIT_GRAVITY  &&  MODEL != HYDRO )
#     error : ERROR : currently UNSPLIT_GRAVITY is only supported in HYDRO !!
#  endif

#  ifdef SERIAL
   int NRank = 1;
#  else
   int NRank, MPI_Thread_Status, MPI_Init_Status;

   MPI_Initialized( &MPI_Init_Status );
   if ( MPI_Init_Status == false )  Aux_Error( ERROR_INFO, "MPI_Init has not been called !!\n" );

   MPI_Query_thread( &MPI_Thread_Status );
   if ( OPT__OVERLAP_MPI  &&  MPI_Thread_Status == MPI_THREAD_SINGLE )
      Aux_Error( ERROR_INFO, "option \"%s\" is NOT supported if the level of MPI thread support == %s\n",
                 "OPT__OVERLAP_MPI", "MPI_THREAD_SINGLE" );

   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
#  endif

   if ( MPI_NRank != NRank )
      Aux_Error( ERROR_INFO, "MPI_NRank (%d) != MPI_Comm_size (%d) !!\n", MPI_NRank, NRank );

   if ( NX0_TOT[0] <= 0  ||  NX0_TOT[1] <= 0  ||  NX0_TOT[2] <= 0 )
      Aux_Error( ERROR_INFO, "incorrect number of base-level grids --> NX0_TOT[0/1/2] = [%d,%d,%d] !!\n",
                 NX0_TOT[0], NX0_TOT[1], NX0_TOT[2] );

   if ( NX0_TOT[0]%PS2 != 0  ||  NX0_TOT[1]%PS2 != 0  ||  NX0_TOT[2]%PS2 != 0 )
      Aux_Error( ERROR_INFO, "number of base-level patches in each direction must be \"%s\" !!\n",
                 "a multiple of TWO" );

   if ( END_STEP < 0  &&  OPT__INIT != INIT_RESTART )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d\" [>=0] !!\n", "END_STEP", END_STEP );

   if ( END_T < 0.0  &&  OPT__INIT != INIT_RESTART )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %14.7e\" [>=0] !!\n", "END_T", END_T );

#  ifdef LOAD_BALANCE
   if ( OPT__INIT != INIT_RESTART )
#  endif
   if ( NX0_TOT[0]%(PS2*MPI_NRank_X[0]) != 0  ||  NX0_TOT[1]%(PS2*MPI_NRank_X[1]) != 0  ||
        NX0_TOT[2]%(PS2*MPI_NRank_X[2]) != 0 )
      Aux_Error( ERROR_INFO, "number of base-level patches in each direction in one rank must be \"%s\" !!\n",
                 "a multiple of TWO" );

   if ( OPT__INIT != INIT_STARTOVER  &&  OPT__INIT != INIT_RESTART  &&  OPT__INIT != INIT_UM )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__INIT = %d\" [1/2/3] !!\n", OPT__INIT );

#  ifdef LOAD_BALANCE
   if ( OPT__INIT != INIT_RESTART )
#  endif
   if ( MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank )
      Aux_Error( ERROR_INFO, "MPI_NRank_X[0]*MPI_NRank_X[1]*MPI_NRank_X[2] != MPI_NRank !!\n" );

   if ( MAX_LEVEL < 0  ||  MAX_LEVEL > NLEVEL-1 )
      Aux_Error( ERROR_INFO, "MAX_LEVEL (%d) is not within the accepted range [0 ... NLEVEL-1] !!\n", MAX_LEVEL );

   if ( FLAG_BUFFER_SIZE > PATCH_SIZE )   Aux_Error( ERROR_INFO, "FLAG_BUFFER_SIZE > PATCH_SIZE !!\n" );

   if ( REGRID_COUNT <= 0 )   Aux_Error( ERROR_INFO, "REGRID_COUNT (%d) <= 0 !!\n", REGRID_COUNT );

   if ( OPT__OUTPUT_MODE != OUTPUT_CONST_STEP  &&  OPT__OUTPUT_MODE != OUTPUT_CONST_DT  &&
        OPT__OUTPUT_MODE != OUTPUT_USE_TABLE )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_MODE = %d\" [1/2/3] !!\n", OPT__OUTPUT_MODE );

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_STEP  &&  OUTPUT_STEP <= 0 )
      Aux_Error( ERROR_INFO, "OUTPUT_STEP <= 0 !!\n" );

   if ( OPT__OUTPUT_MODE == OUTPUT_CONST_DT  &&  OUTPUT_DT <= 0 )
      Aux_Error( ERROR_INFO, "OUTPUT_DT <= 0 !!\n" );

   if ( OPT__RESTART_HEADER != RESTART_HEADER_CHECK  &&  OPT__RESTART_HEADER != RESTART_HEADER_SKIP )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__RESTART_HEADER = %d\" [0/1] !!\n",
                 OPT__RESTART_HEADER );

   if ( OPT__OUTPUT_TOTAL != OUTPUT_TOTAL_NONE  &&  OPT__OUTPUT_TOTAL != OUTPUT_FORMAT_HDF5  &&
        OPT__OUTPUT_TOTAL != OUTPUT_FORMAT_CBINARY )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_TOTAL = %d\" [0/1/2] !!\n", OPT__OUTPUT_TOTAL );

#  ifndef SUPPORT_HDF5
   if ( OPT__OUTPUT_TOTAL == OUTPUT_FORMAT_HDF5 )
      Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for OPT__OUTPUT_TOTAL == 1 !!\n" );
#  endif

   if ( OPT__OUTPUT_PART != OUTPUT_PART_NONE  &&  OPT__OUTPUT_PART != OUTPUT_DIAG  &&
        OPT__OUTPUT_PART != OUTPUT_XY  &&  OPT__OUTPUT_PART != OUTPUT_YZ  &&  OPT__OUTPUT_PART != OUTPUT_XZ  &&
        OPT__OUTPUT_PART != OUTPUT_X   &&  OPT__OUTPUT_PART != OUTPUT_Y   &&  OPT__OUTPUT_PART != OUTPUT_Z )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__OUTPUT_PART = %d\" [0 ~ 6] !!\n", OPT__OUTPUT_PART );

   if (  ( OPT__OUTPUT_PART == OUTPUT_YZ  ||  OPT__OUTPUT_PART == OUTPUT_Y  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_X < 0.0  ||  OUTPUT_PART_X >= amr->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_X (out of range [0<=X<%lf]) !!\n", amr->BoxSize[0] );

   if (  ( OPT__OUTPUT_PART == OUTPUT_XZ  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Z )  &&
         ( OUTPUT_PART_Y < 0.0  ||  OUTPUT_PART_Y >= amr->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Y (out of range [0<=Y<%lf]) !!\n", amr->BoxSize[1] );

   if (  ( OPT__OUTPUT_PART == OUTPUT_XY  ||  OPT__OUTPUT_PART == OUTPUT_X  ||  OPT__OUTPUT_PART == OUTPUT_Y )  &&
         ( OUTPUT_PART_Z < 0.0  ||  OUTPUT_PART_Z >= amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect OUTPUT_PART_Z (out of range [0<=Z<%lf]) !!\n", amr->BoxSize[2] );

   if (  OPT__OUTPUT_PART == OUTPUT_DIAG  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "option \"%s\" only works with CUBIC domain !!\n",
                 "OPT__OUTPUT_PART == 7 (OUTPUT_DIAG)" );

   if (  OPT__OUTPUT_BASEPS  &&  ( NX0_TOT[0] != NX0_TOT[1] || NX0_TOT[0] != NX0_TOT[2] )  )
      Aux_Error( ERROR_INFO, "option \"%s\" only works with CUBIC domain !!\n", "OPT__OUTPUT_BASEPS" );

   if (  OPT__INIT == INIT_UM  &&  ( OPT__UM_START_LEVEL < 0 || OPT__UM_START_LEVEL > NLEVEL-1 )  )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__UM_START_LEVEL = %d\" [0 ... NLEVEL-1] !!\n",
                 OPT__UM_START_LEVEL );

   if (  OPT__INIT == INIT_UM  &&  ( OPT__UM_START_NVAR < 1 || OPT__UM_START_NVAR > NCOMP )  )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__UM_START_NVAR = %d\" [1 ... NCOMP] !!\n",
                 OPT__UM_START_NVAR );

   if ( OPT__INIT == INIT_UM  &&  OPT__UM_START_NVAR != 1  &&  OPT__UM_FACTOR_5OVER3 )
      Aux_Error( ERROR_INFO, "OPT__UM_FACTOR_5OVER3 only works when OPT__UM_START_NVAR == 1 !!\n" );

   if ( OPT__CK_REFINE  &&  !OPT__FLAG_RHO )
      Aux_Error( ERROR_INFO, "currently the check \"%s\" must work with \"%s\" !!\n",
                 "OPT__CK_REFINE", "OPT__FLAG_RHO == 1" );

#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   if (  ( OPT__FLAG_LOHNER_DENS || OPT__FLAG_LOHNER_ENGY || OPT__FLAG_LOHNER_PRES )  &&   Flu_ParaBuf < 2  )
      Aux_Error( ERROR_INFO, "Lohner error estimator does NOT work when Flu_ParaBuf (%d) < 2 !!\n", Flu_ParaBuf );
#  elif ( MODEL == ELBDM )
   if (  OPT__FLAG_LOHNER_DENS  &&  Flu_ParaBuf < 2  )
      Aux_Error( ERROR_INFO, "Lohner error estimator does NOT work when Flu_ParaBuf (%d) < 2 !!\n", Flu_ParaBuf );
#  else
#  error : ERROR : unsupported MODEL !!
#  endif

   if (  OPT__FLAG_LOHNER_FORM != LOHNER_FLASH1     &&  OPT__FLAG_LOHNER_FORM != LOHNER_FLASH2  &&
         OPT__FLAG_LOHNER_FORM != LOHNER_FORM_INV1  &&  OPT__FLAG_LOHNER_FORM != LOHNER_FORM_INV2  )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n", "OPT__FLAG_LOHNER_FORM", OPT__FLAG_LOHNER_FORM );

#  ifdef GPU
   if ( OPT__GPUID_SELECT < 0  &&  OPT__GPUID_SELECT != -1  &&  OPT__GPUID_SELECT != -2 )
   {
#     ifdef LAOHU
      if ( OPT__GPUID_SELECT != -3 )
#     endif
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n", "OPT__GPUID_SELECT", OPT__GPUID_SELECT );
   }
#  endif

   if ( INIT_SUBSAMPLING_NCELL < 0 )   Aux_Error( ERROR_INFO, "INIT_SUBSAMPLING_NCELL < 0 !!\n" );

#  ifdef SERIAL
   if ( MPI_NRank != 1 )   Aux_Error( ERROR_INFO, "\"MPI_NRank != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[0] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[0] != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[1] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[1] != 1\" in the serial code !!\n" );

   if ( MPI_NRank_X[2] != 1 )    Aux_Error( ERROR_INFO, "\"MPI_NRank_X[2] != 1\" in the serial code !!\n" );
#  endif // #ifdef SERIAL

#  ifndef OVERLAP_MPI
   if ( OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT turned on in the Makefile for the option \"%s\" !!\n",
                 "OVERLAP_MPI", "OPT__OVERLAP_MPI" );
#  endif

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC    &&  OPT__BC_FLU[f] != BC_FLU_OUTFLOW  &&
        OPT__BC_FLU[f] != BC_FLU_REFLECTING  &&  OPT__BC_FLU[f] != BC_FLU_USER        )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__BC_FLU[%d] = %d\" [1/2/3/4] !!\n", f, OPT__BC_FLU[f] );

#  if ( MODEL != HYDRO )
   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_OUTFLOW  ||  OPT__BC_FLU[f] == BC_FLU_REFLECTING )
      Aux_Error( ERROR_INFO, "outflow and reflecting boundary conditions (OPT__BC_FLU=2/3) only work with HYDRO !!\n" );
#  endif

   if (  ( OPT__BC_FLU[0] == BC_FLU_PERIODIC || OPT__BC_FLU[1] == BC_FLU_PERIODIC || OPT__BC_FLU[2] == BC_FLU_PERIODIC ||
           OPT__BC_FLU[3] == BC_FLU_PERIODIC || OPT__BC_FLU[4] == BC_FLU_PERIODIC || OPT__BC_FLU[5] == BC_FLU_PERIODIC   ) &&
         ( OPT__BC_FLU[0] != BC_FLU_PERIODIC || OPT__BC_FLU[1] != BC_FLU_PERIODIC || OPT__BC_FLU[2] != BC_FLU_PERIODIC ||
           OPT__BC_FLU[3] != BC_FLU_PERIODIC || OPT__BC_FLU[4] != BC_FLU_PERIODIC || OPT__BC_FLU[5] != BC_FLU_PERIODIC   )   )
      Aux_Error( ERROR_INFO, "currently the periodic BC cannot be mixed with non-periodic BC. !!\n" );

   if ( INT_MONO_COEFF < 1.0  ||  INT_MONO_COEFF > 4.0 )
      Aux_Error( ERROR_INFO, "INT_MONO_COEFF (%14.7e) is not within the correct range [1...4] !!\n", INT_MONO_COEFF );

#  ifndef TIMING
   if ( OPT__TIMING_MPI )  Aux_Error( ERROR_INFO, "OPT__TIMING_MPI only works when TIMING is on !!\n" );
#  endif

#  ifndef INDIVIDUAL_TIMESTEP
   if ( OPT__INT_TIME )    Aux_Error( ERROR_INFO, "OPT__INT_TIME only works when INDIVIDUAL_TIMESTEP is on !!\n" );
#  endif

   if ( OPT__PATCH_COUNT < 0  ||  OPT__PATCH_COUNT > 2 )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__PATCH_COUNT = %d\" [0/1/2] !!\n", OPT__PATCH_COUNT );

   if ( OPT__REUSE_MEMORY < 0  ||  OPT__REUSE_MEMORY > 2 )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__REUSE_MEMORY = %d\" [0/1/2] !!\n", OPT__REUSE_MEMORY );

   if ( OPT__MEMORY_POOL  &&  !OPT__REUSE_MEMORY )
      Aux_Error( ERROR_INFO, "please turn on OPT__REUSE_MEMORY for OPT__MEMORY_POOL !!\n" );


// general warnings
// =======================================================================================
#  ifdef OVERLAP_MPI
#     warning : NOTE : make sure to link with multi-thread supported MPI and FFTW for the option "OVERLAP_MPI"
#  endif

#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {
      if ( OMP_NTHREAD != omp_get_num_threads() )
         Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) != omp_get_num_threads (%d) at MPI_Rank %d !!\n",
                      OMP_NTHREAD, omp_get_num_threads(), MPI_Rank );
   }
#  endif


   if ( MPI_Rank == 0 ) {

   if ( !OPT__OUTPUT_TOTAL  &&  !OPT__OUTPUT_PART  &&  !OPT__OUTPUT_TEST_ERROR  &&  !OPT__OUTPUT_BASEPS )
#  ifdef PARTICLE
   if ( !OPT__OUTPUT_PAR_TEXT )
#  endif
      Aux_Message( stderr, "WARNING : all output options are turned off --> no data will be output !!\n" );

   if ( OPT__CK_REFINE )
      Aux_Message( stderr, "WARNING : \"%s\" check may fail due to the proper-nesting constraint !!\n",
                   "OPT__CK_REFINE" );

   if ( !OPT__INIT_RESTRICT )
      Aux_Message( stderr, "WARNING : option OPT__INIT_RESTRICT is NOT turned on !!\n" );

#  ifdef TIMING_SOLVER
   Aux_Message( stderr, "WARNING : option \"TIMING_SOLVER\" will disable the concurrent execution\n" );
   Aux_Message( stderr, "          between GPU and CPU and hence will decrease the overall performance !!\n" );
#  endif

   if ( OPT__RESTART_HEADER == RESTART_HEADER_SKIP )
   {
      Aux_Message( stderr, "WARNING : to skip the header check during restart, you must make sure that \n" );
      Aux_Message( stderr, "          everything is set correctly in both Input__Parameter and Makefile !!\n" );
   }

   if ( OPT__CK_REFINE )
      Aux_Message( stderr, "WARNING : currently the check \"%s\" only works with \"%s\" !!\n",
                   "OPT__CK_REFINE", "OPT__FLAG_RHO == 1" );

   bool Flag = ( OPT__FLAG_RHO  ||  OPT__FLAG_RHO_GRADIENT  ||  OPT__FLAG_LOHNER_DENS  ||  OPT__FLAG_USER );
#  if ( MODEL == HYDRO )
   Flag |= OPT__FLAG_PRES_GRADIENT;
   Flag |= OPT__FLAG_LOHNER_ENGY;
   Flag |= OPT__FLAG_LOHNER_PRES;
#  endif
#  if ( MODEL == ELBDM )
   Flag |= OPT__FLAG_ENGY_DENSITY;
#  endif
#  ifdef PARTICLE
   Flag |= OPT__FLAG_NPAR_PATCH;
   Flag |= OPT__FLAG_NPAR_CELL;
#  endif

   if ( !Flag )
   {
      Aux_Message( stderr, "WARNING : all flag criteria are turned off --> no refinement will be " );
      Aux_Message( stderr, "performed !!\n" );
   }

   if ( OPT__OVERLAP_MPI )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" is still experimental and is not fully optimized !!\n",
                   "OPT__OVERLAP_MPI" );

#     ifdef OPENMP
      omp_set_nested( true );

      if ( !omp_get_nested() )
         Aux_Message( stderr, "WARNING : OpenMP nested parallelism is NOT supported for the option \"%s\" !!\n",
                      "OPT__OVERLAP_MPI" );

      omp_set_nested( false );
#     else
      Aux_Message( stderr, "WARNING : OpenMP is NOT turned on for the option \"%s\" !!\n", "OPT__OVERLAP_MPI" );
#     endif

   } // if ( OPT__OVERLAP_MPI )

   if (  ( OPT__BC_FLU[0] == BC_FLU_USER || OPT__BC_FLU[1] == BC_FLU_USER || OPT__BC_FLU[2] == BC_FLU_USER ||
           OPT__BC_FLU[3] == BC_FLU_USER || OPT__BC_FLU[4] == BC_FLU_USER || OPT__BC_FLU[5] == BC_FLU_USER   ) &&
         ( OPT__BC_FLU[0] != BC_FLU_USER || OPT__BC_FLU[1] != BC_FLU_USER || OPT__BC_FLU[2] != BC_FLU_USER ||
           OPT__BC_FLU[3] != BC_FLU_USER || OPT__BC_FLU[4] != BC_FLU_USER || OPT__BC_FLU[5] != BC_FLU_USER   )   )
      Aux_Message( stderr, "WARNING : corner cells may not be well defined when mixing user-defined BC with others !!\n" );

   if ( OPT__TIMING_BARRIER )
      Aux_Message( stderr, "WARNING : option \"%s\" may deteriorate performance (especially if %s is on) ...\n",
                   "OPT__TIMING_BARRIER", "OPT__OVERLAP_MPI" );

   if ( OPT__TIMING_BARRIER  &&  !OPT__TIMING_BALANCE )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" is on, but the time waiting for other ranks will NOT be included in individual timers ...\n",
                   "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "          --> the sum of individual timer may be less than the total elapsed time due to load imbalance ...\n" );
   }

#  ifdef TIMING
   if ( !OPT__TIMING_BARRIER )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" is off for TIMING\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "          --> Some timing results (especially MPI and particle routines) may be less accurate due to load imbalance ...\n" );
   }
#  endif

   if (  ( OPT__TIMING_BALANCE || OPT__TIMING_MPI )  &&  !OPT__TIMING_BARRIER  )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" is off for OPT__TIMING_BALANCE/OPT__TIMING_MPI\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "          --> Some timing results (especially MPI and particle routines) may be less accurate due to load imbalance ...\n" );
   }

#  ifdef PARTICLE
   if ( OPT__TIMING_BALANCE )
   {
      Aux_Message( stderr, "WARNING : option \"%s\" does NOT work well for particle routines\n", "OPT__TIMING_BARRIER" );
      Aux_Message( stderr, "          --> Because many particle routines call MPI_Barrier implicitly\n" );
      Aux_Message( stderr, "              (so as Gra_AdvanceDt when PARTICLE is on)\n" );
   }
#  endif

   if (  OPT__INIT == INIT_UM  &&  OPT__UM_FACTOR_5OVER3  )
   {
#     if ( MODEL != HYDRO  &&  MODEL != MHD  &&  MODEL != ELBDM )
      Aux_Message( stderr, "WARNING : OPT__UM_FACTOR_5OVER3 has no effect in the current model !!\n" );
#     endif

#     ifndef COMOVING
      Aux_Message( stderr, "WARNING : COMOVING is NOT defined for the option OPT__UM_FACTOR_5OVER3 !!\n" );
#     endif

      Aux_Message( stderr, "REMINDER : please make sure that \"background density ~= 1.0\" for OPT__UM_FACTOR_5OVER3\n" );
   }

#  if ( GRA_GHOST_SIZE == 0  &&  defined STORE_POT_GHOST )
   Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless when GRA_GHOST_SIZE == 0 !!\n" );
#  endif

#  if ( !defined GRAVITY  &&  defined STORE_POT_GHOST )
   Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless when GRAVITY is off !!\n" );
#  endif

   } // if ( MPI_Rank == 0 )



// load balance
// =======================================================================================
#ifdef LOAD_BALANCE

// errors
// ------------------------------
#  ifdef SERIAL
#     error : ERROR : options LOAD_BALANCE and SERIAL should NOT be turned on at the same time !!
#  endif

#  if ( LOAD_BALANCE != HILBERT )
#     error : ERROR : currently GAMER only supports "LOAD_BALANCE == HILBERT" !!
#  endif

// for sending fluid data fixed by coarse-fine fluxes correctly
   if ( OPT__FIXUP_FLUX  &&  Flu_ParaBuf >= PATCH_SIZE )
      Aux_Error( ERROR_INFO, "we must have \"%s\" for \"%s\" in LOAD_BALANCE !!\n",
                 "Flu_ParaBuf < PATCH_SIZE", "OPT__FIXUP_FLUX" );

// ensure that the variable "PaddedCr1D" will not overflow
   const int Padded              = 1<<NLEVEL;
   const int BoxNScale_Padded[3] = { amr->BoxScale[0]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[1]/PATCH_SIZE + 2*Padded,
                                     amr->BoxScale[2]/PATCH_SIZE + 2*Padded };
   if (  (double)BoxNScale_Padded[0]*(double)BoxNScale_Padded[1]*(double)BoxNScale_Padded[2]-1.0 > (double)__ULONG_MAX__  )
   {
      Aux_Message( stderr, "ERROR : PaddedCr1D can overflow !!\n" );
      Aux_Message( stderr, "    --> The currently maximum allowed 1-D resolution is about 2^24 (for PS1=8)\n" );
      Aux_Message( stderr, "    --> Please either set NLEVEL to a smaller number or disable LOAD_BALANCE !!\n" );
      MPI_Exit();
   }


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( NX0_TOT[0] != NX0_TOT[1]  ||  NX0_TOT[0] != NX0_TOT[2] )
   {
      Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-cubic simulation box\n" );
      Aux_Message( stderr, "          (NX0_TOT[0] != NX0_TOT[1] and/or NX0_TOT[0] != NX0_TOT[2]) !!\n" );
   }

   for (int d=0; d<3; d++)
   {
      if ( NX0_TOT[d] & (NX0_TOT[d]-1) )
      {
         Aux_Message( stderr, "WARNING : LOAD_BALANCE has NOT been fully optimized for non-power-of-two " );
         Aux_Message( stderr, "simulation box (NX0_TOT[%d] = %d) !!\n", d, NX0_TOT[d] );
      }
   }

   } // if ( MPI_Rank == 0 )

#else // #ifdef LOAD_BALANCE ... else ...

   if ( OPT__OVERLAP_MPI  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : currently the option \"%s\" is useful only in LOAD_BALANCE !!\n",
                   "OPT__OVERLAP_MPI" );

#endif // #ifdef LOAD_BALANCE ... else ...



// comoving frame
// =======================================================================================
#ifdef COMOVING

// errors
// ------------------------------
   if ( OMEGA_M0 < 0.0  ||  OMEGA_M0 > 1.0 )
      Aux_Error( ERROR_INFO, "incorrect OMEGA_M0 (0.0 <= OMEGA_M0 <= 1.0) !!\n" );

   if ( HUBBLE0 <= 0.0 )
      Aux_Error( ERROR_INFO, "HUBBLE0 = %14.7e <= 0.0 !!\n", HUBBLE0 );

   if ( A_INIT > 1.0  ||  A_INIT < 0.0 )
      Aux_Error( ERROR_INFO, "incorrect A_INIT (0.0 <= A_INIT <= 1.0) !!\n" );

#  if   ( MODEL == HYDRO )
   if ( fabs(GAMMA-5.0/3.0) > 1.e-4 )
      Aux_Error( ERROR_INFO, "GAMMA must be equal to 5.0/3.0 in cosmological simuluations !!\n" );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif


// warnings
// ------------------------------
#  ifndef GRAVITY
   if ( MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : option \"%s\" is useless if the option \"%s\" is NOT turned on !!\n",
                   "COMOVING", "GRAVITY" );
#  endif

#endif // COMOVING



// fluid solver in all models
// =======================================================================================

// errors
// ------------------------------
   if ( Flu_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Flu_ParaBuf, PATCH_SIZE );

   if ( GPU_NSTREAM < 1 )  Aux_Error( ERROR_INFO, "GPU_NSTREAM < 1 !!\n" );

   if ( FLU_GPU_NPGROUP % GPU_NSTREAM != 0 )    Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP%%GPU_NSTREAM != 0 !!\n" );

   if ( OPT__FLU_INT_SCHEME != INT_MINMOD3D  &&  OPT__FLU_INT_SCHEME != INT_MINMOD1D  &&
        OPT__FLU_INT_SCHEME != INT_VANLEER   &&  OPT__FLU_INT_SCHEME != INT_CQUAD     &&
        OPT__FLU_INT_SCHEME != INT_QUAD      &&  OPT__FLU_INT_SCHEME != INT_CQUAR     &&
        OPT__FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );

   if ( OPT__REF_FLU_INT_SCHEME != INT_MINMOD3D  &&  OPT__REF_FLU_INT_SCHEME != INT_MINMOD1D  &&
        OPT__REF_FLU_INT_SCHEME != INT_VANLEER   &&  OPT__REF_FLU_INT_SCHEME != INT_CQUAD     &&
        OPT__REF_FLU_INT_SCHEME != INT_QUAD      &&  OPT__REF_FLU_INT_SCHEME != INT_CQUAR     &&
        OPT__REF_FLU_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );

   if ( OPT__FIXUP_FLUX  &&  !amr->WithFlux )
      Aux_Error( ERROR_INFO, "%s is turned on but amr->WithFlux is off !!\n", "OPT__FIXUP_FLUX" );

#  if ( NLEVEL > 1 )
   int Trash_RefFlu, NGhost_RefFlu;
   Int_Table( OPT__REF_FLU_INT_SCHEME, Trash_RefFlu, NGhost_RefFlu );
   if ( Flu_ParaBuf < NGhost_RefFlu )
      Aux_Error( ERROR_INFO, "Flu_ParaBuf (%d) < NGhost_RefFlu (%d) --> refinement will fail !!\n",
                 Flu_ParaBuf, NGhost_RefFlu );
#  endif

   if ( OPT__RESET_FLUID  &&   OPT__OVERLAP_MPI )
      Aux_Error( ERROR_INFO, "\"%s\" is NOT supported for \"%s\" !!\n", "OPT__OVERLAP_MPI", "OPT__RESET_FLUID" );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( OPT__CK_FLUX_ALLOCATE  &&  !amr->WithFlux )
      Aux_Message( stderr, "WARNING : option \"%s\" is useless since no flux is required !!\n",
                   OPT__CK_FLUX_ALLOCATE );

   if ( DT__FLUID < 0.0  ||  DT__FLUID > 1.0 )
      Aux_Message( stderr, "WARNING : DT__FLUID (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__FLUID );

   if ( DT__FLUID_INIT < 0.0  ||  DT__FLUID_INIT > 1.0 )
      Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__FLUID_INIT );

   if ( OPT__RESET_FLUID  &&   OPT__INIT == INIT_UM )
      Aux_Message( stderr, "WARNING : \"%s\" will NOT be applied to the input uniform data !!\n", "OPT__RESET_FLUID" );

   } // if ( MPI_Rank == 0 )



// fluid solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( NCOMP != 5 )
#     error : ERROR : NCOMP != 5 in HYDRO !!
#  endif

#  if ( FLU_SCHEME != RTVD  &&  FLU_SCHEME != MHM  &&  FLU_SCHEME != MHM_RP  &&  FLU_SCHEME != CTU  &&  \
        FLU_SCHEME != WAF )
#     error : ERROR : unsupported hydro scheme in the makefile !!
#  endif

#  if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == RTVD || FLU_SCHEME == WAF )  )
#     error : ERROR : RTVD and WAF do not support UNSPLIT_GRAVITY !!
#  endif

#  if ( defined LR_SCHEME  &&  LR_SCHEME != PLM  &&  LR_SCHEME != PPM )
#     error : ERROR : unsupported data reconstruction scheme (PLM/PPM) !!
#  endif

#  if ( defined RSOLVER  &&  RSOLVER != EXACT  &&  RSOLVER != ROE  &&  RSOLVER != HLLE  &&  RSOLVER != HLLC )
#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!
#  endif

#  if ( NPASSIVE < 0 )
#     error : ERROR : incorrect number of NPASSIVE !!
#  endif

#  if ( defined CHECK_INTERMEDIATE  &&  CHECK_INTERMEDIATE != EXACT  &&  CHECK_INTERMEDIATE != HLLE  &&  \
        CHECK_INTERMEDIATE != HLLC )
#     error : ERROR : unsupported option in CHECK_INTERMEDIATE (EXACT/HLLE/HLLC) !!
#  endif

#  if ( defined MIN_PRES_DENS  &&  defined MIN_PRES )
#    error : ERROR : MIN_PRES_DENS and MIN_PRES in the file "CUFLU.h" cannot be turned on at the same time !!
#  endif

#  if ( NPASSIVE < 0 )
   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "currently the passively advected scalars can only work with the periodic BC !!\n" );
#  endif

   if ( OPT__CK_NEGATIVE < 0  ||  OPT__CK_NEGATIVE > 3 )
      Aux_Error( ERROR_INFO, "unsupported parameter \"%s = %d\" !!\n", "OPT__CK_NEGATIVE", OPT__CK_NEGATIVE );

   if ( OPT__CORR_UNPHY )
   {
      if ( OPT__CORR_UNPHY_SCHEME != RSOLVER_ROE  &&  OPT__CORR_UNPHY_SCHEME != RSOLVER_HLLC  &&
           OPT__CORR_UNPHY_SCHEME != RSOLVER_HLLE )
         Aux_Error( ERROR_INFO, "unsupported parameter \"%s = %d\" !!\n", "OPT__CORR_UNPHY_SCHEME", OPT__CORR_UNPHY_SCHEME );

#     if ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == WAF )
         Aux_Error( ERROR_INFO, "RTVD and WAF fluid schemes do not support the option \"OPT__CORR_UNPHY\" !!\n" );
#     endif
   }


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( FLU_SCHEME == MHM  &&  LR_SCHEME == PPM )
#     warning : WARNING : PPM is not recommended for the MHM scheme !!
      Aux_Message( stderr, "WARNING : PPM is not recommended for the MHM scheme !!\n" );
#  endif

#  if ( FLU_SCHEME == MHM_RP  &&  LR_SCHEME == PPM )
#     warning : WARNING : PPM is not recommended for MHM_RP scheme !!
      Aux_Message( stderr, "WARNING : PPM is not recommended for the MHM_RP scheme !!\n" );
#  endif

#  if ( defined RSOLVER  &&  RSOLVER == EXACT )
#     warning : WARNING : exact RSOLVER is not recommended since the vacuum solution has not been implemented
      Aux_Message( stderr, "WARNING : exact Riemann solver is not recommended since the vacuum solution " );
      Aux_Message( stderr,           "has not been implemented !!\n" );
#  endif

#  if ( defined CHAR_RECONSTRUCTION  &&  defined GRAVITY )
#     warning : WARNING : option "CHAR_RECONSTRUCTION" is less robust and can cause negative density/pressure !!
      Aux_Message( stderr, "WARNING : option \"CHAR_RECONSTRUCTION\" is less robust and can cause negative " );
      Aux_Message( stderr,           "density/pressure !!\n" );
#  endif

#  if ( !defined MIN_PRES  &&  !defined MIN_PRES_DENS )
#     warning : WARNING : options "MIN_PRES and MIN_PRES_DENS" are both turned off --> negative pressure may happen !!
      Aux_Message( stderr, "WARNING : option \"MIN_PRES and MIN_PRES_DENS\" are both turned off --> negative " );
      Aux_Message( stderr,           "pressure may happen !!\n" );
#  endif

#  ifdef MIN_PRES
   Aux_Message( stderr, "WARNING : MIN_PRES (%14.7e) is on --> please make sure that this value is reasonable !!\n",
                MIN_PRES );
#  endif

#  ifdef MIN_PRES_DENS
   Aux_Message( stderr, "WARNING : MIN_PRES_DENS (%14.7e) is on --> please make sure that this value is reasonable !!\n",
                MIN_PRES_DENS );
#  endif

   if ( !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_FLUX" );

   if ( !OPT__FIXUP_RESTRICT )
      Aux_Message( stderr, "WARNING : option \"%s\" is disabled in HYDRO !!\n", "OPT__FIXUP_RESTRICT" );

   if ( OPT__CK_FLUX_ALLOCATE  &&  !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option %s is useless since %s is off !!\n",
                   OPT__CK_FLUX_ALLOCATE, OPT__FIXUP_FLUX );

   } // if ( MPI_Rank == 0 )


// check for MHM/MHM_RP/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

// errors
// ------------------------------
#  if ( LR_SCHEME == PPM )
   if ( OPT__LR_LIMITER == EXTPRE )
      Aux_Error( ERROR_INFO, "currently the PPM reconstruction does not support the \"%s\" limiter\n",
                 "extrema-preserving" );
#  endif

   if ( OPT__LR_LIMITER != VANLEER  &&  OPT__LR_LIMITER != GMINMOD  &&  OPT__LR_LIMITER != ALBADA  &&
        OPT__LR_LIMITER != EXTPRE   &&  OPT__LR_LIMITER != VL_GMINMOD )
      Aux_Error( ERROR_INFO, "unsupported data reconstruction limiter (OPT__LR_IMITER = %d) !!\n",
                 OPT__LR_LIMITER );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if (  ( OPT__LR_LIMITER == GMINMOD || OPT__LR_LIMITER == EXTPRE || OPT__LR_LIMITER == VL_GMINMOD )  &&
         ( MINMOD_COEFF < 1.0 || MINMOD_COEFF > 2.0 )  )
      Aux_Message( stderr, "WARNING : MinMod limiter coefficient (%14.7e) is not within the range [1...2] !!\n",
                   MINMOD_COEFF );

   if ( OPT__LR_LIMITER == EXTPRE  &&  EP_COEFF < 1.0 )
      Aux_Message( stderr, "WARNING : coefficient of the extrema-preserving limiter (%14.7e) < 1.0 !!\n",
                   EP_COEFF );

   } // if ( MPI_Rank == 0 )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )


// check for MHM/CTU
// ------------------------------
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 2 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 2  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 2", "MHM/CTU scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM/CTU scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == CTU )


// check for MHM_RP
// ------------------------------
#  if ( FLU_SCHEME == MHM_RP )

#  if ( LR_SCHEME == PLM )
   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER == EXTPRE  &&  FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PLM reconstruction + EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE < 3 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );

   if ( OPT__LR_LIMITER != EXTPRE  &&  FLU_GHOST_SIZE > 3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 3", "MHM_RP scheme + PLM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PLM )

#  if ( LR_SCHEME == PPM )
   if ( FLU_GHOST_SIZE < 4 )
      Aux_Error( ERROR_INFO, "please set \"%s\" for \"%s\" !!\n",
                 "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );

   if ( FLU_GHOST_SIZE > 4  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : please set \"%s\" in \"%s\" for higher performance !!\n",
                   "FLU_GHOST_SIZE = 4", "MHM_RP scheme + PPM reconstruction + non-EXTPRE limiter" );
#  endif // #if ( LR_SCHEME == PPM )

#  endif // #if ( FLU_SCHEME == MHM_RP )


// check for WAF
// ------------------------------
#  if ( FLU_SCHEME == WAF )

#  if ( RSOLVER == HLLE  ||  RSOLVER == HLLC )
#     error : ERROR : currently the WAF scheme does not support HLLE/HLLC Riemann solvers
#  endif

#  if ( FLU_GHOST_SIZE != 2 )
#     error : ERROR : please set FLU_GHOST_SIZE = 2 for the WAF scheme !!
#  endif

   if ( OPT__WAF_LIMITER != WAF_SUPERBEE  &&  OPT__WAF_LIMITER != WAF_VANLEER  &&
        OPT__WAF_LIMITER != WAF_ALBADA    &&  OPT__WAF_LIMITER != WAF_MINBEE      )
      Aux_Error( ERROR_INFO, "unsupported WAF flux limiter (%d) !!\n", OPT__WAF_LIMITER );
#  endif // if ( FLU_SCHEME == WAF )


// check for RTVD
// ------------------------------
#  if ( FLU_SCHEME == RTVD )

#  if ( FLU_GHOST_SIZE != 3 )
#     error : ERROR : please set FLU_GHOST_SIZE = 3 for the relaxing TVD scheme !!
#  endif

#  endif // if ( FLU_SCHEME == RTVD )



// fluid solver in MHD
// =======================================================================================
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

// errors
// ------------------------------

// warnings
// ------------------------------



// fluid solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )

// errors
// ------------------------------
#  if ( NCOMP != 3 )
#     error : ERROR : NCOMP != 3 in ELBDM !!
#  endif

#  if ( FLU_NIN != 2 )
#     error : ERROR : FLU_NIN != 2 in ELBDM !!
#  endif

#  if ( FLU_NOUT != 3 )
#     error : ERROR : FLU_NOUT != 3 in ELBDM !!
#  endif

#  ifdef QUARTIC_SELF_INTERACTION
#  ifndef GRAVITY
#     error : ERROR : currently QUARTIC_SELF_INTERACTION must work with GRAVITY !!
#  endif

#  ifdef COMOVING
#     error : ERROR : currently QUARTIC_SELF_INTERACTION does not work with COMOVING yet !!
#  endif
#  endif // ifdef QUARTIC_SELF_INTERACTION

   if ( ELBDM_MASS <= 0.0 )
      Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "ELBDM_MASS", ELBDM_MASS );

   if ( ELBDM_PLANCK_CONST <= 0.0 )
      Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "ELBDM_PLANCK_CONST", ELBDM_PLANCK_CONST );

   if ( ELBDM_ETA <= 0.0 )
      Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "ELBDM_ETA", ELBDM_ETA );

   if ( OPT__INT_PHASE  &&  OPT__FLU_INT_SCHEME == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" when OPT__INT_PHASE is on !!\n",
                 "OPT__FLU_INT_SCHEME", OPT__FLU_INT_SCHEME );

   if ( OPT__INT_PHASE  &&  OPT__FLU_INT_SCHEME == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme \"%s = %d\" when OPT__INT_PHASE is on !!\n",
                 "OPT__REF_FLU_INT_SCHEME", OPT__REF_FLU_INT_SCHEME );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_REFLECTING  ||  OPT__BC_FLU[f] == BC_FLU_OUTFLOW )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__BC_FLU[%d] = %d\" [1/4] !!\n", f, OPT__BC_FLU[f] );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( !ELBDM_TAYLOR3_AUTO  &&  ELBDM_TAYLOR3_COEFF < 1.0/8.0 )
      Aux_Message( stderr, "WARNING : ELBDM_TAYLOR3_COEFF (%13.7e) < 0.125 is unconditionally unstable !!\n",
                   ELBDM_TAYLOR3_COEFF );

#  ifdef LAPLACIAN_4TH
   const double dt_fluid_max = 3.0*M_PI/16.0;
#  else
   const double dt_fluid_max = 0.25*M_PI;
#  endif
   if ( DT__FLUID > dt_fluid_max )
      Aux_Message( stderr, "WARNING : DT__FLUID (%13.7e) > %13.7e is unconditionally unstable (even with %s) !!\n",
                   DT__FLUID, dt_fluid_max, "ELBDM_TAYLOR3_AUTO" );

   if ( DT__FLUID_INIT > dt_fluid_max )
      Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%13.7e) > %13.7e is unconditionally unstable (even with %s) !!\n",
                   DT__FLUID_INIT, dt_fluid_max, "ELBDM_TAYLOR3_AUTO" );

   if ( !ELBDM_TAYLOR3_AUTO )
   {
//    stability limit for ELBDM_TAYLOR3_COEFF == 1.0/6.0
#     ifdef LAPLACIAN_4TH
      const double dt_fluid_max_normal = SQRT(27.0)*M_PI/32.0;
#     else
      const double dt_fluid_max_normal = SQRT(3.0)*M_PI/8.0;
#     endif

      if ( DT__FLUID > dt_fluid_max_normal  &&  ELBDM_TAYLOR3_COEFF <= 1.0/6.0 )
      {
         Aux_Message( stderr, "WARNING : DT__FLUID (%13.7e) > stability limit (%13.7e) for ELBDM_TAYLOR3_COEFF <= 1/6\n",
                      DT__FLUID, dt_fluid_max_normal );
         Aux_Message( stderr, "          --> Please either (a) set ELBDM_TAYLOR3_COEFF (%13.7e) > 1/6\n",
                      ELBDM_TAYLOR3_COEFF );
         Aux_Message( stderr, "                            (b) set DT__FLUID smaller (c) turn on ELBDM_TAYLOR3_AUTO\n" );
      }

      if ( DT__FLUID_INIT > dt_fluid_max_normal  &&  ELBDM_TAYLOR3_COEFF <= 1.0/6.0 )
      {
         Aux_Message( stderr, "WARNING : DT__FLUID_INIT (%13.7e) > stability limit (%13.7e) for ELBDM_TAYLOR3_COEFF <= 1/6\n",
                      DT__FLUID_INIT, dt_fluid_max_normal );
         Aux_Message( stderr, "          --> Please either (a) set ELBDM_TAYLOR3_COEFF (%13.7e) > 1/6\n",
                      ELBDM_TAYLOR3_COEFF );
         Aux_Message( stderr, "                            (b) set DT__FLUID_INIT smaller (c) turn on ELBDM_TAYLOR3_AUTO\n" );
      }
   }

   if ( DT__PHASE < 0.0  ||  DT__PHASE > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PHASE (%13.7e) is not within the normal range [0...1] !!\n",
                   DT__PHASE );

   if ( OPT__CK_FLUX_ALLOCATE  &&  !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option %s is useless since %s is off !!\n",
                   OPT__CK_FLUX_ALLOCATE, OPT__FIXUP_FLUX );

#  ifdef CONSERVE_MASS
   if ( !OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option \"%s\" is disabled in ELBDM even though CONSERVE_MASS is on !!\n",
                   "OPT__FIXUP_FLUX" );
#  else
   if ( OPT__FIXUP_FLUX )
      Aux_Message( stderr, "WARNING : option %s is useless in ELBDM if CONSERVE_MASS is off !!\n", OPT__FIXUP_FLUX );
#  endif
   } // if ( MPI_Rank == 0 )

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL



// Poisson and Gravity solvers
// =======================================================================================
#ifdef GRAVITY

// errors
// ------------------------------
#  if ( POT_SCHEME != SOR  &&  POT_SCHEME != MG )
#     error : ERROR : unsupported Poisson solver in the makefile (SOR/MG) !!
#  endif

#  if ( POT_GHOST_SIZE <= GRA_GHOST_SIZE )
      #error : ERROR : POT_GHOST_SIZE <= GRA_GHOST_SIZE !!
#  endif

#  if ( POT_GHOST_SIZE < 1 )
#     error : ERROR : POT_GHOST_SIZE < 1 !!
#  endif

#  ifdef GPU
#     if ( POT_GHOST_SIZE > 5 )
#        error : ERROR : POT_GHOST_SIZE must <= 5 for the GPU Poisson solver !!
#     endif

   if (  ( OPT__GRAVITY_TYPE == GRAVITY_SELF || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&  PATCH_SIZE != 8  )
      Aux_Error( ERROR_INFO, "PATCH_SIZE must == 8 for the GPU Poisson solver !!\n" );
#  endif // GPU

#  ifndef LOAD_BALANCE
   if ( NX0_TOT[2]%MPI_NRank != 0 )
   {
      Aux_Message( stderr, "ERROR : NX0_TOT[2] %% MPI_NRank != 0 !!\n" );
      Aux_Message( stderr, "--> All MPI processes must have the same number of cells in the z direction for " );
      Aux_Message( stderr,     "the slab decomposition in FFTW 2.1.5 if LOAD_BALANCE is off !!\n" );
      MPI_Exit();
   }
#  endif

   if ( Pot_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Pot_ParaBuf, PATCH_SIZE );

   if ( Rho_ParaBuf > PATCH_SIZE )
      Aux_Error( ERROR_INFO, "Rho_ParaBuf (%d) > PATCH_SIZE (%d) !!\n", Rho_ParaBuf, PATCH_SIZE );

#  if ( POT_SCHEME == SOR )
   if ( SOR_MIN_ITER < 3 )    Aux_Error( ERROR_INFO, "SOR_MIN_ITER < 3 !!\n" );
#  endif

   if ( POT_GPU_NPGROUP % GPU_NSTREAM != 0 )
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP %% GPU_NSTREAM != 0 !!\n" );

   if ( OPT__POT_INT_SCHEME != INT_CQUAD  &&  OPT__POT_INT_SCHEME != INT_QUAD )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__POT_INT_SCHEME", OPT__POT_INT_SCHEME );

   if ( OPT__RHO_INT_SCHEME != INT_MINMOD3D  &&  OPT__RHO_INT_SCHEME != INT_MINMOD1D  &&
        OPT__RHO_INT_SCHEME != INT_VANLEER   &&  OPT__RHO_INT_SCHEME != INT_CQUAD     &&
        OPT__RHO_INT_SCHEME != INT_QUAD      &&  OPT__RHO_INT_SCHEME != INT_CQUAR     &&
        OPT__RHO_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__RHO_INT_SCHEME", OPT__RHO_INT_SCHEME );

   if ( OPT__GRA_INT_SCHEME != INT_MINMOD3D  &&  OPT__GRA_INT_SCHEME != INT_MINMOD1D  &&
        OPT__GRA_INT_SCHEME != INT_VANLEER   &&  OPT__GRA_INT_SCHEME != INT_CQUAD     &&
        OPT__GRA_INT_SCHEME != INT_QUAD      &&  OPT__GRA_INT_SCHEME != INT_CQUAR     &&
        OPT__GRA_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__GRA_INT_SCHEME", OPT__GRA_INT_SCHEME );

   if ( OPT__REF_POT_INT_SCHEME != INT_MINMOD3D  &&  OPT__REF_POT_INT_SCHEME != INT_MINMOD1D  &&
        OPT__REF_POT_INT_SCHEME != INT_VANLEER   &&  OPT__REF_POT_INT_SCHEME != INT_CQUAD     &&
        OPT__REF_POT_INT_SCHEME != INT_QUAD      &&  OPT__REF_POT_INT_SCHEME != INT_CQUAR     &&
        OPT__REF_POT_INT_SCHEME != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported input parameter \"%s = %d\" !!\n",
                 "OPT__REF_POT_INT_SCHEME", OPT__REF_POT_INT_SCHEME );

#  if ( NLEVEL > 1 )
   int Trash_RefPot, NGhost_RefPot;
   Int_Table( OPT__REF_POT_INT_SCHEME, Trash_RefPot, NGhost_RefPot );
   if ( Pot_ParaBuf < NGhost_RefPot )
      Aux_Error( ERROR_INFO, "Pot_ParaBuf (%d) < NGhost_RefPot (%d) --> refinement will fail !!\n",
                 Pot_ParaBuf, NGhost_RefPot );
#  endif

   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   if (  ( OPT__BC_FLU[0] == BC_FLU_PERIODIC && OPT__BC_POT != BC_POT_PERIODIC )  ||
         ( OPT__BC_FLU[0] != BC_FLU_PERIODIC && OPT__BC_POT == BC_POT_PERIODIC )    )
      Aux_Error( ERROR_INFO, "periodic BC must be applied to both fluid and self-gravity solvers at the same time !!\n" );

   if ( OPT__BC_POT != BC_POT_PERIODIC  &&  OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__BC_POT = %d\" [1/2] !!\n", OPT__BC_POT );

   if ( OPT__GRAVITY_TYPE != GRAVITY_SELF  &&  OPT__GRAVITY_TYPE != GRAVITY_EXTERNAL  &&  OPT__GRAVITY_TYPE != GRAVITY_BOTH )
      Aux_Error( ERROR_INFO, "unsupported option \"%s = %d\" [1/2/3] !!\n", "OPT__GRAVITY_TYPE", OPT__GRAVITY_TYPE );

   if (  OPT__EXTERNAL_POT  &&  ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  )
   {
      if ( MPI_Rank == 0 )
      {
         Aux_Message( stderr, "ERROR : OPT__EXTERNAL_POT does not work with \"OPT__GRAVITY_TYPE == 2/3 (EXTERNAL/BOTH)\" !!\n" );
         Aux_Message( stderr, "        --> HYDRO : please use OPT__GRAVITY_TYPE = 2/3 only\n" );
         Aux_Message( stderr, "            ELBDM : please use OPT__EXTERNAL_POT only\n" );
      }

      MPI_Exit();
   }


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  if ( POT_SCHEME == MG  &&  PATCH_SIZE <= 8 )
   {
      Aux_Message( stderr, "WARNING : multigrid scheme gives lower performance than SOR for " );
      Aux_Message( stderr, "PATCH_SIZE <= 8 and hence is not recommended !!\n" );
   }
#  endif

   if ( DT__GRAVITY < 0.0  ||  DT__GRAVITY > 1.0 )
      Aux_Message( stderr, "WARNING : DT__GRAVITY (%14.7e) is not within the normal range [0...1] !!\n",
                   DT__GRAVITY );

   if (  ( OPT__GRAVITY_TYPE == GRAVITY_SELF || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&  OPT__BC_POT == BC_POT_ISOLATED  )
   {
      Aux_Message( stderr, "WARNING : currently the patches adjacent to the simulation boundary are NOT allowed to be\n" );
      Aux_Message( stderr, "          refined if the self-gravity with the isolated BC is chosen !!\n" );
   }

   if ( OPT__EXTERNAL_POT  &&  OPT__OUTPUT_POT )
      Aux_Message( stderr, "WARNING : currently the output potential does NOT include the external potential !!\n" );

   } // if ( MPI_Rank == 0 )



// gravity solver in HYDRO
// =======================================================================================
#  if   ( MODEL == HYDRO )

// errors
// ------------------------------
#  if ( GRA_GHOST_SIZE < 1 )
#     error : ERROR : GRA_GHOST_SIZE must >= 1
#  endif

   if ( OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 1 )
      Aux_Error( ERROR_INFO, "option \"%s\" requires \"%s\" !!\n",
                 "OPT__GRA_P5_GRADIENT", "GRA_GHOST_SIZE == 2" );

#  ifdef UNSPLIT_GRAVITY
   if ( OPT__GRA_P5_GRADIENT &&  USG_GHOST_SIZE == 1 )
      Aux_Error( ERROR_INFO, "option \"%s\" requires \"%s\" for UNSPLIT_GRAVITY !!\n",
                 "OPT__GRA_P5_GRADIENT", "USG_GHOST_SIZE == 2" );
#  endif

   if ( OPT__EXTERNAL_POT )   Aux_Error( ERROR_INFO, "OPT__EXTERNAL_POT is NOT supported in HYDRO --> use external gravity !!\n" );


// warnings
// ------------------------------
   if ( MPI_Rank == 0 ) {

#  ifndef STORE_POT_GHOST
   if ( !OPT__GRA_P5_GRADIENT  &&  GRA_GHOST_SIZE == 2 )
   {
      Aux_Message( stderr, "WARNING : \"%s\" is useless if the option \"%s\" is NOT turned on !!\n",
                   "GRA_GHOST_SIZE == 2", "OPT__GRA_P5_GRADIENT" );
   }
#  endif

   if ( GRA_GHOST_SIZE > 2 )  Aux_Message( stderr, "WARNING : \"GRA_GHOST_SIZE > 2\" !?\n" );

   } // if ( MPI_Rank == 0 )



// gravity solver in MHD
// =======================================================================================
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

// errors
// ------------------------------

// warnings
// ------------------------------



// gravity solver in ELBDM
// =======================================================================================
#  elif ( MODEL == ELBDM )

// errors
// ------------------------------
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      Aux_Error( ERROR_INFO, "ELBDM does NOT support external gravity (OPT__GRAVITY_TYPE == 2/3) --> use external potential !!\n" );


// warnings
// ------------------------------
#  if ( GRA_GHOST_SIZE != 0  &&  !defined STORE_POT_GHOST )
#     warning : WARNING : GRA_GHOST_SIZE != 0 in ELBDM when STORE_POT_GHOST is off !!
#  endif


#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#endif // GRAVITY


// particle
// =======================================================================================
#ifdef PARTICLE

// errors
// ------------------------------
#  ifndef GRAVITY
#     error : ERROR : currently PARTICLE must work with GRAVITY !!
#  endif

#  ifdef COMOVING
#     error : ERROR : currently PARTICLE dost NOT support COMOVING !!
#  endif

#  if ( !defined SERIAL  &&  !defined LOAD_BALANCE )
#     error : ERROR : PARTICLE must work with either SERIAL or LOAD_BALANCE !!
#  endif

   if ( OPT__INIT != INIT_RESTART )
   {
      if ( amr->Par->Init == PAR_INIT_BY_RESTART )    Aux_Error( ERROR_INFO, "PAR_INIT == RESTART but OPT__INIT != RESTART !!\n" );

      if ( amr->Par->NPar_Active_AllRank < 0 )
         Aux_Error( ERROR_INFO, "Total number of particles in all MPI ranks = %ld < 0 !!\n",
                    amr->Par->NPar_Active_AllRank );

      if ( amr->Par->NPar_AcPlusInac < 0  ||  amr->Par->NPar_AcPlusInac > amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "Incorrect total number of particles in MPI rank %d = %ld !!\n",
                    MPI_Rank, amr->Par->NPar_AcPlusInac );
   }

   if ( amr->Par->Init < 1  ||  amr->Par->Init > 3 )
      Aux_Error( ERROR_INFO, "unsupported option \"amr->Par->Init = %d\" [1/2/3] !!\n", amr->Par->Init );

   if ( amr->Par->Interp < 1  ||  amr->Par->Interp > 3 )
      Aux_Error( ERROR_INFO, "unsupported option \"amr->Par->Interp = %d\" [1/2/3] !!\n", amr->Par->Interp );

   if ( amr->Par->Integ < 1  ||  amr->Par->Integ > 2 )
      Aux_Error( ERROR_INFO, "unsupported option \"amr->Par->Integ = %d\" [1/2] !!\n", amr->Par->Integ );

#  ifndef STORE_PAR_ACC
   if ( amr->Par->SyncDump )
      Aux_Error( ERROR_INFO, "please turn on STORE_PAR_ACC for amr->Par->SyncDump (PAR_SYNC_DUMP) !!\n" );
#  endif

#  ifndef STORE_POT_GHOST
   if ( amr->Par->ImproveAcc )
      Aux_Error( ERROR_INFO, "please turn on STORE_POT_GHOST for amr->Par->ImproveAcc (PAR_IMPROVE_ACC) !!\n" );
#  endif

   if ( amr->Par->ImproveAcc  &&  amr->Par->Interp == 1 )
      Aux_Error( ERROR_INFO, "amr->Par->ImproveAcc does NOT work with amr->Par->Interp == 1 (NGP) !!\n" );

   if ( DT__PARVEL < 0.0 )
      Aux_Error( ERROR_INFO, "DT__PARVEL (%14.7e) is not within the normal range [>=0] !!\n", DT__PARVEL );

#  ifdef STORE_PAR_ACC
   if ( DT__PARACC < 0.0 )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) is not within the normal range [>=0] !!\n", DT__PARACC );
#  else
   if ( DT__PARACC != 0.0 )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) is NOT supported when STORE_PAR_ACC is off !!\n", DT__PARACC );
#  endif

   if ( OPT__FLAG_NPAR_PATCH < 0  ||  OPT__FLAG_NPAR_PATCH > 2 )
      Aux_Error( ERROR_INFO, "unsupported option \"OPT__FLAG_NPAR_PATCH = %d\" [0/1/2] !!\n", OPT__FLAG_NPAR_PATCH );

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC )
   for (int d=0; d<3; d++)
   {
      if ( NX0_TOT[d]/PS2 == 1 )
         Aux_Error( ERROR_INFO, "\"%s\" does NOT work for NX0_TOT[%d] = 2*PATCH_SIZE when periodic BC is adopted !!\n",
                    "Par_MassAssignment", d );
   }

   if ( OPT__PARTICLE_COUNT < 0  ||  OPT__PARTICLE_COUNT > 2 )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__PARTICLE_COUNT = %d\" [0/1/2] !!\n", OPT__PARTICLE_COUNT );

   if ( OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_NONE  &&  OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_PAR_ONLY  &&
        OPT__OUTPUT_PAR_DENS != PAR_OUTPUT_DENS_TOTAL )
      Aux_Error( ERROR_INFO, "incorrect option \"OPT__OUTPUT_PAR_DENS = %d\" [0/1/2] !!\n", OPT__OUTPUT_PAR_DENS );


// warning
// ------------------------------
   if ( MPI_Rank == 0 ) {

   if ( DT__PARVEL > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PARVEL (%13.7e) is not within the normal range [<=1.0] !!\n", DT__PARVEL );

   if ( DT__PARACC > 1.0 )
      Aux_Message( stderr, "WARNING : DT__PARACC (%13.7e) is not within the normal range [<=1.0] !!\n", DT__PARACC );

   if ( OPT__OVERLAP_MPI )
      Aux_Message( stderr, "WARNING : PARTICLE does not support OPT__OVERLAP_MPI !!\n" );

#  ifdef STORE_POT_GHOST
   if ( !amr->Par->ImproveAcc )
      Aux_Message( stderr, "WARNING : STORE_POT_GHOST is useless if amr->Par->ImproveAcc is off !!\n" );
#  endif

   if ( OPT__GRA_P5_GRADIENT )
      Aux_Message( stderr, "WARNING : option \"%s\" is not supported for updating particles !!\n",
                 "OPT__GRA_P5_GRADIENT" );

   } // if ( MPI_Rank == 0 )

#else // #ifdef PARTICLE

#  ifdef STORE_POT_GHOST
   Aux_Message( stderr, "WARNING : currently STORE_POT_GHOST is useless if PARTICLE is off !!\n" );
#  endif

#endif // PARTICLE


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Aux_Check_Parameter ... done\n" );

} // FUNCTION : Aux_Check_Parameter
