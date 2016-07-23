#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize GAMER
//-------------------------------------------------------------------------------------------------------
void Init_GAMER( int *argc, char ***argv )
{

// initialize MPI
   Init_MPI( argc, argv );


// initialize the AMR structure
   amr = new AMR_t;


// initialize the particle structure 
#  ifdef PARTICLE
   amr->Par = new Particle_t();
#  endif


// load parameters from the file "Input__Parameter"
   Init_Load_Parameter();


// set default parameters
   Init_SetDefaultParameter();


// initialize OpenMP settings
#  ifdef OPENMP
   Init_OpenMP();
#  endif


#  ifdef COMOVING
// initialize the scale factor for the cosmological simulation (it will be overwritten during restart)
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Time[lv] = A_INIT;

      amr->FluSgTime[lv][ amr->FluSg[lv] ] = Time[lv];
#     ifdef GRAVITY
      amr->PotSgTime[lv][ amr->PotSg[lv] ] = Time[lv];
#     endif
   }



// reset the gravitational constant
#  ifdef GRAVITY
   NEWTON_G = double( 3.0*OMEGA_M0/8.0/M_PI );

   if ( MPI_Rank == 0 )   
      Aux_Message( stdout, "NOTE : gravitational constant is reset to %13.7e in the comological simulations\n",
                   NEWTON_G );
#  endif
#  endif // #ifdef COMOVING
   

// initialize parameters for the parallelization (rectangular domain decomposition)
   Init_Parallelization();


#  ifdef GRAVITY
// initialize FFTW
   Init_FFTW();
#  endif


// initialize the test problem parameters
   Init_TestProb();


// initialize the external potential and acceleration parameters
// (must AFTER Init_TestProb but BEFORE CUAPI_Set_Default_GPU_Parameter)
#  ifdef GRAVITY
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      Init_ExternalAcc();

   if ( OPT__EXTERNAL_POT )
      Init_ExternalPot();
#  endif


// set the GPU ID and several GPU parameters
#  ifdef GPU
#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_Set_Default_GPU_Parameter( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP );
#  endif


// initialize the array recording the previous physical time as an arbitrary "negative" number
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Time_Prev[lv] = -__FLT_MAX__;

      amr->FluSgTime[lv][ 1-amr->FluSg[lv] ] = Time_Prev[lv];
#     ifdef GRAVITY
      amr->PotSgTime[lv][ 1-amr->PotSg[lv] ] = Time_Prev[lv];
#     endif
   }


// verify the input parameters
   Aux_Check_Parameter();


// initialize the timer function
#  ifdef TIMING
   Aux_CreateTimer();
#  endif


// load the tables of the flag criteria from the input files "Input__FLAG_XXX"
   Init_Load_FlagCriteria();


// load the dump table from the input file "Input__DumpTable"
//###NOTE: unit has not been converted into internal unit
   if ( OPT__OUTPUT_MODE == OUTPUT_USE_TABLE )
#  ifdef PARTICLE
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_TEST_ERROR || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PARTICLE )
#  else
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_TEST_ERROR || OPT__OUTPUT_BASEPS )
#  endif
   Init_Load_DumpTable();


// allocate memory for several global arrays
   Init_MemAllocate();


// initialize particles
#  ifdef PARTICLE
   switch ( amr->Par->Init )
   {
      case PAR_INIT_BY_FUNCTION:    Par_Init_Function();    break;

      case PAR_INIT_BY_RESTART:

      case PAR_INIT_BY_FILE:

      default : Aux_Error( ERROR_INFO, "unsupported particle initialization (%s = %d) !!\n", 
                           "PAR_INIT", (int)amr->Par->Init );
   }

   Par_Aux_InitCheck();

#  ifndef SERIAL
   if ( amr->Par->Init != PAR_INIT_BY_RESTART )    Par_Init_RedistributeByRectangular();
#  endif
#  endif // #ifdef PARTICLE


// initialize the AMR structure and fluid field
   switch ( OPT__INIT )
   {
      case INIT_STARTOVER:    Init_StartOver();    break;

      case INIT_RESTART :     Init_Restart();      break;

      case INIT_UM :          Init_UM();           break;

      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__INIT", OPT__INIT );
   }


// get the total number of patches at all ranks
   for (int lv=0; lv<NLEVEL; lv++)     Mis_GetTotalPatchNumber( lv );


// improve load balance    
#  ifdef LOAD_BALANCE
   LB_Init_LoadBalance( OPT__INIT == INIT_RESTART );

// fill up the data for patches that are not leaf patches (for RESTART only)
   if ( OPT__INIT == INIT_RESTART )
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
      Flu_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], NULL_INT, NULL_INT, _FLU );

      LB_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_RESTRICT, _FLU, NULL_INT );

      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_YES );
   }
#  endif // #ifdef LOAD_BALANCE


#  ifdef GRAVITY
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
//    initialize the k-space Green's function for the isolated BC.
      if ( OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();


//    evaluate the average density if it is not set yet (may already be set in Init_Restart)
      if ( AveDensity <= 0.0 )   Poi_GetAverageDensity();


//    evaluate the gravitational potential
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" ); 

      for (int lv=0; lv<NLEVEL; lv++)     
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

         Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf, USELB_YES );

         Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false );

         if ( lv > 0 )  
         Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=0; lv<NLEVEL; lv++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" ); 
   } // if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
#  endif // #ifdef GARVITY


// initialize particle acceleration
#  if ( defined PARTICLE  &&  defined STORE_PAR_ACC )
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" ); 

   const bool StoreAcc_Yes    = true;
   const bool UseStoredAcc_No = false;

   for (int lv=0; lv<NLEVEL; lv++)
   Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" ); 
#  endif


// initialize the array "MinDtInfo_Fluid" since the kernel "CUFLU_GetMaxCFL" will NOT work during initialization
   real MinDtVar_AllLv_Fluid[NLEVEL][NCOMP];
   if ( OPT__ADAPTIVE_DT )     
   {
#     if   ( MODEL == HYDRO )
      Hydro_GetMaxCFL( MinDtInfo_Fluid, MinDtVar_AllLv_Fluid );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
//    nothing to do here

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   }

// initialize the array "MinDtInfo_Gravity" since the kernel "XXX" will NOT work during initialization
#  ifdef GRAVITY
   if ( OPT__ADAPTIVE_DT )     
   {
#     if   ( MODEL == HYDRO )
      Hydro_GetMaxAcc( MinDtInfo_Gravity );

#     elif ( MODEL == MHD )
#     error : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      real MinDtVar_AllLv_Gravity[NLEVEL][3];
      ELBDM_GetMaxPot( MinDtInfo_Gravity, MinDtVar_AllLv_Gravity );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL
   }
#  endif

// initialize the array "MinDtInfo_Phase" since the kernel "XXX" will NOT work during initialization
#  if ( MODEL == ELBDM )
   real MinDtVar_AllLv_Phase[NLEVEL][3];

   if ( OPT__ADAPTIVE_DT )     
      ELBDM_GetMaxPhaseDerivative( MinDtInfo_Phase, MinDtVar_AllLv_Phase );
#  endif


} // FUNCTION : Init_GAMER
