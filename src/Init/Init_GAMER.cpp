#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY
extern void (*Init_ExternalAcc_Ptr)();
extern void (*Init_ExternalPot_Ptr)();
#endif




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


// load runtime parameters
   Init_Load_Parameter();


// set code units
   Init_Unit();


// reset parameters --> must be called after Init_Unit()
   Init_ResetParameter();


// initialize OpenMP settings
#  ifdef OPENMP
   Init_OpenMP();
#  endif


// initialize yt inline analysis
#  ifdef SUPPORT_LIBYT
   YT_Init( *argc, *argv );
#  endif


// initialize parameters for the parallelization (rectangular domain decomposition)
   Init_Parallelization();


#  ifdef GRAVITY
// initialize FFTW
   Init_FFTW();
#  endif


// initialize the test problem parameters
   Init_TestProb();


// initialize the external potential and acceleration parameters
// --> must be called AFTER Init_TestProb() but BEFORE CUAPI_Set_Default_GPU_Parameter()
// --> these function pointers point to "Init_ExternalAcc()" and "Init_ExternalPot()" by default
//     but may be overwritten by various test problem initializers
#  ifdef GRAVITY
   if (  ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&  Init_ExternalAcc_Ptr != NULL  )
      Init_ExternalAcc_Ptr();

   if ( OPT__EXTERNAL_POT  &&  Init_ExternalPot_Ptr != NULL )
      Init_ExternalPot_Ptr();
#  endif


// initialize settings for the passive variables
// --> must be called BEFORE CUAPI_Set_Default_GPU_Parameter()
   Init_PassiveVariable();


// set the GPU ID and several GPU parameters
// --> must be called AFTER Init_ExternalAcc(), Init_ExternalPot(), and Init_PassiveVariable()
#  ifdef GPU
#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_Set_Default_GPU_Parameter( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP );
#  endif


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
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PAR_TEXT )
#  else
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS )
#  endif
   Init_Load_DumpTable();


// initialize memory pool
   if ( OPT__MEMORY_POOL )    Init_MemoryPool();


// allocate memory for several global arrays
   Init_MemAllocate();


// initialize particles
#  ifdef PARTICLE
   switch ( amr->Par->Init )
   {
      case PAR_INIT_BY_FUNCTION:    Par_Init_ByFunction();  break;

      case PAR_INIT_BY_RESTART:                             break;   // nothing to do here for the restart mode

      case PAR_INIT_BY_FILE:        Par_Init_ByFile();      break;

      default : Aux_Error( ERROR_INFO, "unsupported particle initialization (%s = %d) !!\n",
                           "PAR_INIT", (int)amr->Par->Init );
   }

   if ( amr->Par->Init != PAR_INIT_BY_RESTART )
   {
      Par_Aux_InitCheck();

#     ifndef SERIAL
      Par_LB_Init_RedistributeByRectangular();
#     endif
   }
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

// we don't have to redistribute all patches during the RESTART process since we already did that in Init_Restart()
// --> but note that Init_Restart() does NOT consider load-balance weighting of particles
// we don't have enough information to calculate the load-balance weighting of particles when
// calling LB_Init_LoadBalance() for the first time
// --> must disable particle weighting (by setting ParWeight==0.0) first
// must not reset load-balance variables (i.e., must adopt ResetLB_No) when calling LB_Init_LoadBalance() for the first time
// since we MUST NOT overwrite IdxList_Real and IdxList_Real_IdxList set during the restart process
   const double ParWeight_Zero   = 0.0;
   const bool   Redistribute_Yes = true;
   const bool   Redistribute_No  = false;
   const bool   ResetLB_Yes      = true;
   const bool   ResetLB_No       = false;

   LB_Init_LoadBalance( (OPT__INIT==INIT_RESTART)?Redistribute_No:Redistribute_Yes, ParWeight_Zero, ResetLB_No );

// redistribute patches again if we want to take into account the load-balance weighting of particles
#  ifdef PARTICLE
   if ( amr->LB->Par_Weight > 0.0 )
   LB_Init_LoadBalance( Redistribute_Yes, amr->LB->Par_Weight, ResetLB_Yes );
#  endif


// record the initial weighted load-imbalance factor
   if ( OPT__RECORD_LOAD_BALANCE )  LB_EstimateLoadImbalance();


// fill up the data for patches that are not leaf patches (for RESTART only)
// --> It's for bitwise consistency between load-balance and non-load-balance runs
// --> Should be deprecated (or removed) after adding the makefile option "BITWISE_REPRODUCIBILITY",
//     which will always apply data restriction before dumping data
   if ( OPT__INIT == INIT_RESTART )
   for (int lv=NLEVEL-2; lv>=0; lv--)
   {
      Flu_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], NULL_INT, NULL_INT, _TOTAL );

      LB_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_RESTRICT, _TOTAL, NULL_INT );

      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, Flu_ParaBuf, USELB_YES );
   }
#  endif // #ifdef LOAD_BALANCE


#  ifdef GRAVITY
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
//    initialize the k-space Green's function for the isolated BC.
      if ( OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();


//    evaluate the initial average density if it is not set yet (may already be set in Init_Restart)
      if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();


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
   if ( OPT__ADAPTIVE_DT )
   {
#     if   ( MODEL == HYDRO )
      real MinDtVar_AllLv_Fluid[NLEVEL][5];
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
