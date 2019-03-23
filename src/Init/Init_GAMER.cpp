#include "GAMER.h"

extern void (*Init_User_Ptr)();
#ifdef PARTICLE
extern void (*Par_Init_ByFunction_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                        real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                        real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                        real *AllAttribute[PAR_NATT_TOTAL] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize GAMER
//-------------------------------------------------------------------------------------------------------
void Init_GAMER( int *argc, char ***argv )
{

// initialize MPI
#  ifdef SERIAL
   MPI_Rank  = 0;
   MPI_NRank = 1;
#  else
   Init_MPI( argc, argv );
#  endif


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


// initialize Grackle
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_ACTIVATE )    Grackle_Init();
#  endif


// initialize parameters for the parallelization (rectangular domain decomposition)
   Init_Parallelization();


#  ifdef GRAVITY
// initialize FFTW
   Init_FFTW();
#  endif


// initialize the test problem parameters
   Init_TestProb();


// initialize all fields and particle attributes
// --> Init_Field() must be called BEFORE CUAPI_Set_Default_GPU_Parameter()
   Init_Field();
#  ifdef PARTICLE
   Par_Init_Attribute();
#  endif


// initialize the external potential and acceleration parameters
// --> must be called AFTER Init_TestProb()
#  ifdef GRAVITY
   Init_ExternalAccPot();
#  endif


// set the GPU ID and several GPU parameters
// --> must be called AFTER Init_Field()
#  ifdef GPU
#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
#  ifndef SUPPORT_GRACKLE
   int CHE_GPU_NPGROUP = NULL_INT;
#  endif
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_Set_Default_GPU_Parameter( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP, CHE_GPU_NPGROUP );
#  endif


// verify the input parameters
   Aux_Check_Parameter();


// initialize the timer function
#  ifdef TIMING
   Aux_CreateTimer();
#  endif


// load the tables of the flag criteria from the input files "Input__Flag_XXX"
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
      case PAR_INIT_BY_FUNCTION:
         if ( Par_Init_ByFunction_Ptr == NULL )    Aux_Error( ERROR_INFO, "Par_Init_ByFunction_Ptr == NULL !!\n" );
         Par_Init_ByFunction_Ptr( amr->Par->NPar_Active, amr->Par->NPar_Active_AllRank,
                                  amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                  amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ, amr->Par->Time,
                                  amr->Par->Attribute );
         break;

      case PAR_INIT_BY_RESTART:
         break;   // nothing to do here for the restart mode

      case PAR_INIT_BY_FILE:
         Par_Init_ByFile();
         break;

      default:
         Aux_Error( ERROR_INFO, "unsupported particle initialization (%s = %d) !!\n",
                    "PAR_INIT", (int)amr->Par->Init );
   }

   if ( amr->Par->Init != PAR_INIT_BY_RESTART )    Par_Aux_InitCheck();
#  endif // #ifdef PARTICLE


// initialize the AMR structure and fluid field
   switch ( OPT__INIT )
   {
      case INIT_BY_FUNCTION :    Init_ByFunction();   break;

      case INIT_BY_RESTART :     Init_ByRestart();    break;

      case INIT_BY_FILE :        Init_ByFile();       break;

      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__INIT", OPT__INIT );
   }


// user-defined initialization
   if ( Init_User_Ptr != NULL )  Init_User_Ptr();


// record the initial weighted load-imbalance factor
#  ifdef LOAD_BALANCE
   if ( OPT__RECORD_LOAD_BALANCE )  LB_EstimateLoadImbalance();
#  endif


#  ifdef GRAVITY
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
//    initialize the k-space Green's function for the isolated BC.
      if ( OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();


//    evaluate the initial average density if it is not set yet (may already be set in Init_ByRestart)
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

} // FUNCTION : Init_GAMER
