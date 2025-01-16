#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize GAMER
//
// Note        :  1. Function pointers "Init_User_Ptr" and "Init_User_AfterPoisson_Ptr" may be set by a
//                   test problem initializer
//
// Parameter   :  argc, argv: Command line arguments
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


// reset parameters
// --> must be called after Init_Unit()
   Init_ResetParameter();


// load the tables of the flag criteria from the input files "Input__Flag_XXX"
   Init_Load_FlagCriteria();


// load the dump table from the input file "Input__DumpTable"
   if ( OPT__OUTPUT_MODE == OUTPUT_USE_TABLE )
#  ifdef PARTICLE
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS || OPT__OUTPUT_PAR_MODE )
#  else
   if ( OPT__OUTPUT_TOTAL || OPT__OUTPUT_PART || OPT__OUTPUT_USER || OPT__OUTPUT_BASEPS )
#  endif
   Init_Load_DumpTable();


// initialize OpenMP settings
#  ifdef OPENMP
   Init_OpenMP();
#  endif


// initialize GPU
// --> must be called before Init_ExtAccPot() and EoS_Init()
#  ifdef GPU
   CUAPI_SetDevice( OPT__GPUID_SELECT );

   CUAPI_SetCache();
#  endif // #ifdef GPU


// initialize yt inline analysis
#  ifdef SUPPORT_LIBYT
   YT_Init( *argc, *argv );
#  endif


// initialize Grackle
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_ACTIVATE )    Grackle_Init();
#  endif


#  ifdef SUPPORT_FFTW
// initialize FFTW
   Init_FFTW();
#  endif


// initialize the test problem parameters
   Init_TestProb();


// initialize parameters for the parallelization
// --> call it after Init_TestProb() since Init_TestProb() may overwrite the total number of particles
   Init_Parallelization();


// initialize all fields and particle attributes
// --> Init_Field() must be called before CUAPI_SetConstMemory()
   Init_Field();
#  ifdef PARTICLE
   Par_Init_Attribute();
#  endif


// initialize the external potential and acceleration parameters
// --> must be called after Init_TestProb()
#  ifdef GRAVITY
   Init_ExtAccPot();
#  endif


// initialize the EoS routines
#  if ( MODEL == HYDRO )
   EoS_Init();
#  endif


// initialize the source-term routines
// --> must be called before memory allocation
   Src_Init();


// initialize the feedback routines
#  ifdef FEEDBACK
   FB_Init();
#  endif


// initialize the microphysics
   Microphysics_Init();


// initialize the user-defined derived fields
   if ( OPT__OUTPUT_USER_FIELD )
   {
      if ( Init_DerivedField_User_Ptr != NULL )
         Init_DerivedField_User_Ptr();

      else
         Aux_Error( ERROR_INFO, "Init_DerivedField_User_Ptr == NULL for OPT__OUTPUT_USER_FIELD !!\n" );
   }


// set GPU constant memory
// --> must be called after Init_Field() and Init_ExtAccPot()
#  ifdef GPU
   CUAPI_SetConstMemory();
#  endif


// verify the input parameters
   Aux_Check_Parameter();


// initialize the timer function
#  ifdef TIMING
   Aux_CreateTimer();
#  endif


// initialize memory pool
   if ( OPT__MEMORY_POOL )    Init_MemoryPool();


// allocate memory for several CPU/GPU global arrays
#  ifdef GPU
   CUAPI_MemAllocate();
#  endif

   Init_MemAllocate();


// load the external potential table
// --> before Init_ByFunction() so that the test problem initializer can access
//     the external potential table if required
// --> after Init_MemAllocate() to allocate the potential table array first
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_TABLE )   Init_LoadExtPotTable();
#  endif


// initialize variables for dumping particle attributes mapped from mesh quantities
#  ifdef PARTICLE
   if ( OPT__OUTPUT_PAR_MESH )   Par_Init_Attribute_Mesh();
#  endif


// initialize particles
#  ifdef PARTICLE
   switch ( amr->Par->Init )
   {
      case PAR_INIT_BY_FUNCTION:
         if ( Par_Init_ByFunction_Ptr != NULL )
            Par_Init_ByFunction_Ptr( amr->Par->NPar_Active, amr->Par->NPar_Active_AllRank,
                                     amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                     amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ, amr->Par->Time,
                                     amr->Par->Type, amr->Par->AttributeFlt, amr->Par->AttributeInt );
         else
            Aux_Error( ERROR_INFO, "Par_Init_ByFunction_Ptr == NULL for PAR_INIT = 1 !!\n" );
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


// initialize the local AMR structure and grid fields
   switch ( OPT__INIT )
   {
      case INIT_BY_FUNCTION :    Init_ByFunction();   break;

      case INIT_BY_RESTART :     Init_ByRestart();    break;

      case INIT_BY_FILE :        Init_ByFile();       break;

      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "OPT__INIT", OPT__INIT );
   }


// initialize the global AMR structure if required
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   delete GlobalTree;   // in case it has been allocated already
   GlobalTree = new LB_GlobalTree;
#  endif


// ensure B field consistency on the shared interfaces between sibling patches
#  if ( MODEL == HYDRO  &&  defined MHD )
   if ( OPT__SAME_INTERFACE_B )
   for (int lv=0; lv<NLEVEL; lv++)  MHD_SameInterfaceB( lv );
#  endif


// user-defined initialization (before the Poisson solver)
   if ( Init_User_Ptr != NULL )  Init_User_Ptr();


// record the initial weighted load-imbalance factor
#  ifdef LOAD_BALANCE
   if ( OPT__RECORD_LOAD_BALANCE )  LB_EstimateLoadImbalance();
#  endif


#  ifdef GRAVITY
   if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
   {
#     ifdef SUPPORT_FFTW
//    initialize the k-space Green's function for the isolated BC.
      if ( OPT__SELF_GRAVITY  &&  OPT__BC_POT == BC_POT_ISOLATED )    Init_GreenFuncK();
#     endif


//    evaluate the initial average density if it is not set yet (may already be set in Init_ByRestart)
      if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();


//    evaluate the gravitational potential
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

         Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Rho_ParaBuf, USELB_YES );

         Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false, true );

         if ( lv > 0 )
         Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=0; lv<NLEVEL; lv++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" );
   } // if ( OPT__SELF_GRAVITY  ||  OPT__EXT_POT )
#  endif // #ifdef GARVITY


#  ifdef PARTICLE

// initialize particle acceleration
#  if ( defined MASSIVE_PARTICLES  &&  defined STORE_PAR_ACC )
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" );

   const bool StoreAcc_Yes    = true;
   const bool UseStoredAcc_No = false;

   for (int lv=0; lv<NLEVEL; lv++)
      Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" );
#  endif // #if ( defined MASSIVE_PARTICLES  &&  defined STORE_PAR_ACC )

#  ifdef TRACER
// initialize tracer particles
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Initializing tracer particles" );

   const bool MapOnly_Yes = true;

   for (int lv=0; lv<NLEVEL; lv++)
      Par_UpdateTracerParticle( lv, Time[lv], NULL_REAL, MapOnly_Yes );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Initializing tracer particles" );
#  endif // #ifdef TRACER

#  endif // #ifdef PARTICLE


// user-defined initialization (after the Poisson solver)
   if ( Init_User_AfterPoisson_Ptr != NULL )    Init_User_AfterPoisson_Ptr();


// initialize source-term fields (e.g., cooling time)
// --> necessary for, for example, estimating the source-term time-step at the first step
// --> must ensure each source term does not modify any fluid field (e.g., gas internal energy) when dt=0.0
//     --> source terms violating this criterion (e.g., deleptonization) must be temporarily disabled before calling Src_AdvanceDt()
   if ( OPT__INIT != INIT_BY_RESTART )
   {
      const bool OverlapMPI_No   = false;
      const bool Overlap_Sync_No = false;

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Initializing source-term fields" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

//       adopt dt=0.0 to prevent any update of the fluid fields
         Src_AdvanceDt( lv, Time[lv], Time[lv], 0.0, amr->FluSg[lv], amr->MagSg[lv], OverlapMPI_No, Overlap_Sync_No );

//       must exchange all source-term fields since they are currently implemeted as passive scalars,
//       which will be advected by the fluid solver (and then be overwritten by Src_AdvanceDt())
//       --> here we exchange all _TOTAL and _MAG fields just for simplicity
         Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=0; lv<NLEVEL; lv++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Initializing source-term fields" );
   } // if ( OPT__INIT != INIT_BY_RESTART )

} // FUNCTION : Init_GAMER
