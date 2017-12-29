#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
       bool    Merger_Coll;               // (true/false) --> test (cluster merger / single cluster)
static char    Merger_File_Prof1[1000];   // profile table of cluster 1
static char    Merger_File_Prof2[1000];   // profile table of cluster 2
       char    Merger_File_Par1 [1000];   // particle file of cluster 1
       char    Merger_File_Par2 [1000];   // particle file of cluster 2
       double  Merger_Coll_D;             // initial distance between two clusters
       double  Merger_Coll_B;             // impact parameter
       double  Merger_Coll_BulkVel1;      // bulk velocity of cluster 1 (on the left  side)
       double  Merger_Coll_BulkVel2;      // bulk velocity of cluster 2 (on the right side)

static double *Merger_Prof1 = NULL;       // radial profiles [gas mass density/gas pressure/radius] of cluster 1
static double *Merger_Prof2 = NULL;       // radial profiles [gas mass density/gas pressure/radius] of cluster 2
static int     Merger_NBin1;              // number of radial bins of cluster 1
static int     Merger_NBin2;              // number of radial bins of cluster 2
// =======================================================================================


// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_Merger( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                 real *ParPassive[PAR_NPASSIVE] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC for this test !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined PARTICLE )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format (some handy constants are defined in TestProb.h):
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",        &VARIABLE,               DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************
   ReadPara->Add( "Merger_Coll",            &Merger_Coll,             true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_File_Prof1",       Merger_File_Prof1,       Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par1",        Merger_File_Par1,        Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Prof2",       Merger_File_Prof2,       Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par2",        Merger_File_Par2,        Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_Coll_D",          &Merger_Coll_D,          -1.0,             0.0,           NoMax_double   );
   ReadPara->Add( "Merger_Coll_B",          &Merger_Coll_B,          -1.0,             0.0,           NoMax_double   );
   ReadPara->Add( "Merger_Coll_BulkVel1",   &Merger_Coll_BulkVel1,   -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_BulkVel2",   &Merger_Coll_BulkVel2,   -1.0,             NoMin_double,  NoMax_double   );

   ReadPara->Read( FileName );

   delete ReadPara;

// convert to code units
   if ( Merger_Coll )
   {
      Merger_Coll_D        *= Const_kpc          / UNIT_L;
      Merger_Coll_B        *= Const_kpc          / UNIT_L;
      Merger_Coll_BulkVel1 *= (Const_km/Const_s) / UNIT_V;
      Merger_Coll_BulkVel2 *= (Const_km/Const_s) / UNIT_V;
   }


// (2) load the radial profiles
   if ( OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
      const int  NCol         = 3;           // total number of columns to load
      const int  Col[NCol]    = {2, 6, 7};   // target columns: (density, pressure, radius)
      double *Table_D, *Table_P, *Table_R;

//    cluster 1
      Merger_NBin1 = Aux_LoadTable( Merger_Prof1, Merger_File_Prof1, NCol, Col, RowMajor_No, AllocMem_Yes );

//    convert to code units (assuming the input units are cgs)
      Table_D = Merger_Prof1 + 0*Merger_NBin1;
      Table_P = Merger_Prof1 + 1*Merger_NBin1;
      Table_R = Merger_Prof1 + 2*Merger_NBin1;

      for (int b=0; b<Merger_NBin1; b++)
      {
         Table_D[b] /= UNIT_D;
         Table_P[b] /= UNIT_P;
         Table_R[b] /= UNIT_L;
      }

//    cluster 2
      if ( Merger_Coll ) {
      Merger_NBin2 = Aux_LoadTable( Merger_Prof2, Merger_File_Prof2, NCol, Col, RowMajor_No, AllocMem_Yes );

//    convert to code units (assuming the input units are cgs)
      Table_D = Merger_Prof2 + 0*Merger_NBin2;
      Table_P = Merger_Prof2 + 1*Merger_NBin2;
      Table_R = Merger_Prof2 + 2*Merger_NBin2;

      for (int b=0; b<Merger_NBin2; b++)
      {
         Table_D[b] /= UNIT_D;
         Table_P[b] /= UNIT_P;
         Table_R[b] /= UNIT_L;
      }
      } // if ( Merger_Coll )
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID   = %d\n",           TESTPROB_ID );
      Aux_Message( stdout, "   profile file 1   = %s\n",           Merger_File_Prof1 );
      if ( Merger_Coll )
      Aux_Message( stdout, "   profile file 2   = %s\n",           Merger_File_Prof2 );
      Aux_Message( stdout, "   particle file 1  = %s\n",           Merger_File_Par1 );
      if ( Merger_Coll )
      Aux_Message( stdout, "   particle file 2  = %s\n",           Merger_File_Par2 );
      Aux_Message( stdout, "   test mode        = %s\n",          (Merger_Coll)? "merging cluster":"single cluster" );
      if ( Merger_Coll ) {
      Aux_Message( stdout, "   initial distance = % 14.7e kpc\n",  Merger_Coll_D*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   impact parameter = % 14.7e kpc\n",  Merger_Coll_B*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   bulk velocity 1  = %+14.7e km/s\n", Merger_Coll_BulkVel1*UNIT_V/(Const_km/Const_s) );
      Aux_Message( stdout, "   bulk velocity 2  = %+14.7e km/s\n", Merger_Coll_BulkVel2*UNIT_V/(Const_km/Const_s) ); }
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const double  BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double *Table_D1     = Merger_Prof1 + 0*Merger_NBin1;    // density  table of cluster 1
   const double *Table_P1     = Merger_Prof1 + 1*Merger_NBin1;    // pressure table of cluster 1
   const double *Table_R1     = Merger_Prof1 + 2*Merger_NBin1;    // radius   table of cluster 1
   const double *Table_D2     = Merger_Prof2 + 0*Merger_NBin2;    // density  table of cluster 2
   const double *Table_P2     = Merger_Prof2 + 1*Merger_NBin2;    // pressure table of cluster 2
   const double *Table_R2     = Merger_Prof2 + 2*Merger_NBin2;    // radius   table of cluster 2

   if ( Merger_Coll )
   {
      const double ClusterCenter1[3] = { BoxCenter[0]-0.5*Merger_Coll_D, BoxCenter[1]-0.5*Merger_Coll_B, BoxCenter[2] };
      const double ClusterCenter2[3] = { BoxCenter[0]+0.5*Merger_Coll_D, BoxCenter[1]+0.5*Merger_Coll_B, BoxCenter[2] };

      double r1, r2, Dens1, Dens2, Pres1, Pres2, Vel;

//    for each cell, we sum up the density and pressure from each halo and then calculate the weighted velocity
      r1    = sqrt( SQR(x-ClusterCenter1[0]) + SQR(y-ClusterCenter1[1]) + SQR(z-ClusterCenter1[2]) );
      r2    = sqrt( SQR(x-ClusterCenter2[0]) + SQR(y-ClusterCenter2[1]) + SQR(z-ClusterCenter2[2]) );
      Dens1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_D1, r1 );
      Dens2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_D2, r2 );
      Pres1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_P1, r1 );
      Pres2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_P2, r2 );
      Vel   = ( Merger_Coll_BulkVel1*Dens1 + Merger_Coll_BulkVel2*Dens2 ) / ( Dens1 + Dens2 );

      if ( Dens1 == NULL_REAL  ||  Pres1 == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e for cluster 1 (probably outside the input table) !!\n", r1 );

      if ( Dens2 == NULL_REAL  ||  Pres2 == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e for cluster 2 (probably outside the input table) !!\n", r2 );

      fluid[DENS] = Dens1 + Dens2;
      fluid[MOMX] = fluid[DENS]*Vel;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = ( Pres1 + Pres2 ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   }

   else
   {
      double r, Dens, Pres;

      r    = sqrt( SQR(x-BoxCenter[0]) + SQR(y-BoxCenter[1]) + SQR(z-BoxCenter[2]) );
      Dens = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_D1, r );
      Pres = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_P1, r );

      if ( Dens == NULL_REAL  ||  Pres == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (probably outside the input table)!!\n", r );

      fluid[DENS] = Dens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Pres / ( GAMMA - 1.0 );
   } // if ( Merger_Coll ) ... else ...

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_ClusterMerger
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_ClusterMerger()
{

   delete [] Merger_Prof1;
   delete [] Merger_Prof2;

} // FUNCTION : End_ClusterMerger
#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ClusterMerger_vs_Flash
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ClusterMerger_vs_Flash()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined PARTICLE )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = End_ClusterMerger;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
   Par_Init_ByFunction_Ptr  = Par_Init_ByFunction_Merger;
#  endif // if ( MODEL == HYDRO  &&  defined PARTICLE )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ClusterMerger_vs_Flash
