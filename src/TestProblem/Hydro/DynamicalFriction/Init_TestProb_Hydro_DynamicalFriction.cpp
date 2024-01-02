#include "GAMER.h"
#include "TestProb.h"
#include "Par_EquilibriumIC.h"
#include "string"
#include <stdio.h>
#include "global_var.h"

using namespace std;

// negligibly small uniform density and energy
double GC_SmallGas;

// parameter for center setting
double FixCenter;
double SearchRadius;


// GC position
double GC_xx = 100000000.0;
double GC_yy = 100000000.0;
double GC_zz = 100000000.0;

// declare the potential minimum last step
double min_pot_last[3] ;


// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_GC( const long NPar_ThisRank, const long NPar_AllRank,
                                   real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                   real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );
#endif

// external potential routines
void Init_ExtPot_GC();

void Mis_UserWorkBeforeNextLevel_Find_GC( const int lv, const double TimeNew, const double TimeOld, const double dt );

void User_Output();
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Aux_Record_User_GC();


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


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
     Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  else
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  10;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// load run-time parameters
   const char* FileName = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "GC_SmallGas",             &GC_SmallGas,           1e-10,          0.,               NoMax_double      );

   ReadPara->Add( "GC_POSX",                 &GC_xx,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_POSY",                 &GC_yy,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_POSZ",                 &GC_zz,               NoDef_double,  NoMin_double,     NoMax_double      );
   
   ReadPara->Add( "FIX_CENTER",              &FixCenter,           NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "SEARCH_RADIUS",           &SearchRadius,        NoDef_double,  NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );
   if ( MPI_Rank == 0)
   {
   Aux_Message(stdout, "Setparameter(): %7.5f %7.5f %7.5f \n",GC_xx,GC_yy,GC_zz);
   Aux_Message(stdout, "Setparameter(): %7.1f %7.1f \n",FixCenter,SearchRadius);
   }
   delete ReadPara;

} // FUNCTION : SetParameter


//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_User_GC
// Description :  Update the GC position and Halo position in each step
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_User_GC()
{

// set the center finding method base on Input__TestProb
if (FixCenter == 1.0){

if ( MPI_Rank == 0)
{
	char Filename[MAX_STRING];
	sprintf( Filename, "%s", "Record__RESULT_Center.txt" );
	FILE *File = fopen( Filename, "a" );
        if (Time[0]==0.0){
        fprintf(File, "%15s%15s%15s%15s\n", "Time", "CenterX", "CenterY", "CenterZ");
	}
	fprintf(File, "%13.7f\t%13.7e\t%13.7e\t%13.7e\n",
           Time[0]  ,amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2]);
	fclose ( File );	
}
}else{
// 1. Find the minimum potential position
// Extrema.Radius    = HUGE_NUMBER; // entire domain
Extrema_t Extrema;
Extrema.Field     = _POTE;
Extrema.Radius = SearchRadius*amr->dh[MAX_LEVEL]; //the cell width of the finest level of resolution in the AMR grid multiply by the searchRadius setting

if (Time[0]==0.0)
{
  Extrema.Center[0] = amr->BoxCenter[0];
  Extrema.Center[1] = amr->BoxCenter[1];
  Extrema.Center[2] = amr->BoxCenter[2];
  min_pot_last[0] = amr->BoxCenter[0];
  min_pot_last[1] = amr->BoxCenter[1];
  min_pot_last[2] = amr->BoxCenter[2];
}
else{
Extrema.Center[0] = min_pot_last[0];
Extrema.Center[1] = min_pot_last[1];
Extrema.Center[2] = min_pot_last[2];
}
Aux_FindExtrema( &Extrema, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );

min_pot_last[0] = Extrema.Coord[0];
min_pot_last[1] = Extrema.Coord[1];
min_pot_last[2] = Extrema.Coord[2];

// 2. write them into a .txt file

if ( MPI_Rank == 0)
{
	char Filename[MAX_STRING];
	sprintf( Filename, "%s", "Record__RESULT_Center.txt" );
	FILE *File = fopen( Filename, "a" );
        if (Time[0]==0.0){
        fprintf(File, "%15s%15s%15s%15s\n", "Time", "CenterX", "CenterY", "CenterZ");
	}
	fprintf(File, "%13.7f\t%13.7e\t%13.7e\t%13.7e\n",
           Time[0]  ,Extrema.Coord[0], Extrema.Coord[1], Extrema.Coord[2]);
	fclose ( File );
	
	char Filename_[MAX_STRING];
	sprintf( Filename_, "%s", "Record__Previous_Center.txt" );
	FILE *File_ = fopen( Filename_, "a" );
        if (Time[0]==0.0){
        fprintf(File_, "%15s%15s%15s%15s%15s\n", "Time", "CenterX", "CenterY", "CenterZ","Search R");
	}
	fprintf(File_, "%13.7f\t%13.7e\t%13.7e\t%13.7e\t%13.7e\n",
           Time[0]  ,min_pot_last[0], min_pot_last[1], min_pot_last[2],Extrema.Radius);
	fclose ( File_ );
}
};
// 2. Find the GC's position
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!! This routine assume there is only one GC particle !!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double GC_x=99999.9;
double GC_y=99999.9;
double GC_z=99999.9;

for (long p=0; p<amr->Par->NPar_AcPlusInac; p++) {
if ( amr->Par->Mass[p] < PTYPE_GC )  continue;
   if ( amr->Par->Type[p] == PTYPE_GC) {
      GC_x = amr->Par->PosX[p]; 
      GC_y = amr->Par->PosY[p]; 
      GC_z = amr->Par->PosZ[p];
      
      char Filename[MAX_STRING];
      sprintf( Filename, "%s", "Record__RESULT_GC_position.txt" );
      FILE *File = fopen( Filename, "a" );
      if (Time[0]==0.0){
        fprintf(File, "%15s%15s%15s%15s\n", "Time", "GCX", "GCY", "GCZ");
      }
      fprintf(File, "%13.7f\t%13.7e\t%13.7e\t%13.7e\n", Time[0],GC_x,GC_y,GC_z );
      fclose ( File );	
   }
}

} // FUNCTION : Aux_Record_User_GC

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

   fluid[DENS] = GC_SmallGas;
   fluid[MOMX] = 0;
   fluid[MOMY] = 0;
   fluid[MOMZ] = 0;
#  ifdef GRAVITY
   fluid[ENGY] = GC_SmallGas;
#  endif

// just set all passive scalars as zero
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_DynamicalFriction
// Description :  Test problem initializer for dynamical friction problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_DynamicalFriction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();
   Init_Function_User_Ptr  = SetGridIC;
   Aux_Record_User_Ptr     = Aux_Record_User_GC;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_GC;
#  endif
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr         = Init_ExtPot_GC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GC
