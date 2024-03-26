#include "GAMER.h"
#include "TestProb.h"
#include "Par_EquilibriumIC.h"
#include "string"
#include <stdio.h>

using namespace std;

// negligibly small uniform density and energy
double GC_SmallGas;

// parameter for center setting
bool FixCenter;
double SearchRadius;

extern double GC_MASS;
extern double GC_R;
extern bool PURE_TABLE;
extern double SEARCH_RADIUS;
extern double Halo_Rho0;
extern double Halo_Rs;
extern double Halo_Rt;

extern double Halo_Profile_Param_a;
extern double Halo_Profile_Param_b;
extern double Halo_Profile_Param_c;


static char HaloType[MAX_STRING];
static char TableName[MAX_STRING];
extern const char* HaloType_g = nullptr;
extern const char* TableName_g = nullptr;
// declare the potential minimum last step
double min_pot_last[3] ;


// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_DynamicalFriction( const long NPar_ThisRank, const long NPar_AllRank,
                                   real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                   real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );


#endif

// testing the subtraction
// void Output_UserWorkBeforeOutput_DF();


void Mis_UserWorkBeforeNextLevel_Find_GC( const int lv, const double TimeNew, const double TimeOld, const double dt );

void AdjustGCPotential(const bool add, const double GC_x, const double GC_y, const double GC_z);

void User_Output();
// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Aux_Record_User_DynamicalFriction();

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
   const double End_T_Default    =  14000;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   // load run-time parameters
   const char* FileName = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "GC_SmallGas",             &GC_SmallGas,           1e-10,         0.,               NoMax_double      );

   ReadPara->Add( "GC_initial_R",            &GC_R,                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "GC_MASS",                 &GC_MASS,               NoDef_double,  NoMin_double,     NoMax_double      );

   ReadPara->Add( "SEARCH_RADIUS",           &SearchRadius,          NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "FIX_CENTER",              &FixCenter,             Useless_bool,  Useless_bool,     Useless_bool      );
   ReadPara->Add( "PURE_TABLE",              &PURE_TABLE,            Useless_bool,  Useless_bool,     Useless_bool      );
   ReadPara->Add( "Density_Table_Name",      TableName,              "None",        Useless_str,      Useless_str       );

   ReadPara->Add( "HALO_TYPE",               HaloType,               "None",        Useless_str,      Useless_str       );
   ReadPara->Add( "HALO_RHO_0",              &Halo_Rho0,             NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "HALO_Rs",                 &Halo_Rs,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "HALO_Rt",                 &Halo_Rt,               NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Halo_Profile_Param_a",    &Halo_Profile_Param_a,  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Halo_Profile_Param_b",    &Halo_Profile_Param_b,  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Halo_Profile_Param_c",    &Halo_Profile_Param_c,  NoDef_double,  NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );
   if ( MPI_Rank == 0)
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  GC small Gas        = %13.5e\n", GC_SmallGas );
      Aux_Message( stdout, "  GC initial radius   = %13.5e\n", GC_R );
      Aux_Message( stdout, "  GC Mass             = %13.5e\n", GC_MASS );
      Aux_Message( stdout, "  Fix Center          = %d\n",     FixCenter );
      Aux_Message( stdout, "  Search Radius       = %13.5e\n",  SearchRadius );
      Aux_Message( stdout, "  Pure Table          = %d\n",     PURE_TABLE );
      Aux_Message( stdout, "  Density Table Name  = %s\n",     TableName );
      Aux_Message( stdout, "  Halo Type           = %s\n",     HaloType );
      Aux_Message( stdout, "  HALO_RHO_0          = %13.5e\n", Halo_Rho0 );
      Aux_Message( stdout, "  HALO_Rs             = %13.5e\n", Halo_Rs );
      Aux_Message( stdout, "  HALO_Rt             = %13.5e\n", Halo_Rt );
      Aux_Message( stdout, "  Halo_Profle_Param_a = %13.5f\n", Halo_Profile_Param_a );
      Aux_Message( stdout, "  Halo_Profle_Param_b = %13.5f\n", Halo_Profile_Param_b );
      Aux_Message( stdout, "  Halo_Profle_Param_c = %13.5f\n", Halo_Profile_Param_c );
      Aux_Message( stdout, "=============================================================================\n" );
   }
   delete ReadPara;
   // Make the string be global parameters
   HaloType_g = HaloType;
   TableName_g = TableName;

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  AdjustGCPotential()
// Description :  Adjust the Potential
//
// Note        :  
// Parameter   :  add : the booldean value that determines the action (subtraction or add back)
//                GC_x,GC_y,GC_z : the current GC's position
//-------------------------------------------------------------------------------------------------------
void AdjustGCPotential(const bool add, const double GC_x, const double GC_y, const double GC_z)
{

   double sign = add ? 1.0 : -1.0;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const double dh = amr->dh[lv];
#  pragma omp for schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         const double *EdgeL              = amr->patch[0][lv][PID]->EdgeL;
         real   (*PotPtr)[PS1][PS1] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
         const double x0                  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0                  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0                  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;
         for (int k=0; k<PS1; k++)
         {
            const double z  = z0 + k*dh;
            const double dz = z - GC_z;
            for (int j=0; j<PS1; j++)
            {
               const double y  = y0 + j*dh;
               const double dy = y - GC_y;
               for (int i=0; i<PS1; i++)
               {
                  const double x  = x0 + i*dh;
                  const double dx = x - GC_x;
                  //const double r      = sqrt( SQR(dx) + SQR(dy) + SQR(dz) + SQR(dh) );
                  const double r      = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
                  const double GC_pot = -NEWTON_G*GC_MASS/r;
         
                  PotPtr[k][j][i] += sign * (real)GC_pot;
               }
            }
          } 
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)
}

//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_User_DynamicalFriction
// Description :  Update the GC position and Halo position in each step
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//                4. Notice that this function can only get "ONE" GC, multiple GCs are not supported.
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_User_DynamicalFriction()
{

// A. Find the GC's position
   double GC_x_local = 0.0, GC_y_local = 0.0, GC_z_local = 0.0;
   double GC_x, GC_y, GC_z;



   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++) {
   if ( amr->Par->Mass[p] < (real)0.0 )  continue;
      if ( amr->Par->Type[p] == PTYPE_GC)
      {
         GC_x_local = amr->Par->PosX[p];
         GC_y_local = amr->Par->PosY[p];
         GC_z_local = amr->Par->PosZ[p];

         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_GC_position" );
         FILE *File = fopen( Filename, "a" );
         if ( Time[0]==0.0 )
         {
           fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " GCX", " GCY", " GCZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n", Time[0],GC_x_local,GC_y_local,GC_z_local );
         fclose ( File );
      }
   }

// A. Broadcast the GC's position to all rank
   MPI_Allreduce(&GC_x_local, &GC_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&GC_y_local, &GC_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&GC_z_local, &GC_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   
   MPI_Bcast(&GC_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(&GC_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
   MPI_Barrier(MPI_COMM_WORLD);

//   Aux_Message(stdout, "[%02d] GC's position is %15.6e %15.6e %15.6e\n",MPI_Rank,GC_x,GC_y,GC_z);




   // B. Set the center finding method base on Input__TestProb
   if ( FixCenter )
   {
      if ( MPI_Rank == 0 )
         {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_Center" );
         FILE *File = fopen( Filename, "a" );
         if ( Time[0]==0.0 )
         {
            fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " CenterX", " CenterY", " CenterZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n",Time[0]  ,amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2]);
         fclose ( File );
      }
   }
   else
   {


   // B-1. subtract GC potential
      AdjustGCPotential(false,GC_x,GC_y,GC_z);     
      
       
   // B-2. Find the minimum potential position
      Extrema_t Extrema;
      Extrema.Field     = _POTE;
      Extrema.Radius    = SearchRadius*amr->dh[MAX_LEVEL];

      if ( Time[0]==0.0 )
      {
         Extrema.Center[0] = amr->BoxCenter[0];
         Extrema.Center[1] = amr->BoxCenter[1];
         Extrema.Center[2] = amr->BoxCenter[2];
         min_pot_last[0] = amr->BoxCenter[0];
         min_pot_last[1] = amr->BoxCenter[1];
         min_pot_last[2] = amr->BoxCenter[2];
      }
      else
      {
         Extrema.Center[0] = min_pot_last[0];
         Extrema.Center[1] = min_pot_last[1];
         Extrema.Center[2] = min_pot_last[2];
      }

      Aux_FindExtrema( &Extrema, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );

      min_pot_last[0] = Extrema.Coord[0];
      min_pot_last[1] = Extrema.Coord[1];
      min_pot_last[2] = Extrema.Coord[2];

   // B-3. write them into a file

      if ( MPI_Rank == 0 )
      {
         char Filename[MAX_STRING];
         sprintf( Filename, "%s", "Record__Result_Center" );
         FILE *File = fopen( Filename, "a" );
         if (Time[0]==0.0)
         {
            fprintf(File, "%15s\t%15s\t%15s\t%15s\n", "#          Time", " CenterX", " CenterY", " CenterZ");
         }
         fprintf(File, "%15.7f\t%15.7e\t%15.7e\t%15.7e\n",
            Time[0]  ,Extrema.Coord[0], Extrema.Coord[1], Extrema.Coord[2]);
         fclose ( File );
      }

    // B-4 Add back the GC's potential
      AdjustGCPotential(true,GC_x,GC_y,GC_z);     
   };

} // FUNCTION : Aux_Record_User_DynamicalFriction



//void Output_UserWorkBeforeOutput_DF(){
//// A. Find the GC's position
//
//   double GC_x_local = 0.0, GC_y_local = 0.0, GC_z_local = 0.0;
//   double GC_x, GC_y, GC_z;
//
//
//
//   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++) {
//   if ( amr->Par->Mass[p] < (real)0.0 )  continue;
//      if ( amr->Par->Type[p] == PTYPE_GC)
//      {
//         GC_x_local = amr->Par->PosX[p];
//         GC_y_local = amr->Par->PosY[p];
//         GC_z_local = amr->Par->PosZ[p];
//
//      }
//   }
//
//// A. Broadcast the GC's position to all rank
//   MPI_Allreduce(&GC_x_local, &GC_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//   MPI_Allreduce(&GC_y_local, &GC_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//   MPI_Allreduce(&GC_z_local, &GC_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
//   
//   MPI_Bcast(&GC_x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//   MPI_Bcast(&GC_y,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
//   MPI_Bcast(&GC_z,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   
//   MPI_Barrier(MPI_COMM_WORLD);
//
////   Aux_Message(stdout, "[%02d] GC's position is %15.6e %15.6e %15.6e\n",MPI_Rank,GC_x,GC_y,GC_z);
//
//
//
//   AdjustGCPotential(false,GC_x,GC_y,GC_z);     
//   if (MPI_Rank == 0)
//   {
//      Aux_Message(stdout,"Successfully subtracted the GC's potential!");
//   }
//
//}














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
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = GC_SmallGas;

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
//   Output_UserWorkBeforeOutput_Ptr = Output_UserWorkBeforeOutput_DF;
   Aux_Record_User_Ptr     = Aux_Record_User_DynamicalFriction;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_DynamicalFriction;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GC
