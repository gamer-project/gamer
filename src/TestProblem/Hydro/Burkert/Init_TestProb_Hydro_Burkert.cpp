#include "GAMER.h"
#include "TestProb.h"
#include"math.h"
#include <gsl/gsl_integration.h>

/***Deep Copy***/
// problem-specific global variables
// =======================================================================================
       int    Burkert_RSeed;        // random seed for setting particle position and velocity
       double Burkert_Rho0;         // peak density
       double Burkert_R0;           // scale radius
       double Burkert_MaxR;         // maximum radius for particles
       bool   Burkert_Collision;    // true/false --> two colliding Burkert clouds/single Burkert cloud
       double Burkert_Collision_D;  // distance between two colliding Burkert clouds
       double Burkert_Center[3];    // central coordinates
       double Burkert_BulkVel[3];   // bulk velocity
       double Burkert_GasMFrac;     // gas mass fraction
       int    Burkert_MassProfNBin; // number of radial bins in the mass profile table
static bool   Burkert_AddColor;     // assign different colors to different clouds for Burkert_Collision

static double Burkert_FreeT;        // free-fall time at Burkert_R0

static FieldIdx_t Burkert_Idx_Cloud0 = Idx_Undefined;    // field indices for Burkert_AddColor
static FieldIdx_t Burkert_Idx_Cloud1 = Idx_Undefined;
// =======================================================================================
/***Deep Copy***/
/***gsl function***/

double de_pressure (double x, void * parameters) {
   
  double dens=Burkert_Rho0 /((1+x)*(1+x*x));
  double g=4*M_PI*NEWTON_G*Burkert_R0*Burkert_Rho0*(-0.5*atan(x)+0.25*log(x*x+1)+0.5*log(x+1))/(x*x);
  double de_pressure = dens*g;
  return de_pressure;
}
/***gsl function***/
// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_Burkert( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *AllAttribute[PAR_NATT_TOTAL] );
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
/***Deep Copy***/
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
   } // if ( MPI_Rank == 0 )x


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );
/***Deep Copy***/
} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
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
/***Deep Copy***/
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Burkert_RSeed",        &Burkert_RSeed,         123,           0,                NoMax_int         );
   ReadPara->Add( "Burkert_Rho0",         &Burkert_Rho0,          1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Burkert_R0",           &Burkert_R0,            0.1,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Burkert_MaxR",         &Burkert_MaxR,          0.375,         Eps_double,       NoMax_double      );
   ReadPara->Add( "Burkert_Collision",    &Burkert_Collision,     false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Burkert_Collision_D",  &Burkert_Collision_D,   1.5,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_CenterX",      &Burkert_Center[0],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_CenterY",      &Burkert_Center[1],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_CenterZ",      &Burkert_Center[2],     NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_BulkVelX",     &Burkert_BulkVel[0],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_BulkVelY",     &Burkert_BulkVel[1],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_BulkVelZ",     &Burkert_BulkVel[2],    0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Burkert_GasMFrac",     &Burkert_GasMFrac,      0.5,           Eps_double,       1.0               );
   ReadPara->Add( "Burkert_MassProfNBin", &Burkert_MassProfNBin,  1000,          2,                NoMax_int         );
   ReadPara->Add( "Burkert_AddColor",     &Burkert_AddColor,      false,         Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   for (int d=0; d<3; d++)
      if ( Burkert_Center[d] == NoDef_double )  Burkert_Center[d] = 0.5*amr->BoxSize[d];

   if ( !Burkert_Collision  &&  Burkert_AddColor )
      Aux_Error( ERROR_INFO, "\"Burkert_AddColor\" must work with \"Burkert_Collision\" !!\n" );

// (1-3) check and reset the runtime parameters
   if ( Burkert_AddColor  &&  NCOMP_PASSIVE_USER != 2 )
      Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_USER to 2 for \"Burkert_AddColor\" !!\n" );

#  ifndef PARTICLE
   if ( Burkert_GasMFrac != 1.0 )
   {
      Burkert_GasMFrac = 1.0;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : \"Burkert_GasMFrac\" is reset to 1.0 since PARTICLE is disabled !!\n" );
   }
#  endif

// (2) set the problem-specific derived parameters
#  ifdef GRAVITY
   Burkert_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Burkert_Rho0) );
#  endif


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = (Burkert_Collision) ? 50.0 : 20.0*Burkert_FreeT;

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
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     Burkert_RSeed );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", Burkert_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", Burkert_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", Burkert_MaxR );
      Aux_Message( stdout, "  test mode                                 = %s\n",    (Burkert_Collision)?
                                                                                     "colliding clouds":"single cloud" );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "  central coordinate [%d]                   = %14.7e\n", d, Burkert_Center[d] );
      if ( Burkert_Collision ) {
      Aux_Message( stdout, "  initial distance between two clouds       = %13.7e\n", Burkert_Collision_D );
      Aux_Message( stdout, "  assign colors to different clouds         = %d\n",     Burkert_AddColor ); }
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "  bulk velocity [%d]                        = %14.7e\n", d, Burkert_BulkVel[d] );
#     if ( MODEL == HYDRO )
      Aux_Message( stdout, "  gas mass fraction                         = %13.7e\n", Burkert_GasMFrac );
#     endif
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     Burkert_MassProfNBin );
      Aux_Message( stdout, "  free-fall time at the scale radius        = %13.7e\n", Burkert_FreeT );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );
/***Deep Copy***/
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

  // gas share the same density profile as particles (except for different total masses)
   const double TotM    = 4.0/3.0*M_PI*CUBE(Burkert_R0)*Burkert_Rho0;
   const double GasRho0 = Burkert_Rho0*Burkert_GasMFrac;
   const double PresBg  = 0.0;   // background pressure (set to 0.0 by default)

   double r2, a2,a3, Dens;

   /***calculating pressure at r=0***/
   double pressure_0;
   
   /***gsl workspace***/
   
   gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
   
   double result, error;

        gsl_function F_test;
         F_test.function = &de_pressure;

         gsl_integration_qags (&F_test, 0,10000000, 0, 1e-7, 1000,
                        w, &result, &error);
         pressure_0=result;
    /***gsl workspace***/
  
  
   if ( Burkert_Collision )
   {
      const double Coll_Offset = 0.5*Burkert_Collision_D/sqrt(3.0);
      double Center[3];

      fluid[DENS] = 0.0;
      fluid[ENGY] = 0.0;
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;

      for (int t=-1; t<=1; t+=2)
      {
         for (int d=0; d<3; d++)    Center[d] = Burkert_Center[d] + Coll_Offset*(double)t;

         r2   = SQR(x-Center[0]) + SQR(y-Center[1]) + SQR(z-Center[2]);
         a2   = r2 / SQR(Burkert_R0);
         a3   = pow(a2,0.5);
          /***Main Change***/
         Dens = GasRho0 /((1+a3)*(1+a3*a3));
         /***Main Change***/

         fluid[DENS] += Dens;
#        ifdef GRAVITY
         /***Main Change***//***Not Yet***/
         /***gsl integration***/

         gsl_function F;
         F.function = &de_pressure;

         gsl_integration_qags (&F, a3,10000, 0, 1e-7, 1000,
                        w, &result, &error);

         result=pressure_0-result;
         fluid[ENGY] =(  result + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
         /*fluid[ENGY] =(  NEWTON_G*TotM*GasRho0 / ( 6.0*Burkert_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];*/
         /***gsl integration***/
         /*fluid[ENGY] = (-4*M_PI*NEWTON_G*GasRho0*Burkert_R0*Burkert_R0)*(atan(a3)/(2*a3)-log(a3*a3+1)/(4*a3)-log(a3+1)/(2*a3)+atan(a3)/2+log(a3*a3+1)/4-log(a3+1)/2)
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];*/
         /*fluid[ENGY] =(  NEWTON_G*TotM*GasRho0 / ( 6.0*Burkert_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];*/
         /***Main Change***/
#        endif

         if ( Burkert_AddColor )
         fluid[ (t==-1)?Burkert_Idx_Cloud0:Burkert_Idx_Cloud1 ] = Dens;
      }

      fluid[MOMX]  = fluid[DENS]*Burkert_BulkVel[0];
      fluid[MOMY]  = fluid[DENS]*Burkert_BulkVel[1];
      fluid[MOMZ]  = fluid[DENS]*Burkert_BulkVel[2];
      fluid[ENGY] += 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   }

   else
   {
      r2   = SQR(x-Burkert_Center[0]) + SQR(y-Burkert_Center[1]) + SQR(z-Burkert_Center[2]);
      a2   = r2 / SQR(Burkert_R0);
      a3   = pow(a2,0.5);
      /***Main Change***/
      Dens = GasRho0 /((1+a3)*(1+a3*a3));
      /***Main Change***/

      fluid[DENS] = Dens;
      fluid[MOMX] = fluid[DENS]*Burkert_BulkVel[0];
      fluid[MOMY] = fluid[DENS]*Burkert_BulkVel[1];
      fluid[MOMZ] = fluid[DENS]*Burkert_BulkVel[2];
#     ifdef GRAVITY
      /***Main Change***//***Not Yet***/
      /***gsl integration***/

         gsl_function F;
         F.function = &de_pressure;

         gsl_integration_qags (&F, 0,a3,10000, 1e-7, 1000,
                        w, &result, &error);
         result=pressure_0-result;
         fluid[ENGY] =(  result + PresBg  ) / ( GAMMA - 1.0 )+ 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
         /*fluid[ENGY] =(  NEWTON_G*TotM*GasRho0 / ( 6.0*Burkert_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];*/
         /***gsl integration***/
      //fluid[ENGY] =(  NEWTON_G*TotM*GasRho0 / ( 6.0*Burkert_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 );
      /***Main Change***/
#     endif

//    just set all passive scalars as zero
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;
   } // if ( Burkert_Collision ) ... else ...
gsl_integration_workspace_free (w);/***essential! must be added!***/

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_Burkert
// Description :  Add the problem-specific fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_Burkert()
{

   if ( Burkert_AddColor )
   {
      Burkert_Idx_Cloud0 = AddField( "Cloud0", NORMALIZE_YES );
      Burkert_Idx_Cloud1 = AddField( "Cloud1", NORMALIZE_YES );
   }
} // FUNCTION : AddNewField_Burkert
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Burkert
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
/***Deep Copy***/
void Init_TestProb_Hydro_Burkert()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Init_Field_User_Ptr      = AddNewField_Burkert;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr  = Par_Init_ByFunction_Burkert;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Burkert
/***Deep Copy***/