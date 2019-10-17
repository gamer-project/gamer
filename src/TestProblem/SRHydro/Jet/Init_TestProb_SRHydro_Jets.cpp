#include "GAMER.h"
#include "TestProb.h"

#if ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================

static int      Jet_NJet;                    // number of jets (1/2)
static double   Jet_BgDens;                  // ambient density
static double   Jet_BgVel[3];                // ambient 4-velocity
static double  *Jet_Radius        = NULL;    // radius of the cylinder-shape jet source
static double  *Jet_HalfHeight    = NULL;    // half height of the cylinder-shape jet source
static double  *Jet_SrcVel        = NULL;    // jet 4-velocity
static double  *Jet_SrcDens       = NULL;    // jet density
static double (*Jet_Vec)[3]       = NULL;    // jet orientation vector (x,y,z) (NOT necessary to be a unit vector)
static double (*Jet_CenOffset)[3] = NULL;    // jet central coordinates offset

static double   Jet_BgTemp;                  // ambient temperature
static double  *Jet_SrcTemp       = NULL;    // jet temperature
static double (*Jet_Cen)[3]       = NULL;    // jet central coordinates
static double  *Jet_WaveK         = NULL;    // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double  *Jet_MaxDis        = NULL;    // maximum distance between the cylinder-shape jet source and the jet center
                                             // --> ( Jet_Radius^2 + Jet_HalfHeight^2 )^0.5

// background variables
static bool     Jet_HSE           = false;         // hydrostatic equilibrium background
static double   Jet_HSE_Dx        = NULL_REAL;     // for Jet_HSE: distance between box and cluster centre (along x)
static double   Jet_HSE_Dy        = NULL_REAL;     // for Jet_HSE: distance between box and cluster centre (along y)
static double   Jet_HSE_Dz        = NULL_REAL;     // for Jet_HSE: distance between box and cluster centre (along z)
static char     Jet_HSE_BgTable_File[MAX_STRING];  // for Jet_HSE: filename of the background gas table
static double *Jet_HSE_BgTable_Data = NULL;        // for Jet_HSE: background gas table [radius/density/temperature]
static int     Jet_HSE_BgTable_NBin;               // for Jet_HSE: number of bins in Jet_HSE_BgTable_Data[]

// =======================================================================================

void SRHydro_Pri2Con_Double (const double In[], double Out[], const double Gamma);


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
#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
#  endif


#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
// load the number of jets first
   ReadPara->Add( "Jet_NJet",          &Jet_NJet,             -1,             1,                2                 );

// allocate memory for the variables whose size depends on the number of jets
   Jet_Radius     = new double [Jet_NJet];
   Jet_HalfHeight = new double [Jet_NJet];
   Jet_SrcVel     = new double [Jet_NJet];
   Jet_SrcDens    = new double [Jet_NJet];
   Jet_SrcTemp    = new double [Jet_NJet];
   Jet_Vec        = new double [Jet_NJet][3];
   Jet_CenOffset  = new double [Jet_NJet][3];
   Jet_Cen        = new double [Jet_NJet][3];
   Jet_WaveK      = new double [Jet_NJet];
   Jet_MaxDis     = new double [Jet_NJet];

// load the input parameters of the hydrostatic equilibrium setup
   ReadPara->Add( "Jet_HSE",             &Jet_HSE,               false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Jet_HSE_Dx",          &Jet_HSE_Dx,            0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_HSE_Dy",          &Jet_HSE_Dy,            0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_HSE_Dz",          &Jet_HSE_Dz,            0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_HSE_BgTable_File", Jet_HSE_BgTable_File,  Useless_str,   Useless_str,      Useless_str       );

// load the input parameters independent of the number of jets
   ReadPara->Add( "Jet_BgDens",        &Jet_BgDens,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_BgVel_x",       &Jet_BgVel[0],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_y",       &Jet_BgVel[1],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_z",       &Jet_BgVel[2],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgTemp",        &Jet_BgTemp,           -1.0,           Eps_double,       NoMax_double      );

// load the input parameters of each jet
   ReadPara->Add( "Jet0_Radius",       &Jet_Radius    [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_HalfHeight",   &Jet_HalfHeight[0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcVel",       &Jet_SrcVel    [0],    -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_SrcDens",      &Jet_SrcDens   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcTemp",      &Jet_SrcTemp   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_Vec_x",        &Jet_Vec       [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_y",        &Jet_Vec       [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_z",        &Jet_Vec       [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_x",  &Jet_CenOffset [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_y",  &Jet_CenOffset [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_z",  &Jet_CenOffset [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );

   if ( Jet_NJet == 2 ) {
   ReadPara->Add( "Jet1_Radius",       &Jet_Radius    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_HalfHeight",   &Jet_HalfHeight[1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcVel",       &Jet_SrcVel    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcDens",      &Jet_SrcDens   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcTemp",      &Jet_SrcTemp   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_Vec_x",        &Jet_Vec       [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_y",        &Jet_Vec       [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_z",        &Jet_Vec       [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_x",  &Jet_CenOffset [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_y",  &Jet_CenOffset [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_z",  &Jet_CenOffset [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   }

   ReadPara->Read( FileName );

   delete ReadPara;


// (1-2) check runtime parameters
   if ( Jet_HSE )
   {
   } // if ( Jet_HSE )

   else
   {
#     ifdef GRAVITY
         Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#     endif
   } // if ( Jet_HSE ) ... else ...

// (1-2) convert to code unit

   Jet_BgDens           *= 1.0       / UNIT_D;
   Jet_BgTemp           *= Const_GeV / ( Const_mp * SQR(Const_c) );

   for (int d=0; d<3; d++)
   Jet_BgVel[d]         *= 1.0       / UNIT_V;

   for (int n=0; n<Jet_NJet; n++) {
   Jet_Radius    [n]    *= Const_kpc / UNIT_L;
   Jet_HalfHeight[n]    *= Const_kpc / UNIT_L;
   Jet_SrcVel    [n]    *= 1.0       / UNIT_V;
   Jet_SrcTemp   [n]    *= Const_GeV / ( Const_mp * SQR(Const_c) );
   Jet_SrcDens   [n]    *= 1.0       / UNIT_D;

   for (int d=0; d<3; d++)
   Jet_CenOffset [n][d] *= Const_kpc / UNIT_L;
   }

// (2) set the problem-specific derived parameters

   for (int n=0; n<Jet_NJet; n++)
   {
      Jet_WaveK  [n] = 0.5*M_PI/Jet_HalfHeight[n];
      Jet_MaxDis [n] = sqrt( SQR(Jet_HalfHeight[n]) + SQR(Jet_Radius[n]) );

      for (int d=0; d<3; d++)    Jet_Cen[n][d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[n][d];
   }

// (3) load the hydrostatic equilibrium table
   if ( Jet_HSE  &&  OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
      const int  NCol         = 3;           // total number of columns to load
      const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)
      double *Table_R, *Table_D, *Table_T;

      Jet_HSE_BgTable_NBin = Aux_LoadTable( Jet_HSE_BgTable_Data, Jet_HSE_BgTable_File, NCol, Col, RowMajor_No, AllocMem_Yes );

//    convert to code units
      Table_R = Jet_HSE_BgTable_Data + 0*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + 1*Jet_HSE_BgTable_NBin;
      Table_T = Jet_HSE_BgTable_Data + 2*Jet_HSE_BgTable_NBin;

      for (int b=0; b<Jet_HSE_BgTable_NBin; b++)
      {
         Table_R[b] *= Const_kpc / UNIT_L;
         Table_D[b] *= 1.0       / UNIT_D;
         Table_T[b] *= Const_GeV / ( Const_mp * SQR(Const_c) );
      }
   } // if ( Jet_HSE  &&  OPT__INIT != INIT_BY_RESTART )

// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 50.0*Const_Myr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID      = %d\n",              TESTPROB_ID                                        );
      Aux_Message( stdout, "  Jet_BgDens           = % 14.7e g/cm^3\n",  Jet_BgDens*UNIT_D                                  );
      Aux_Message( stdout, "  Jet_BgTemp           = % 14.7e GeV\n",     Jet_BgTemp*Const_mp*SQR(Const_c)/Const_GeV         );
      Aux_Message( stdout, "  Jet_BgVel[x]         = % 14.7e c\n",       Jet_BgVel[0]*UNIT_V                                );
      Aux_Message( stdout, "  Jet_BgVel[y]         = % 14.7e c\n",       Jet_BgVel[1]*UNIT_V                                );
      Aux_Message( stdout, "  Jet_BgVel[z]         = % 14.7e c\n",       Jet_BgVel[2]*UNIT_V                                );
      Aux_Message( stdout, "  Jet_NJet             = %d\n",              Jet_NJet                                           );
      if ( Jet_HSE ) {
      Aux_Message( stdout, "  Jet_HSE              = %d\n",              Jet_HSE                                            );
      Aux_Message( stdout, "  Jet_HSE_Dx           = %e kpc\n",          Jet_HSE_Dx*Const_kpc                               );
      Aux_Message( stdout, "  Jet_HSE_Dy           = %e kpc\n",          Jet_HSE_Dy*Const_kpc                               );
      Aux_Message( stdout, "  Jet_HSE_Dz           = %e kpc\n",          Jet_HSE_Dz*Const_kpc                               );
      Aux_Message( stdout, "  Jet_HSE_BgTable_File = %s\n",              Jet_HSE_BgTable_File                               ); }
      for (int n=0; n<Jet_NJet; n++) {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Jet # %d:\n", n );
      Aux_Message( stdout, "  Jet_Radius           = % 14.7e kpc\n",     Jet_Radius    [n]*UNIT_L/Const_kpc                 );
      Aux_Message( stdout, "  Jet_HalfHeight       = % 14.7e kpc\n",     Jet_HalfHeight[n]*UNIT_L/Const_kpc                 );
      Aux_Message( stdout, "  Jet_SrcVel           = % 14.7e c\n",       Jet_SrcVel    [n]*UNIT_V                           );
      Aux_Message( stdout, "  Jet_SrcDens          = % 14.7e g/cm^3\n",  Jet_SrcDens   [n]*UNIT_D                           );
      Aux_Message( stdout, "  Jet_SrcTemp          = % 14.7e GeV\n",     Jet_SrcTemp   [n]*Const_mp*SQR(Const_c)/Const_GeV  );
      Aux_Message( stdout, "  Jet_Vec[x]           = % 14.7e\n",         Jet_Vec       [n][0]                               );
      Aux_Message( stdout, "  Jet_Vec[y]           = % 14.7e\n",         Jet_Vec       [n][1]                               );
      Aux_Message( stdout, "  Jet_Vec[z]           = % 14.7e\n",         Jet_Vec       [n][2]                               );
      Aux_Message( stdout, "  Jet_CenOffset[x]     = % 14.7e kpc\n",     Jet_CenOffset [n][0]*UNIT_L/Const_kpc              );
      Aux_Message( stdout, "  Jet_CenOffset[y]     = % 14.7e kpc\n",     Jet_CenOffset [n][1]*UNIT_L/Const_kpc              );
      Aux_Message( stdout, "  Jet_CenOffset[z]     = % 14.7e kpc\n",     Jet_CenOffset [n][2]*UNIT_L/Const_kpc              ); }
      Aux_Message( stdout, "=============================================================================\n" );

	  printf("Jet_BgDens=%e, Jet_BgTemp=%e, Jet_Radius=%e, Jet_HalfHeight=%e, Jet_SrcDens=%e, Jet_SrcTemp=%e", 
					   Jet_BgDens, Jet_BgTemp, Jet_Radius[0], Jet_HalfHeight[0], Jet_SrcDens[0], Jet_SrcTemp[0]);
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

// variables for Jet_HSE
   const double *Table_R=NULL, *Table_D=NULL, *Table_T=NULL;
   double dx, dy, dz, dr;
   const int  NCol         = 3;           // total number of columns to load
   const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)


// variables for jet
   double Pri4Vel[NCOMP_FLUID];

   Pri4Vel[1] = Jet_BgVel[0];
   Pri4Vel[2] = Jet_BgVel[1];
   Pri4Vel[3] = Jet_BgVel[2];

   if ( Jet_HSE )
   {
      Table_R = Jet_HSE_BgTable_Data + (Col[0]-1)*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + (Col[1]-1)*Jet_HSE_BgTable_NBin;
      Table_T = Jet_HSE_BgTable_Data + (Col[2]-1)*Jet_HSE_BgTable_NBin;

      dx = x - amr->BoxCenter[0] + Jet_HSE_Dx;
      dy = y - amr->BoxCenter[1] + Jet_HSE_Dy;
      dz = z - amr->BoxCenter[2] + Jet_HSE_Dz;

      dr = sqrt( dx*dx + dy*dy + dz*dz );

      Pri4Vel[0] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, dr );
      Pri4Vel[4] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, dr ) * Pri4Vel[0];

//      if ( dr <= 0.0574 )
//      {
//        Pri4Vel[0] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, dr );
//        Pri4Vel[4] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, dr ) * Pri4Vel[0];
//      }
//      else
//      {
//        Pri4Vel[0] = 1e-138;
//        Pri4Vel[4] = 1e-138;
//      }
   }
   else
   {
	  Pri4Vel[0] = Jet_BgDens;
	  Pri4Vel[4] = Jet_BgTemp * Jet_BgDens;
   }


// cast double to real
#  ifndef FLOAT8
   double Out[NCOMP_FLUID];

   SRHydro_Pri2Con_Double(Pri4Vel, Out, GAMMA);

   fluid [0] = (real) Out[0];
   fluid [1] = (real) Out[1];
   fluid [2] = (real) Out[2];
   fluid [3] = (real) Out[3];
   fluid [4] = (real) Out[4];
#  else
   SRHydro_Pri2Con_Double(Pri4Vel, fluid, GAMMA);
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Jets
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Jets()
{

   delete [] Jet_Radius;
   delete [] Jet_HalfHeight;
   delete [] Jet_SrcVel;
   delete [] Jet_SrcDens;
   delete [] Jet_Vec;
   delete [] Jet_CenOffset;
   delete [] Jet_SrcTemp;
   delete [] Jet_Cen;
   delete [] Jet_WaveK;
   delete [] Jet_MaxDis;
   delete [] Jet_HSE_BgTable_Data;

} // FUNCTION : End_Jets



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Jets
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_Jets( real fluid[], const double x, const double y, const double z, const double Time,
                                    const int lv, double AuxArray[] )
{
//   return false;
   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   double TempVec[3], _Jet_VecAbs, Jet_SrcVel_xyz[3];
   real   MomSin;


// loop over all cells to add the jet source
   for (int n=0; n<Jet_NJet; n++)
   {
//    distance: jet center to mesh
      for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[n][d];
      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );
/*
      printf("\nJet_Cen[0][0]=%e\n", Jet_Cen[0][0]);
      printf("\nJet_Cen[0][1]=%e\n", Jet_Cen[0][1]);
      printf("\nJet_Cen[0][2]=%e\n", Jet_Cen[0][2]);
      printf("\nJet_MaxDis[0]=%e\n", Jet_MaxDis[0]);
*/     
//    exclude cells far away from the jet source
      if ( Dis_c2m > Jet_MaxDis[n] )   continue;


//    vectors for calculating the distance between cells and the jet sources
      for (int d=0; d<3; d++)    TempVec[d] = Jet_Cen[n][d] + Jet_Vec[n][d];


//    velocity components along different directions
      _Jet_VecAbs = 1.0/sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );
      for (int d=0; d<3; d++)    Jet_SrcVel_xyz[d] = Jet_SrcVel[n] * Jet_Vec[n][d] * _Jet_VecAbs;


//    distance: temporary vector to mesh
      for (int d=0; d<3; d++)    Vec_v2m[d] = r[d] - TempVec[d];
      Dis_v2m = sqrt( SQR(Vec_v2m[0]) + SQR(Vec_v2m[1]) + SQR(Vec_v2m[2]) );


//    distance: jet center to temporary vector
      Dis_c2v = sqrt( SQR(Jet_Vec[n][0]) + SQR(Jet_Vec[n][1]) + SQR(Jet_Vec[n][2]) );


//    check whether or not the target cell is within the jet source
      S      = 0.5*( Dis_c2m + Dis_v2m + Dis_c2v );
      Area   = sqrt( S*(S-Dis_c2m)*(S-Dis_v2m)*(S-Dis_c2v) );
      Jet_dr = 2.0*Area/Dis_c2v;
      Jet_dh = sqrt( Dis_c2m*Dis_c2m - Jet_dr*Jet_dr );

//       reset the fluid variables within the jet source
      if ( Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
      {
//       use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
         MomSin      = sin( Jet_WaveK[n]*Jet_dh );
         MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[n][0] + Vec_c2m[1]*Jet_Vec[n][1] + Vec_c2m[2]*Jet_Vec[n][2] );

         double Pri4Vel[NCOMP_FLUID]
               = { Jet_SrcDens[n], Jet_SrcVel_xyz[0]*MomSin, Jet_SrcVel_xyz[1]*MomSin, Jet_SrcVel_xyz[2]*MomSin, Jet_SrcTemp[n]*Jet_SrcDens[n] };
 

//       cast double to real
#        ifndef FLOAT8
         double Out[NCOMP_FLUID];

         SRHydro_Pri2Con_Double(Pri4Vel, Out, GAMMA);

         fluid [0] = (real) Out[0];
         fluid [1] = (real) Out[1];
         fluid [2] = (real) Out[2];
         fluid [3] = (real) Out[3];
         fluid [4] = (real) Out[4];
#        else
         SRHydro_Pri2Con_Double(Pri4Vel, fluid, GAMMA);
#        endif

//       return immediately since we do NOT allow different jet source to overlap
         return true;

      } // if (  Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
   } // for (int n=0; n<Jet_NJet; n++)


   return false;

} // FUNCTION : Flu_ResetByUser_Jets

bool Flag_User( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{


   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const real (*Dens)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];  // density
   const real (*MomX)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX];  // momentum x
   const real (*MomY)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY];  // momentum y
   const real (*MomZ)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ];  // momentum z
   const real (*Engy)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY];  // total energy
#  ifdef GRAVITY
   const real (*Pot )[PS1][PS1] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;          // potential
#  endif

   bool Flag = false;


// put your favorite flag criteria here
// ##########################################################################################################

   /*
// Example 1 : flag if the velocity exceeds the given threshold
   const real Vel = sqrt(  ( MomX[k][j][i]*MomX[k][j][i] + MomY[k][j][i]*MomY[k][j][i] +
                             MomZ[k][j][i]*MomZ[k][j][i] )  ) / Dens[k][j][i];
   Flag = Vel > Threshold;
   */

// Example 2 : flag if the grid is within a sphere with the radius eqaul to the input "Threshold" and the origin
//             in the center of the simulation box
//   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
//   const double dr[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
//   const double Radius    = sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );


//    Flag |= ( Radius <= amr->dh[lv] * 0.5  * 1.733);
// ##########################################################################################################

//   real Cons[NCOMP_FLUID], Prim[NCOMP_FLUID];
//
//  
//   Cons[DENS] = Dens[k][j][i];
//   Cons[MOMX] = MomX[k][j][i];
//   Cons[MOMY] = MomY[k][j][i];
//   Cons[MOMZ] = MomZ[k][j][i];
//   Cons[ENGY] = Engy[k][j][i];
//
//
//   double Temperature =
//   SRHydro_GetTemperature( Dens[k][j][i], MomX[k][j][i], MomY[k][j][i], MomZ[k][j][i], Engy[k][j][i], GAMMA, MIN_TEMP );
//   Flag |= Temperature >= Threshold;
//   double LorentzFactor;
//   SRHydro_Con2Pri(Cons, Prim);

// ##########################################################################################################

  double StartTime = 600;
  double EndPtVel  = 0.196;  // the velocity of jet head
  double IniX      = 375.0-202.0; // the distance between source and jet head

  double Shift = +EndPtVel*(Time[lv]-StartTime)+IniX;

  const double Center[3]    = { 0.5*amr->BoxSize[0]      , 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
  const double EndPointR[3] = { 0.5*amr->BoxSize[0]+Shift, 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
  const double EndPointL[3] = { 0.5*amr->BoxSize[0]-Shift, 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

  const double dr0[3]     = { Pos[0]-   Center[0], Pos[1]-   Center[1], Pos[2]-   Center[2] };
  const double dr1[3]     = { Pos[0]-EndPointR[0], Pos[1]-EndPointR[1], Pos[2]-EndPointR[2] };
  const double dr2[3]     = { Pos[0]-EndPointL[0], Pos[1]-EndPointL[1], Pos[2]-EndPointL[2] };

  Flag |= sqrt( dr0[1]*dr0[1] + dr0[2]*dr0[2] ) < 1.5;
  Flag &=  abs( dr0[0] ) < Shift;

  double Radius1 = sqrt( dr1[0]*dr1[0] + dr1[1]*dr1[1] + dr1[2]*dr1[2] );
  double Radius2 = sqrt( dr2[0]*dr2[0] + dr2[1]*dr2[1] + dr2[2]*dr2[2] );


  Flag |= ( Radius1 < amr->dh[lv] * 0.5  * 1.733 );
  Flag |= ( Radius2 < amr->dh[lv] * 0.5  * 1.733 );
  Flag |= ( Radius1 < amr->dh[lv] * 0.5  *   1.1 );
  Flag |= ( Radius2 < amr->dh[lv] * 0.5  *   1.1 );

// ##########################################################################################################
//   const double Center[3]    = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
//   const double dr0[3]     = { Pos[0]-   Center[0], Pos[1]-   Center[1], Pos[2]-   Center[2] };
//   double BOX_SIZE_X = 750.0;
//   double BOX_SIZE_Y = 187.50;
//   double BOX_SIZE_Z = 187.50;
//
//   if ( lv==0 )
//   {
//     Flag |=  ( dr0[0] < 100.0 );
//     Flag |=  ( dr0[0] > 700.0 );
//   }
//   else if ( lv == 2 )
//   {
//     Flag |=  ( dr0[0] > BOX_SIZE_X *3.0 / 8.0 + amr->dh[lv] ) & ( dr0[0] < BOX_SIZE_X *3.0 / 8.0               );
//     Flag &=  ( dr0[1] > BOX_SIZE_Y      / 2.0 - amr->dh[lv] ) & ( dr0[1] < BOX_SIZE_Y      / 2.0 + amr->dh[lv] );
//     Flag &=  ( dr0[2] > BOX_SIZE_Z      / 2.0 - amr->dh[lv] ) & ( dr0[2] < BOX_SIZE_Z      / 2.0 + amr->dh[lv] );
//
//     Flag |=  ( dr0[0] > BOX_SIZE_X *5.0 / 8.0               ) & ( dr0[0] < BOX_SIZE_X *5.0 / 8.0 - amr->dh[lv] );
//     Flag &=  ( dr0[1] > BOX_SIZE_Y      / 2.0 - amr->dh[lv] ) & ( dr0[1] < BOX_SIZE_Y      / 2.0 + amr->dh[lv] );
//     Flag &=  ( dr0[2] > BOX_SIZE_Z      / 2.0 - amr->dh[lv] ) & ( dr0[2] < BOX_SIZE_Z      / 2.0 + amr->dh[lv] );
//   }


   return Flag;

} // FUNCTION : Flag_User

#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAcc
// Description :  Set ExtAcc_AuxArray[] used by the external acceration routine ExternalAcc()
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_ExternalAcc_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__GRAVITY_TYPE == 2/3"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAcc()
{

// ExtAcc_AuxArray has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)
// --> by default we set
//     ExtAcc_AuxArray[0] = x coordinate of the external acceleration center
//     ExtAcc_AuxArray[1] = y ...
//     ExtAcc_AuxArray[2] = z ..
//     ExtAcc_AuxArray[3] = gravitational_constant*point_source_mass
//     ExtAcc_AuxArray[4] = soften_length (<=0.0 --> disable)
// --> to change the this default behavior, please edit "GPU_Gravity/CUPOT_ExternalAcc.cu"

   const double M   = 1e-3;
   const double GM  = NEWTON_G*M;
   const double Eps = 0.0;

   ExtAcc_AuxArray[0] = 0.5*amr->BoxSize[0] - Jet_HSE_Dx;
   ExtAcc_AuxArray[1] = 0.5*amr->BoxSize[1] - Jet_HSE_Dy;
   ExtAcc_AuxArray[2] = 0.5*amr->BoxSize[2] - Jet_HSE_Dz;
   ExtAcc_AuxArray[3] = GM;
   ExtAcc_AuxArray[4] = Eps;

} // FUNCTION : Init_ExternalAcc
#endif


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_Jets
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_Jets()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_User;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Jets;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = End_Jets;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_Jets

#endif
