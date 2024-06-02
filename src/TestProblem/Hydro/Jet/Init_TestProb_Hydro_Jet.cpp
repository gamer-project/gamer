#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int      Jet_NJet;                          // number of jets (1/2)
static double   Jet_BgDens;                        // ambient density
static double   Jet_BgTemp;                        // ambient temperature
       double   Jet_BgVel[3];                      // ambient velocity
static double  *Jet_Radius        = NULL;          // radius of the cylinder-shape jet source
static double  *Jet_HalfHeight    = NULL;          // half height of the cylinder-shape jet source
static double  *Jet_SrcVel        = NULL;          // jet velocity
static double  *Jet_SrcDens       = NULL;          // jet density
static double  *Jet_SrcTemp       = NULL;          // jet temperature
static double (*Jet_Vec)[3]       = NULL;          // jet orientation vector (x,y,z) (NOT necessary to be a unit vector)
static double (*Jet_CenOffset)[3] = NULL;          // jet central coordinates offset
static bool     Jet_HSE           = false;         // hydrostatic equilibrium background
       double   Jet_HSE_D         = NULL_REAL;     // for Jet_HSE: distance between box and cluster centre (assuming along y)
       double   Jet_HSE_M200      = NULL_REAL;     // for Jet_HSE: cluster virial mass
       double   Jet_HSE_R200      = NULL_REAL;     // for Jet_HSE: cluster virial radius
       double   Jet_HSE_C200      = NULL_REAL;     // for Jet_HSE: cluster concentration parameter assuming NFW
static char     Jet_HSE_BgTable_File[MAX_STRING];  // for Jet_HSE: filename of the background gas table

static double   Jet_BgEint;                        // ambient internal energy
static double  *Jet_SrcEint       = NULL;          // jet internal energy
static double (*Jet_Cen)[3]       = NULL;          // jet central coordinates
static double  *Jet_WaveK         = NULL;          // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double  *Jet_MaxDis        = NULL;          // maximum distance between the cylinder-shape jet source and the jet center
                                                   // --> ( Jet_Radius^2 + Jet_HalfHeight^2 )^0.5
static double *Jet_HSE_BgTable_Data = NULL;        // for Jet_HSE: background gas table [radius/density/temperature]
static int     Jet_HSE_BgTable_NBin;               // for Jet_HSE: number of bins in Jet_HSE_BgTable_Data[]
// =======================================================================================

#ifdef GRAVITY
void Init_ExtAcc_Jet();
#endif
static void IsolatedBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                        const int GhostSize, const int idx[], const double pos[], const double Time,
                        const int lv, const int TFluVarIdxList[], double AuxArray[] );




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
   ReadPara->Read( FileName );
   delete ReadPara;
   ReadPara = new ReadPara_t; // for loading the remaining parameters

// allocate memory for the variables whose size depends on the number of jets
   Jet_Radius     = new double [Jet_NJet];
   Jet_HalfHeight = new double [Jet_NJet];
   Jet_SrcVel     = new double [Jet_NJet];
   Jet_SrcDens    = new double [Jet_NJet];
   Jet_SrcTemp    = new double [Jet_NJet];
   Jet_Vec        = new double [Jet_NJet][3];
   Jet_CenOffset  = new double [Jet_NJet][3];
   Jet_SrcEint    = new double [Jet_NJet];
   Jet_Cen        = new double [Jet_NJet][3];
   Jet_WaveK      = new double [Jet_NJet];
   Jet_MaxDis     = new double [Jet_NJet];

// load the input parameters independent of the number of jets
   ReadPara->Add( "Jet_BgDens",          &Jet_BgDens,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_BgTemp",          &Jet_BgTemp,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_BgVel_x",         &Jet_BgVel[0],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_y",         &Jet_BgVel[1],          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet_BgVel_z",         &Jet_BgVel[2],          0.0,           NoMin_double,     NoMax_double      );

// load the input parameters of the hydrostatic equilibrium setup
   ReadPara->Add( "Jet_HSE",             &Jet_HSE,               false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Jet_HSE_D",           &Jet_HSE_D,             NoDef_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_HSE_M200",        &Jet_HSE_M200,          NoDef_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_HSE_R200",        &Jet_HSE_R200,          NoDef_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_HSE_C200",        &Jet_HSE_C200,          NoDef_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet_HSE_BgTable_File", Jet_HSE_BgTable_File,  NoDef_str,     Useless_str,      Useless_str       );

// load the input parameters of each jet
   ReadPara->Add( "Jet0_Radius",         &Jet_Radius    [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_HalfHeight",     &Jet_HalfHeight[0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcVel",         &Jet_SrcVel    [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcDens",        &Jet_SrcDens   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_SrcTemp",        &Jet_SrcTemp   [0],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet0_Vec_x",          &Jet_Vec       [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_y",          &Jet_Vec       [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_Vec_z",          &Jet_Vec       [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_x",    &Jet_CenOffset [0][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_y",    &Jet_CenOffset [0][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet0_CenOffset_z",    &Jet_CenOffset [0][2],  NoDef_double,  NoMin_double,     NoMax_double      );

   if ( Jet_NJet == 2 ) {
   ReadPara->Add( "Jet1_Radius",         &Jet_Radius    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_HalfHeight",     &Jet_HalfHeight[1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcVel",         &Jet_SrcVel    [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcDens",        &Jet_SrcDens   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_SrcTemp",        &Jet_SrcTemp   [1],    -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jet1_Vec_x",          &Jet_Vec       [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_y",          &Jet_Vec       [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_Vec_z",          &Jet_Vec       [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_x",    &Jet_CenOffset [1][0],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_y",    &Jet_CenOffset [1][1],  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jet1_CenOffset_z",    &Jet_CenOffset [1][2],  NoDef_double,  NoMin_double,     NoMax_double      );
   }

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) check runtime parameters
   if ( Jet_HSE )
   {
#     ifndef GRAVITY
         Aux_Error( ERROR_INFO, "GRAVITY must be enabled for Jet_HSE !!\n" );
#     endif

#     ifdef GRAVITY
      if ( OPT__SELF_GRAVITY  &&  MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : are you sure about enabling OPT__SELF_GRAVITY for Jet_HSE !?\n" );

      if ( ! OPT__EXT_ACC )
         Aux_Error( ERROR_INFO, "must enable OPT__EXT_ACC for Jet_HSE !!\n" );
#     endif

      if ( Jet_HSE_D == NoDef_double )
         Aux_Error( ERROR_INFO, "\"Jet_HSE_D\" is not set for Jet_HSE !!\n" );

      if ( Jet_HSE_M200 == NoDef_double )
         Aux_Error( ERROR_INFO, "\"Jet_HSE_M200\" is not set for Jet_HSE !!\n" );

      if ( Jet_HSE_R200 == NoDef_double )
         Aux_Error( ERROR_INFO, "\"Jet_HSE_R200\" is not set for Jet_HSE !!\n" );

      if ( Jet_HSE_C200 == NoDef_double )
         Aux_Error( ERROR_INFO, "\"Jet_HSE_C200\" is not set for Jet_HSE !!\n" );
   } // if ( Jet_HSE )

   else
   {
#     ifdef GRAVITY
         Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#     endif
   } // if ( Jet_HSE ) ... else ...

// (1-3) convert to code units
   Jet_BgDens           *= 1.0        / UNIT_D;
   Jet_BgTemp           *= Const_keV  / UNIT_E;
   for (int d=0; d<3; d++)
   Jet_BgVel[d]         *= 1.0        / UNIT_V;

   for (int n=0; n<Jet_NJet; n++) {
   Jet_Radius    [n]    *= Const_kpc  / UNIT_L;
   Jet_HalfHeight[n]    *= Const_kpc  / UNIT_L;
   Jet_SrcVel    [n]    *= 1.0        / UNIT_V;
   Jet_SrcDens   [n]    *= 1.0        / UNIT_D;
   Jet_SrcTemp   [n]    *= Const_keV  / UNIT_E;
   for (int d=0; d<3; d++)
   Jet_CenOffset [n][d] *= Const_kpc  / UNIT_L;
   }

   if ( Jet_HSE ) {
   Jet_HSE_D            *= Const_kpc  / UNIT_L;
   Jet_HSE_M200         *= Const_Msun / UNIT_M;
   Jet_HSE_R200         *= Const_kpc  / UNIT_L;
   }


// (2) set the problem-specific derived parameters
// must initialize EoS first
   EoS_Init();

// assuming EoS requires no passive scalars
   Jet_BgEint = EoS_DensPres2Eint_CPUPtr(  Jet_BgDens,
                                           EoS_DensTemp2Pres_CPUPtr( Jet_BgDens, Jet_BgTemp*UNIT_E/Const_kB, NULL,
                                                                     EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ),
                                           NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table  );

   for (int n=0; n<Jet_NJet; n++)
   {
      Jet_WaveK  [n] = 0.5*M_PI/Jet_HalfHeight[n];
      Jet_MaxDis [n] = sqrt( SQR(Jet_HalfHeight[n]) + SQR(Jet_Radius[n]) );
//    assuming EoS requires no passive scalars
      Jet_SrcEint[n] = EoS_DensPres2Eint_CPUPtr(  Jet_SrcDens[n],
                                                  EoS_DensTemp2Pres_CPUPtr( Jet_SrcDens[n], Jet_SrcTemp[n]*UNIT_E/Const_kB, NULL,
                                                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ),
                                                  NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table  );

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
         Table_T[b] *= Const_keV / UNIT_E;

//       convert temperature to internal energy (assuming EoS requires no passive scalars)
         Table_T[b]  = EoS_DensPres2Eint_CPUPtr(  Table_D[b],
                                                  EoS_DensTemp2Pres_CPUPtr( Table_D[b], Table_T[b]*UNIT_E/Const_kB, NULL,
                                                                            EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ),
                                                  NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table  );
      }
   } // if ( Jet_HSE  &&  OPT__INIT != INIT_BY_RESTART )


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 50.0*Const_Myr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID      = %d\n",             TESTPROB_ID                           );
      Aux_Message( stdout, "  Jet_BgDens           = % 14.7e g/cm^3\n", Jet_BgDens*UNIT_D                     );
      Aux_Message( stdout, "  Jet_BgTemp           = % 14.7e keV\n",    Jet_BgTemp*UNIT_E/Const_keV           );
      Aux_Message( stdout, "  Jet_BgVel[x]         = % 14.7e cm/s\n",   Jet_BgVel[0]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_BgVel[y]         = % 14.7e cm/s\n",   Jet_BgVel[1]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_BgVel[z]         = % 14.7e cm/s\n",   Jet_BgVel[2]*UNIT_V                   );
      Aux_Message( stdout, "  Jet_HSE              = %d\n",             Jet_HSE                               );
      if ( Jet_HSE ) {
      Aux_Message( stdout, "  Jet_HSE_D            = % 14.7e kpc\n",    Jet_HSE_D*UNIT_L/Const_kpc            );
      Aux_Message( stdout, "  Jet_HSE_M200         = % 14.7e Msun\n",   Jet_HSE_M200*UNIT_M/Const_Msun        );
      Aux_Message( stdout, "  Jet_HSE_R200         = % 14.7e kpc\n",    Jet_HSE_R200*UNIT_L/Const_kpc         );
      Aux_Message( stdout, "  Jet_HSE_C200         = % 14.7e\n",        Jet_HSE_C200                          );
      Aux_Message( stdout, "  Jet_HSE_BgTable_File = %s\n",             Jet_HSE_BgTable_File                  ); }
      Aux_Message( stdout, "  Jet_NJet             = %d\n",             Jet_NJet                              );
      for (int n=0; n<Jet_NJet; n++) {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Jet # %d:\n", n );
      Aux_Message( stdout, "     Jet_Radius        = % 14.7e kpc\n",    Jet_Radius    [n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "     Jet_HalfHeight    = % 14.7e kpc\n",    Jet_HalfHeight[n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "     Jet_SrcVel        = % 14.7e cm/s\n",   Jet_SrcVel    [n]*UNIT_V              );
      Aux_Message( stdout, "     Jet_SrcDens       = % 14.7e g/cm^3\n", Jet_SrcDens   [n]*UNIT_D              );
      Aux_Message( stdout, "     Jet_SrcTemp       = % 14.7e keV\n",    Jet_SrcTemp   [n]*UNIT_E/Const_keV    );
      Aux_Message( stdout, "     Jet_Vec[x]        = % 14.7e\n",        Jet_Vec       [n][0]                  );
      Aux_Message( stdout, "     Jet_Vec[y]        = % 14.7e\n",        Jet_Vec       [n][1]                  );
      Aux_Message( stdout, "     Jet_Vec[z]        = % 14.7e\n",        Jet_Vec       [n][2]                  );
      Aux_Message( stdout, "     Jet_CenOffset[x]  = % 14.7e kpc\n",    Jet_CenOffset [n][0]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "     Jet_CenOffset[y]  = % 14.7e kpc\n",    Jet_CenOffset [n][1]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "     Jet_CenOffset[z]  = % 14.7e kpc\n",    Jet_CenOffset [n][2]*UNIT_L/Const_kpc ); }
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

// variables for Jet_HSE
   const double *Table_R=NULL, *Table_D=NULL, *Table_E=NULL;
   double dy;

   if ( Jet_HSE )
   {
      Table_R = Jet_HSE_BgTable_Data + 0*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + 1*Jet_HSE_BgTable_NBin;
      Table_E = Jet_HSE_BgTable_Data + 2*Jet_HSE_BgTable_NBin;
      dy      = y - amr->BoxCenter[1] + Jet_HSE_D - Jet_BgVel[1]*Time;

      fluid[DENS] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, dy );
   }

   else
      fluid[DENS] = Jet_BgDens;

   fluid[MOMX] = fluid[DENS]*Jet_BgVel[0];
   fluid[MOMY] = fluid[DENS]*Jet_BgVel[1];
   fluid[MOMZ] = fluid[DENS]*Jet_BgVel[2];
   fluid[ENGY] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

   if ( Jet_HSE )
      fluid[ENGY] += Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_E, dy );
   else
      fluid[ENGY] += Jet_BgEint;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_Jet
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Jet()
{

   delete [] Jet_Radius;
   delete [] Jet_HalfHeight;
   delete [] Jet_SrcVel;
   delete [] Jet_SrcDens;
   delete [] Jet_SrcTemp;
   delete [] Jet_Vec;
   delete [] Jet_CenOffset;
   delete [] Jet_SrcEint;
   delete [] Jet_Cen;
   delete [] Jet_WaveK;
   delete [] Jet_MaxDis;
   delete [] Jet_HSE_BgTable_Data;

} // FUNCTION : End_Jet



//-------------------------------------------------------------------------------------------------------
// Function    :  IsolatedBC
// Description :  Isolated boundary condition for galaxies, only allow outflow velocities
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void IsolatedBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                 const int GhostSize, const int idx[], const double pos[], const double Time,
                 const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// simply call the IC function
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : IsolatedBC



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Jet
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
//                Emag     : Magnetic energy (MHD only)
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                dt       : Time interval to advance solution
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
int Flu_ResetByUser_Jet( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                         const double dt, const int lv, double AuxArray[] )
{

   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   double TempVec[3], _Jet_VecAbs, Jet_SrcVel_xyz[3];
   real   MomSin;


// loop over all jet sources
   for (int n=0; n<Jet_NJet; n++)
   {
//    distance: jet center to mesh
      for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[n][d];
      Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );


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

      if ( Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
      {
//       reset the fluid variables within the jet source
         fluid[DENS] = Jet_SrcDens[n];

//       use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
         MomSin      = fluid[DENS]*sin( Jet_WaveK[n]*Jet_dh );
         MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[n][0] + Vec_c2m[1]*Jet_Vec[n][1] + Vec_c2m[2]*Jet_Vec[n][2] );
         fluid[MOMX] = MomSin*Jet_SrcVel_xyz[0];
         fluid[MOMY] = MomSin*Jet_SrcVel_xyz[1];
         fluid[MOMZ] = MomSin*Jet_SrcVel_xyz[2];

//       assuming no magnetic field
         fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Jet_SrcEint[n], 0.0 );

//       return immediately since we do NOT allow different jet source to overlap
         return true;

      } // if (  Jet_dh <= Jet_HalfHeight[n]  &&  Jet_dr <= Jet_Radius[n] )
   } // for (int n=0; n<Jet_NJet; n++)


   return false;

} // FUNCTION : Flu_ResetByUser_Jet
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Jet
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Jet()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = IsolatedBC;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Jet;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = End_Jet;
#  ifdef GRAVITY
   if ( OPT__EXT_ACC == EXT_ACC_FUNC )
   Init_ExtAcc_Ptr         = Init_ExtAcc_Jet;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Jet
