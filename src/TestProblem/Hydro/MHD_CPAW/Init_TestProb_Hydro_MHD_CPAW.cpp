#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double MHD_CPAW_Rho0;          // background density
static double MHD_CPAW_A;             // perturbation amplitude
static double MHD_CPAW_P0;            // background pressure
static double MHD_CPAW_B0;            // background magnetic field
static int    MHD_CPAW_Dir;           // wave direction: (0/1/2/3) --> (x/y/z/diagonal)

static double MHD_CPAW_WaveSpeed;     // propagation speed
static double MHD_CPAW_WaveLength;    // wavelength
static double MHD_CPAW_WaveNumber;    // wavenumber
static double MHD_CPAW_WaveFrequency; // wave frequency
// =======================================================================================




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

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  if ( EOS != EOS_GAMMA )
   Aux_Error( ERROR_INFO, "EOS != EOS_GAMMA !!\n" );
#  endif

   if ( OPT__INIT_BFIELD_BYVECPOT != INIT_MAG_BYVECPOT_FUNC )
      Aux_Error( ERROR_INFO, "Must set OPT__INIT_BFIELD_BYVECPOT != %d !!\n", INIT_MAG_BYVECPOT_FUNC );

   if ( MHD_CPAW_Dir == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for MHD_CPAW_Dir = %d !!\n", MHD_CPAW_Dir );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC (i.e., \"OPT__BC_FLU_* = 1\") !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for this test !!\n" );
#     endif

      if ( !OPT__OUTPUT_USER )   Aux_Message( stdout, "WARNING : OPT__OUTPUT_USER is off !!\n" );
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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "MHD_CPAW_Rho0",    &MHD_CPAW_Rho0,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "MHD_CPAW_A",       &MHD_CPAW_A,           -1.0,          0.0,              NoMax_double      );
   ReadPara->Add( "MHD_CPAW_P0",      &MHD_CPAW_P0,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "MHD_CPAW_Dir",     &MHD_CPAW_Dir,          3,            0,                3                 );
   ReadPara->Add( "MHD_CPAW_B0",      &MHD_CPAW_B0,          -1.0,          0.0,              NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters
   MHD_CPAW_WaveLength = ( MHD_CPAW_Dir == 3 ) ? amr->BoxSize[0]/sqrt(3.0) : amr->BoxSize[MHD_CPAW_Dir];
   MHD_CPAW_WaveNumber    = 2.0*M_PI/MHD_CPAW_WaveLength;
   MHD_CPAW_WaveSpeed     = MHD_CPAW_B0/SQRT(MHD_CPAW_Rho0);
   MHD_CPAW_WaveFrequency = MHD_CPAW_WaveNumber*MHD_CPAW_WaveSpeed;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 2.0*M_PI/MHD_CPAW_WaveFrequency;
   const long   End_Step_Default = __INT_MAX__;

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
      Aux_Message( stdout, "  test problem ID        = %d\n",      TESTPROB_ID            );
      Aux_Message( stdout, "  background density     = % 14.7e\n", MHD_CPAW_Rho0          );
      Aux_Message( stdout, "  perturbation amplitude = % 14.7e\n", MHD_CPAW_A             );
      Aux_Message( stdout, "  background pressure    = % 14.7e\n", MHD_CPAW_P0            );
      Aux_Message( stdout, "  background B field     = % 14.7e\n", MHD_CPAW_B0            );
      Aux_Message( stdout, "  direction              = %d\n",      MHD_CPAW_Dir           );
      Aux_Message( stdout, "  wave speed             = % 14.7e\n", MHD_CPAW_WaveSpeed     );
      Aux_Message( stdout, "  wavelength             = % 14.7e\n", MHD_CPAW_WaveLength    );
      Aux_Message( stdout, "  wavenumber             = % 14.7e\n", MHD_CPAW_WaveNumber    );
      Aux_Message( stdout, "  wave frequency         = % 14.7e\n", MHD_CPAW_WaveFrequency );
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

// short names
   const double Rho0       = MHD_CPAW_Rho0;
   const double A          = MHD_CPAW_A;
   const double P0         = MHD_CPAW_P0;
   const double WaveSpeed  = MHD_CPAW_WaveSpeed;
   const double WaveLength = MHD_CPAW_WaveLength;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;

   double kr, WaveForm1, WaveForm2, Gamma, DampForm;
   double Dens, Mom0, Mom1, Mom2, MomX, MomY, MomZ, Pres, Eint, Etot;

   switch ( MHD_CPAW_Dir ) {
      case 0:  kr = WaveK*x;                        break;
      case 1:  kr = WaveK*y;                        break;
      case 2:  kr = WaveK*z;                        break;
      case 3:  kr = WaveK*( x + y + z )/sqrt(3.0);  break;
   }

   Dens = Rho0;
   Pres = P0;

   WaveForm1 =  sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   Mom0      =  Dens * WaveSpeed;
   Mom1      = -Mom0 * A * WaveForm1;
   Mom2      = -Mom0 * A * WaveForm2;

   switch ( MHD_CPAW_Dir ) {
      case 0:  
         MomX = Mom0;                
         MomY = Mom1; 
         MomZ = Mom2; 
         break;
      case 1:  
         MomX = Mom2;     
         MomY = Mom0;            
         MomZ = Mom1; 
         break;
      case 2:  
         MomX = Mom1;    
         MomY = Mom2; 
         MomZ = Mom0;           
         break;
      case 3:  
         MomX = Mom0/sqrt(3.0) + Mom1/sqrt(2.0) + Mom2/sqrt(6.0);  
         MomY = Mom0/sqrt(3.0) - Mom1/sqrt(2.0) + Mom2/sqrt(6.0);           
         MomZ = Mom0/sqrt(3.0) - 2.0*Mom2/sqrt(6.0);       
         break;
   }
   
// compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetAFieldIC
// Description :  Function template to initialize the magnetic vector potential
//
// Note        :  1. Invoked by MHD_Init_BField_ByVecPot_Function() using the function pointer
//                   "Init_BField_ByVecPot_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target refinement level
//                Component : Component of the output magnetic vector potential
//                            --> Supported components: 'x', 'y', 'z'
//                AuxArray  : Auxiliary array
//                            --> Useless since it is currently fixed to NULL
//
// Return      :  "Component"-component of the magnetic vector potential at (x, y, z, Time)
//-------------------------------------------------------------------------------------------------------
double SetAFieldIC( const double x, const double y, const double z, const double Time,
                    const int lv, const char Component, double AuxArray[] )
{

// short names
   const double A          = MHD_CPAW_A;
   const double B0         = MHD_CPAW_B0;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;

   double kr, WaveForm1, WaveForm2;
   double e1x, e1y, e1z, e2x, e2y, e2z;
   double A1, A2, mag_vecpot;
   double xcomp, ycomp, zcomp;

   switch ( MHD_CPAW_Dir ) {
      case 0:  
         kr    =  WaveK*x;
         e1x   =  0.0;           e1y   =  1.0;           e1z   =  0.0;
         e2x   =  0.0;           e2y   =  0.0;           e2z   =  1.0;   
         xcomp =  0.0;           ycomp =  0.0;           zcomp =  y;                     
         break;
      case 1:  
         kr    =  WaveK*y;  
         e1x   =  0.0;           e1y   =  0.0;           e1z   =  1.0;
         e2x   =  1.0;           e1y   =  0.0;           e2z   =  0.0;       
         xcomp =  z;             ycomp =  0.0;           zcomp =  0.0;                                          
         break;
      case 2:  
         kr    =  WaveK*z;       
         e1x   =  1.0;           e1y   =  0.0;           e1z   =  0.0;
         e2x   =  0.0;           e1y   =  1.0;           e2z   =  0.0;      
         xcomp =  0.0;           ycomp = x;              zcomp =  0.0;                                          
         break;
      case 3:  
         kr    =  WaveK*( x + y + z )/sqrt(3.0);  
         e1x   =  1.0/sqrt(2.0); e1y   = -1.0/sqrt(2.0);  e1z   =  0.0;
         e2x   =  1.0/sqrt(6.0); e2y   =  1.0/sqrt(6.0);  e2z   = -2.0/sqrt(6.0);     
         xcomp =  z/sqrt(3.0);    
         ycomp =  x/sqrt(3.0);    
         zcomp =  y/sqrt(3.0);                                          
         break;
   }
   
   WaveForm1 = -sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   A1        =  B0 * A * WaveForm1 / WaveK;
   A2        =  B0 * A * WaveForm2 / WaveK;

   switch ( Component )
   {
      case 'x' :  mag_vecpot = B0*xcomp+A1*e1x+A2*e2x;  break;
      case 'y' :  mag_vecpot = B0*ycomp+A1*e1y+A2*e2y;  break;
      case 'z' :  mag_vecpot = B0*zcomp+A1*e1z+A2*e2z;  break;
      default  :  Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
   }


   return mag_vecpot;

} // FUNCTION : SetAFieldIC

//-------------------------------------------------------------------------------------------------------
// Function    :  CheckBField
// Description :  Check the problem-specific magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void CheckBField( real magnetic[], const double x, const double y, const double z, const double Time,
                    const int lv, double AuxArray[] )
{


// short names
   const double A          = MHD_CPAW_A;
   const double B0         = MHD_CPAW_B0;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;

   double kr, B1, B2, WaveForm1, WaveForm2;

   switch ( MHD_CPAW_Dir ) {
      case 0:  kr = WaveK*x;                        break;
      case 1:  kr = WaveK*y;                        break;
      case 2:  kr = WaveK*z;                        break;
      case 3:  kr = WaveK*( x + y + z )/sqrt(3.0);  break;
   }

   WaveForm1 =  sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   B1        =  B0 * A * WaveForm1;
   B2        =  B0 * A * WaveForm2;

   switch ( MHD_CPAW_Dir ) {
      case 0: magnetic[MAGX] = B0;
              magnetic[MAGY] = B1;
              magnetic[MAGZ] = B2;
              break;

      case 1: magnetic[MAGX] = B2;
              magnetic[MAGY] = B0;
              magnetic[MAGZ] = B1;
              break;

      case 2: magnetic[MAGX] = B1;
              magnetic[MAGY] = B2;
              magnetic[MAGZ] = B0;
              break;

      case 3: magnetic[MAGX] = B0/sqrt(3.0) + B1/sqrt(2.0) + B2/sqrt(6.0);
              magnetic[MAGY] = B0/sqrt(3.0) - B1/sqrt(2.0) + B2/sqrt(6.0); 
              magnetic[MAGZ] = B0/sqrt(3.0) - 2.0*B2/sqrt(6.0);     
              break;
   }

} // FUNCTION : CheckBField

#endif // #ifdef MHD


//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void OutputError()
{

   const char Prefix[100]     = "MHD_CPAW";
   const OptOutputPart_t Part = OUTPUT_X + MHD_CPAW_Dir;

   Output_L1Error( SetGridIC, CheckBField, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_MHD_CPAW
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_MHD_CPAW()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_BField_ByVecPot_User_Ptr = SetAFieldIC;
   Output_User_Ptr               = OutputError;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_MHD_CPAW
