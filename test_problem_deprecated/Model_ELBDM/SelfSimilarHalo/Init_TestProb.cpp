#include "GAMER.h"

#if ( MODEL == ELBDM )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void ELBDM_TestProbSol_SelSimHalo( real *fluid, const double x, const double y, const double z, const double Time );
static void LoadTestProbParameter();
static void LoadProfile( const char *FileName, double **Profile, int &NBin, const bool LoadPhase );
static bool CheckEmptyString( const char *InputString );
static int  CountRow( const char *FileName );


// global variables in the self-similar halo test
// =======================================================================================
char    SelSimHalo_FileName[1000];              // name of the profile table
bool    SelSimHalo_LoadPhase;                   // true  --> load the phase of wave function from the profile table directly
                                                // false --> load velocity and integrate to get phase

double *SelSimHalo_Prof[3] = {NULL,NULL,NULL};  // dimensionless radial profiles [eta/mass_density/phase]
int     SelSimHalo_NBin;                        // number of radial bins in the profile table
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the self-similar halo test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "self-similar halo";

// check
#  if ( MODEL != ELBDM )
#  error : ERROR : only support the ELBDM model !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the self-similar halo test !!
#  endif

#  ifndef COMOVING
#  error : ERROR : "COMOVING must be ON" in the self-similar halo test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = ELBDM_TestProbSol_SelSimHalo;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// load the radial profile
   if ( OPT__INIT != INIT_RESTART )
      LoadProfile( SelSimHalo_FileName, SelSimHalo_Prof, SelSimHalo_NBin, SelSimHalo_LoadPhase );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "   profile file = %s\n", SelSimHalo_FileName  );
      Aux_Message( stdout, "   load phase   = %d\n", SelSimHalo_LoadPhase );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 1.0/3.0;     // z = 2.0
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %ld in the %s test\n", "END_STEP", END_STEP, TestProb );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %13.7e Gyr in the %s test\n", "END_T", END_T*UNIT_T/Const_Gyr, TestProb );
   }

} // FUNCTION : Init_TestProb



//------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_TestProbSol_SelSimHalo
// Description :  Initialize the ELBDM self-similar halo test
//
// Note        :  Invoked by "ELBDM_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void ELBDM_TestProbSol_SelSimHalo( real *fluid, const double x, const double y, const double z, const double Time )
{

   const double  Cen[3]      = { 0.5*amr->BoxSize[0],
                                 0.5*amr->BoxSize[1],
                                 0.5*amr->BoxSize[2] };
   const int     MinIdx      = 0;
   const int     MaxIdx      = SelSimHalo_NBin - 1;
   const double *Prof_Radius = SelSimHalo_Prof[0];
   const double *Prof_Dens   = SelSimHalo_Prof[1];
   const double *Prof_Phase  = SelSimHalo_Prof[2];

   double r, Dens, Phase;

   r     = sqrt( SQR(x-Cen[0]) + SQR(y-Cen[1]) + SQR(z-Cen[2]) );
   Dens  = Mis_InterpolateFromTable( SelSimHalo_NBin, Prof_Radius, Prof_Dens,  r );
   Phase = Mis_InterpolateFromTable( SelSimHalo_NBin, Prof_Radius, Prof_Phase, r );

   if ( Dens == NULL_REAL  ||  Phase == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (probably outside the input table)!!\n", r );

   fluid[DENS] = Dens;
   fluid[REAL] = sqrt( Dens )*cos( Phase );
   fluid[IMAG] = sqrt( Dens )*sin( Phase );

} // FUNCTION : ELBDM_TestProbSol_SelSimHalo



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter
// Description :  Load parameters for the test problem
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter()
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;
   int    tmp_int;


// skip the header
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    SelSimHalo_FileName, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &tmp_int,             string );
   SelSimHalo_LoadPhase = (bool)tmp_int;

   fclose( File );
   if ( input_line != NULL )  free( input_line );

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadProfile
// Description :  Load the radial profiles of eta, dimensionless mass density, and phase, and then convert
//                them to the code units
//
// Note        :  1. Regard "#" as comment
//                2. Assumed format (starting from the 0th column)
//                      eta      : 0th column
//                      density  : 1th column
//                      velocity : 2th column
//                      phase    : 3th column
//                3. Phase is defined as the integration of the dimensionless velocity directly
//                   WITHOUT any unit conversion
//                4. "Profile" array must be deallocated manually
//                5. LoadPhase == true  --> load the phase of wave function from the profile table directly
//                                false --> load velocity and integrate to get phase
//
// Parameter   :  FileName : Name of the target profile table
//                Profile  : Pointer to be allocated to store the profile data
//                NBin     : Number of radial bins
//
// Return      :  Profile, NBin
//-------------------------------------------------------------------------------------------------------
void LoadProfile( const char *FileName, double **Profile, int &NBin, const bool LoadPhase )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char   *Line       = new char [MaxLine];
   char   *FirstChar  = NULL;
   FILE   *File       = NULL;
   double *Velocity  = NULL;
   double  tmp;
   int     NLoad;


// check the input file
   if (  ( File = fopen(FileName, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );


// get the total number of mass bins
   NBin = CountRow( FileName );


// allocate data
   for (int v=0; v<3; v++)    Profile[v] = new double [NBin];

   if ( !LoadPhase )          Velocity   = new double [NBin];


// loop over all rows in the input file
   NLoad = 0;
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
//    skip empty lines
      if (  !CheckEmptyString( Line )  )
      {
         FirstChar = Line;

//       find the first non-empty character
         while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//       skip lines starting with "#"
         if ( *FirstChar != '#' )
         {
            if ( NLoad >= NBin )    Aux_Error( ERROR_INFO, "NLoad (%d) >= NBin (%d) !!\n", NLoad, NBin );

            if ( LoadPhase )
               sscanf( Line, "%lf%lf%lf%lf", Profile[0]+NLoad, Profile[1]+NLoad, &tmp, Profile[2]+NLoad );

            else
               sscanf( Line, "%lf%lf%lf",    Profile[0]+NLoad, Profile[1]+NLoad, Velocity+NLoad );

            NLoad ++;
         }
      }
   } // while ( fgets( Line, MaxLine, File ) != NULL )

   fclose( File );
   delete [] Line;


// verify the number of loaded data
   if ( NLoad != NBin )    Aux_Error( ERROR_INFO, "total number of loaded data (%d) != expect (%d) !!\n", NLoad, NBin );


// integrate velocity to get phase
   if ( !LoadPhase )
   {
      Profile[2][0] = 0.0;

      for (int b=1; b<NBin; b++)
      {
         const int bm1 = b - 1;

         Profile[2][b] = Profile[2][bm1] + 0.5*( Velocity[b] + Velocity[bm1] )*( Profile[0][b] - Profile[0][bm1] );
      }

      delete [] Velocity;
   }


// convert to code units (assuming the input units are dimensionless)
   const double H0    = 100.0*Const_km/Const_s/Const_Mpc*HUBBLE0 * UNIT_T;    // Hubble parameter at z=0 in the code units
   const double Eta2x = pow( 1.5*H0*ELBDM_ETA, -0.5 )*pow( Time[0], -0.25 );  // conversion factor between eta and the comoving distance

   for (int b=0; b<NBin; b++)
   {
      Profile[0][b] *= Eta2x;       // length_in_gamer  = comoving_distance
//    Profile[1][b] *= 1.0;         // density_in_gamer = density_in_self_similar_solution
//    Profile[2][b] *= 1.0;         // phase_in_gamer   = phase_in_self_similar_solution
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LoadProfile



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckEmptyString
// Description :  Check whether the input string is empty
//
// Note        :  Empty string is defined as a string containing only " ", "\n" and "\t"
//
// Return      :  true  : Input string is empty
//                false : Input string is NOT empty
//-------------------------------------------------------------------------------------------------------
bool CheckEmptyString( const char *InputString )
{
   static const char *EmptyChar = " \n\t";

   return strspn( InputString, EmptyChar ) == strlen( InputString );

} // FUNCTION : CheckEmptyString



//-------------------------------------------------------------------------------------------------------
// Function    :  CountRow
// Description :  Count the total number of rows in the input table
//
// Parameter   :  FileName : Name of the input file
//
// Return      :  NRow
//-------------------------------------------------------------------------------------------------------
int CountRow( const char *FileName )
{

   const int MaxLine = 1024;        // maximum number of characters per line

   int   NRow      = 0;             // number of data rows
   int   NEmpty    = 0;             // number of empty rows
   char *FirstChar = NULL;
   char *Line = new char [MaxLine];

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );

// get the number of data rows
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
      if (  !CheckEmptyString( Line )  )
      {
         FirstChar = Line;

//       find the first non-empty character
         while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//       skip lines starting with "#"
         if ( *FirstChar != '#' )   NRow   ++;
         else                       NEmpty ++;
      }

      else                          NEmpty ++;
   }

   if ( NRow < 1 )   Aux_Error( ERROR_INFO, "no data rows are found !!\n" );

   fclose( File );
   delete [] Line;


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "   Input file = %s\n", FileName );
      Aux_Message( stdout, "   Data  rows = %d\n", NRow     );
      Aux_Message( stdout, "   Empty rows = %d\n", NEmpty   );
   }

   return NRow;

} // FUNCTION : CountRow



#endif // #if ( MODEL == ELBDM )
