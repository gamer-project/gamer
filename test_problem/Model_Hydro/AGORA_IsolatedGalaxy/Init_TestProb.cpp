#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined PARTICLE )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void HYDRO_TestProbSol_AGORA( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadVcProf( const char *Filename, double **Profile, int &NBin );
static bool CheckEmptyString( const char *InputString );
static int  CountRow( const char *Filename );


// global variables in the AGORA isolated galaxy test
// =======================================================================================
char    AGORA_VcProf_Filename[1000];      // filename of the circular velocity radial profile
double  AGORA_DiskScaleLength;            // disk scale length
double  AGORA_DiskScaleHeight;            // disk scale height
double  AGORA_DiskTotalMass;              // disk total mass (gas + stars)
double  AGORA_DiskGasMassFrac;            // disk gas mass fraction (disk_gas_mass / disk_total_mass)
double  AGORA_DiskGasTemp;                // disk gas temperature
double  AGORA_HaloGasNumDensH;            // halo atomic hydrogen number density (halo gas mass density / atomic_hydrogen_mass)
double  AGORA_HaloGasTemp;                // halo gas temperature

double  AGORA_DiskGasDens0;               // disk gas mass density coefficient
double  AGORA_HaloGasDens;                // halo gas mass density
double  AGORA_HaloGasPres;                // halo gas pressure
double *AGORA_VcProf[2] = {NULL,NULL};    // circular velocity radial profile [radius, velocity]
int     AGORA_VcProf_NBin;                // number of radial bin in AGORA_VcProf

/*
AgoraRestartHaloMetallicity              = 1e-6
AgoraRestartDiskMetallicity              = 1.0
AgoraRestartMetalFractionByMass          = 0.02041
*/
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the AGORA isolated galaxy test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char *TestProb = "AGORA isolated galaxy";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : only support the HYDRO model !!
#  endif

#  ifndef PARTICLE
#  error : ERROR : "PARTICLE must be ON" in the AGORA isolated galaxy test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the AGORA isolated galaxy test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the AGORA isolated galaxy test !!
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "please turn on \"OPT__UNIT\" and set units properly for this test problem !!\n" );

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "one should not adopt periodic boundary condition for this test problem !!\n" );


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_AGORA;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set global variables
   const bool CheckMinPres_Yes = true;

   AGORA_DiskGasDens0 = AGORA_DiskTotalMass*AGORA_DiskGasMassFrac / ( 4.0*M_PI*SQR(AGORA_DiskScaleLength)*AGORA_DiskScaleHeight );
   AGORA_HaloGasDens  = AGORA_HaloGasNumDensH*Const_mH/UNIT_M;
   AGORA_HaloGasPres  = CPU_Temperature2Pressure( AGORA_HaloGasDens, AGORA_HaloGasTemp, MOLECULAR_WEIGHT, Const_mH/UNIT_M,
                                                  CheckMinPres_Yes, MIN_PRES );


// load the circular velocity radial profile
   if ( OPT__INIT != INIT_RESTART )    LoadVcProf( AGORA_VcProf_Filename, AGORA_VcProf, AGORA_VcProf_NBin );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "   VcProf_Filename = %s\n",             AGORA_VcProf_Filename                     );
      Aux_Message( stdout, "   DiskScaleLength = %13.7e kpc\n",     AGORA_DiskScaleLength * UNIT_L/Const_kpc  );
      Aux_Message( stdout, "   DiskScaleHeight = %13.7e kpc\n",     AGORA_DiskScaleHeight * UNIT_L/Const_kpc  );
      Aux_Message( stdout, "   DiskTotalMass   = %13.7e Msun\n",    AGORA_DiskTotalMass   * UNIT_M/Const_Msun );
      Aux_Message( stdout, "   DiskGasMassFrac = %13.7e\n",         AGORA_DiskGasMassFrac                     );
      Aux_Message( stdout, "   DiskGasTemp     = %13.7e K\n",       AGORA_DiskGasTemp     * UNIT_E/Const_kB   );
      Aux_Message( stdout, "   HaloGasNumDensH = %13.7e cm^{-3}\n", AGORA_HaloGasNumDensH / CUBE(UNIT_L)      );
      Aux_Message( stdout, "   HaloGasTemp     = %13.7e K\n",       AGORA_HaloGasTemp     * UNIT_E/Const_kB   );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 500.0*Const_Myr/UNIT_T;
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

   if ( OPT__OUTPUT_TEST_ERROR )
   {
      OPT__OUTPUT_TEST_ERROR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__OUTPUT_TEST_ERROR", OPT__OUTPUT_TEST_ERROR, TestProb );
   }

   if ( OPT__INIT == INIT_STARTOVER  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = PAR_INIT_BY_FUNCTION !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_TestProbSol_AGORA
// Description :  Initialize the gas disk and halo for the AGORA isolated galaxy test
//
// Note        :  1. We do NOT truncate gas disk. Instead, we ensure pressure balance between gas disk and halo.
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_AGORA( real *fluid, const double x, const double y, const double z, const double Time )
{

   const bool   CheckMinPres_Yes = true;
   const double dx               = x - 0.5*amr->BoxSize[0];
   const double dy               = y - 0.5*amr->BoxSize[1];
   const double dz               = z - 0.5*amr->BoxSize[2];
   const double r                = sqrt( SQR(dx) + SQR(dy) );
   const double h                = fabs( dz );

   double DiskGasDens, DiskGasPres, DiskGasVel;

   DiskGasDens = AGORA_DiskGasDens0 * exp( -r/AGORA_DiskScaleLength ) * exp( -h/AGORA_DiskScaleHeight );
   DiskGasPres = CPU_Temperature2Pressure( DiskGasDens, AGORA_DiskGasTemp, MOLECULAR_WEIGHT, Const_mH/UNIT_M,
                                           CheckMinPres_Yes, MIN_PRES );

// disk component
   if ( DiskGasPres > AGORA_HaloGasPres )
   {
//    perform linear interpolation to get the gas circular velocity at r
      if (  ( DiskGasVel = Mis_InterpolateFromTable(AGORA_VcProf_NBin, AGORA_VcProf[0], AGORA_VcProf[1], r) ) == NULL_REAL  )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e !!\n", r );

      fluid[DENS] = DiskGasDens;
      fluid[MOMX] = -dy/r*DiskGasVel*fluid[DENS];
      fluid[MOMY] = +dx/r*DiskGasVel*fluid[DENS];
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = DiskGasPres / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   }

// halo component
   else
   {
      fluid[DENS] = AGORA_HaloGasDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = AGORA_HaloGasPres / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   } // if ( DiskPres > AGORA_HaloGasPres ) ... else ...

} // FUNCTION : HYDRO_TestProbSol_AGORA



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

   const char Filename[] = "Input__TestProb";

   FILE *File = fopen( Filename, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", Filename );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;
   int    tmp_int;


// skip the header
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    AGORA_VcProf_Filename,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_DiskScaleLength,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_DiskScaleHeight,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_DiskTotalMass,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_DiskGasMassFrac,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_DiskGasTemp,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_HaloGasNumDensH,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &AGORA_HaloGasTemp,           string );


// check
   if ( AGORA_DiskScaleLength <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_DiskScaleLength [%14.7e] <= 0.0 !!\n", AGORA_DiskScaleLength );
   if ( AGORA_DiskScaleHeight <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_DiskScaleHeight [%14.7e] <= 0.0 !!\n", AGORA_DiskScaleHeight );
   if ( AGORA_DiskTotalMass   <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_DiskTotalMass   [%14.7e] <= 0.0 !!\n", AGORA_DiskTotalMass   );
   if ( AGORA_DiskGasMassFrac <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_DiskGasMassFrac [%14.7e] <= 0.0 !!\n", AGORA_DiskGasMassFrac );
   if ( AGORA_DiskGasTemp     <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_DiskGasTemp     [%14.7e] <= 0.0 !!\n", AGORA_DiskGasTemp     );
   if ( AGORA_HaloGasNumDensH <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_HaloGasNumDensH [%14.7e] <= 0.0 !!\n", AGORA_HaloGasNumDensH );
   if ( AGORA_HaloGasTemp     <= 0.0 )    Aux_Error( ERROR_INFO, "AGORA_HaloGasTemp     [%14.7e] <= 0.0 !!\n", AGORA_HaloGasTemp     );


// convert to code units
   AGORA_DiskScaleLength *= Const_kpc  / UNIT_L;
   AGORA_DiskScaleHeight *= Const_kpc  / UNIT_L;
   AGORA_DiskTotalMass   *= Const_Msun / UNIT_M;
   AGORA_DiskGasTemp     *= Const_kB   / UNIT_E;
   AGORA_HaloGasNumDensH *= 1.0        * CUBE(UNIT_L);
   AGORA_HaloGasTemp     *= Const_kB   / UNIT_E;

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadVcProf
// Description :  Load the circular velocity radial profile
//
// Note        :  1. Regard "#" as comment
//                2. Assumed format [radius in kpc] [velocity in km/s]
//                3. This function will allocate memory for the input pointer "Profile"
//                   --> One must deallocate it manually
//
// Parameter   :  Filename : Name of the target profile table
//                Profile  : Pointer to be allocated to store the profile data
//                NBin     : Number of radial bins
//
// Return      :  Profile, NBin
//-------------------------------------------------------------------------------------------------------
void LoadVcProf( const char *Filename, double **Profile, int &NBin )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char  *Line       = new char [MaxLine];
   char  *FirstChar  = NULL;
   FILE  *File       = NULL;
   int    NLoad;


// check the input file
   if (  ( File = fopen(Filename, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", Filename );


// get the total number of radial bins
   NBin = CountRow( Filename );


// allocate data
   for (int v=0; v<2; v++)    Profile[v] = new double [NBin];


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

            sscanf( Line, "%lf%lf", Profile[0]+NLoad, Profile[1]+NLoad );

            NLoad ++;
         }
      }
   } // while ( fgets( Line, MaxLine, File ) != NULL )

   fclose( File );
   delete [] Line;


// verify the number of loaded data
   if ( NLoad != NBin )    Aux_Error( ERROR_INFO, "total number of loaded data (%d) != expect (%d) !!\n", NLoad, NBin );


// convert to code units
   for (int b=0; b<NBin; b++)
   {
      Profile[0][b] *= Const_kpc / UNIT_L;
      Profile[1][b] *= Const_km  / UNIT_V;
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LoadVcProf



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
// Parameter   :  Filename : Name of the input file
//
// Return      :  NRow
//-------------------------------------------------------------------------------------------------------
int CountRow( const char *Filename )
{

   const int MaxLine = 1024;        // maximum number of characters per line

   int   NRow      = 0;             // number of data rows
   int   NEmpty    = 0;             // number of empty rows
   char *FirstChar = NULL;
   char *Line = new char [MaxLine];

   FILE *File = fopen( Filename, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", Filename );

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
      Aux_Message( stdout, "   Input file = %s\n", Filename );
      Aux_Message( stdout, "   Data  rows = %d\n", NRow     );
      Aux_Message( stdout, "   Empty rows = %d\n", NEmpty   );
   }

   return NRow;

} // FUNCTION : CountRow



#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )
