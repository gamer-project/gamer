#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined PARTICLE )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_ClusterMerger( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadProfile( const char *FileName, double **Profile, int &NBin );
static bool CheckEmptyString( const char *InputString );
static int  CountRow( const char *FileName );


// global variables in the Cluster Merger test
// =======================================================================================
char    ClusterMerger_File_Prof1[1000];            // profile table of cluster 1
char    ClusterMerger_File_Prof2[1000];            // profile table of cluster 2
char    ClusterMerger_File_Par1 [1000];            // particle file of cluster 1
char    ClusterMerger_File_Par2 [1000];            // particle file of cluster 2
bool    ClusterMerger_Coll;                        // (true/false) --> test (cluster merger / single cluster)
double  ClusterMerger_Coll_D;                      // initial distance between two clusters
double  ClusterMerger_Coll_B;                      // impact parameter
double  ClusterMerger_Coll_BulkVel1;               // bulk velocity of cluster 1 (on the left  side)
double  ClusterMerger_Coll_BulkVel2;               // bulk velocity of cluster 2 (on the right side)

double *ClusterMerger_Prof1[3] = {NULL,NULL,NULL}; // radial profiles [radius/gas mass density/gas pressure] of cluster 1
double *ClusterMerger_Prof2[3] = {NULL,NULL,NULL}; // radial profiles [radius/gas mass density/gas pressure] of cluster 2
int     ClusterMerger_NBin1;                       // number of radial bins of cluster 1
int     ClusterMerger_NBin2;                       // number of radial bins of cluster 2
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the cluster merger test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char *TestProb = "cluster merger";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : only support the HYDRO model !!
#  endif

#  ifndef PARTICLE
#  error : ERROR : "PARTICLE must be ON" in the cluster merger test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the cluster merger test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the cluster merger test !!
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "please turn on \"OPT__UNIT\" and set units properly for this test problem !!\n" );

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "one should not adopt periodic boundary condition for this test problem !!\n" );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_ClusterMerger;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// load the radial profile
   if ( OPT__INIT != INIT_RESTART )
   {
      LoadProfile( ClusterMerger_File_Prof1, ClusterMerger_Prof1, ClusterMerger_NBin1 );

      if ( ClusterMerger_Coll )
      LoadProfile( ClusterMerger_File_Prof2, ClusterMerger_Prof2, ClusterMerger_NBin2 );
   }


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "   profile file 1   = %s\n",           ClusterMerger_File_Prof1 );
      if ( ClusterMerger_Coll )
      Aux_Message( stdout, "   profile file 2   = %s\n",           ClusterMerger_File_Prof2 );
      Aux_Message( stdout, "   particle file 1  = %s\n",           ClusterMerger_File_Par1 );
      if ( ClusterMerger_Coll )
      Aux_Message( stdout, "   particle file 2  = %s\n",           ClusterMerger_File_Par2 );
      Aux_Message( stdout, "   test mode        = %s\n",          (ClusterMerger_Coll)? "merging cluster":"single cluster" );
      if ( ClusterMerger_Coll ) {
      Aux_Message( stdout, "   initial distance = %+14.7e kpc\n",  ClusterMerger_Coll_D*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   impact parameter = %+14.7e kpc\n",  ClusterMerger_Coll_B*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   bulk velocity 1  = %+14.7e km/s\n", ClusterMerger_Coll_BulkVel1*UNIT_V/(Const_km/Const_s) );
      Aux_Message( stdout, "   bulk velocity 2  = %+14.7e km/s\n", ClusterMerger_Coll_BulkVel2*UNIT_V/(Const_km/Const_s) ); }
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;
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
// Function    :  Par_TestProbSol_ClusterMerger
// Description :  Initialize the background density field as zero for the cluster merger test
//
// Note        :  1. Isothermal beta model for gas
//                   --> Currently beta is fixed to 1
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_ClusterMerger( real *fluid, const double x, const double y, const double z, const double Time )
{

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   if ( ClusterMerger_Coll )
   {
      const double ClusterCenter1[3] = { BoxCenter[0]-0.5*ClusterMerger_Coll_D, BoxCenter[1]-0.5*ClusterMerger_Coll_B, BoxCenter[2] };
      const double ClusterCenter2[3] = { BoxCenter[0]+0.5*ClusterMerger_Coll_D, BoxCenter[1]+0.5*ClusterMerger_Coll_B, BoxCenter[2] };

      double r1, r2, Dens1, Dens2, Pres1, Pres2, Vel;

//    for each cell, we sum up the density and pressure from each halos and then calculate the weighted velocity
      r1    = sqrt( SQR(x-ClusterCenter1[0]) + SQR(y-ClusterCenter1[1]) + SQR(z-ClusterCenter1[2]) );
      r2    = sqrt( SQR(x-ClusterCenter2[0]) + SQR(y-ClusterCenter2[1]) + SQR(z-ClusterCenter2[2]) );
      Dens1 = Mis_InterpolateFromTable( ClusterMerger_NBin1, ClusterMerger_Prof1[0], ClusterMerger_Prof1[1], r1 );
      Dens2 = Mis_InterpolateFromTable( ClusterMerger_NBin2, ClusterMerger_Prof2[0], ClusterMerger_Prof2[1], r2 );
      Pres1 = Mis_InterpolateFromTable( ClusterMerger_NBin1, ClusterMerger_Prof1[0], ClusterMerger_Prof1[2], r1 );
      Pres2 = Mis_InterpolateFromTable( ClusterMerger_NBin2, ClusterMerger_Prof2[0], ClusterMerger_Prof2[2], r2 );
      Vel   = ( ClusterMerger_Coll_BulkVel1*Dens1 + ClusterMerger_Coll_BulkVel2*Dens2 ) / ( Dens1 + Dens2 );

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
      Dens = Mis_InterpolateFromTable( ClusterMerger_NBin1, ClusterMerger_Prof1[0], ClusterMerger_Prof1[1], r );
      Pres = Mis_InterpolateFromTable( ClusterMerger_NBin1, ClusterMerger_Prof1[0], ClusterMerger_Prof1[2], r );

      if ( Dens == NULL_REAL  ||  Pres == NULL_REAL )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (probably outside the input table)!!\n", r );

      fluid[DENS] = Dens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Pres / ( GAMMA - 1.0 );
   } // if ( ClusterMerger_Coll ) ... else ...

} // FUNCTION : Par_TestProbSol_ClusterMerger



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
   sscanf( input_line, "%s%s",    ClusterMerger_File_Prof1,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    ClusterMerger_File_Prof2,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    ClusterMerger_File_Par1,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    ClusterMerger_File_Par2,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &tmp_int,                     string );
   ClusterMerger_Coll = (bool)tmp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Coll_D,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Coll_B,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Coll_BulkVel1, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Coll_BulkVel2, string );


// check
   if ( ClusterMerger_Coll_D < 0.0 )   Aux_Error( ERROR_INFO, "ClusterMerger_Coll_D [%14.7e] < 0.0 !!\n", ClusterMerger_Coll_D );
   if ( ClusterMerger_Coll_B < 0.0 )   Aux_Error( ERROR_INFO, "ClusterMerger_Coll_B [%14.7e] < 0.0 !!\n", ClusterMerger_Coll_B );


// convert to code units
   ClusterMerger_Coll_D        *= Const_kpc          / UNIT_L;
   ClusterMerger_Coll_B        *= Const_kpc          / UNIT_L;
   ClusterMerger_Coll_BulkVel1 *= (Const_km/Const_s) / UNIT_V;
   ClusterMerger_Coll_BulkVel2 *= (Const_km/Const_s) / UNIT_V;

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadProfile
// Description :  Load the radial profiles of radius, gas mass density, and pressure
//
// Note        :  1. Regard "#" as comment
//                2. Assumed format (starting from the 0th column)
//                      radius   : 7th column
//                      density  : 2th column
//                      pressure : 6th column
//                3. Assumed input units: cgs
//
// Parameter   :  FileName : Name of the target profile table
//                Profile  : Pointer to be allocated to store the profile data
//                NBin     : Number of radial bins
//
// Return      :  Profile, NBin
//-------------------------------------------------------------------------------------------------------
void LoadProfile( const char *FileName, double **Profile, int &NBin )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char  *Line       = new char [MaxLine];
   char  *FirstChar  = NULL;
   FILE  *File       = NULL;
   double tmp;
   int    NLoad;


// check the input file
   if (  ( File = fopen(FileName, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );


// get the total number of mass bins
   NBin = CountRow( FileName );


// allocate data
   for (int v=0; v<3; v++)    Profile[v] = new double [NBin];


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

            sscanf( Line, "%lf%lf%lf%lf%lf%lf%lf%lf",
                    &tmp, &tmp, Profile[1]+NLoad, &tmp, &tmp, &tmp, Profile[2]+NLoad, Profile[0]+NLoad );

            NLoad ++;
         }
      }
   } // while ( fgets( Line, MaxLine, File ) != NULL )

   fclose( File );
   delete [] Line;


// verify the number of loaded data
   if ( NLoad != NBin )    Aux_Error( ERROR_INFO, "total number of loaded data (%d) != expect (%d) !!\n", NLoad, NBin );


// convert to code units (assuming the input units are cgs)
   for (int b=0; b<NBin; b++)
   {
      Profile[0][b] /= UNIT_L;
      Profile[1][b] /= UNIT_D;
      Profile[2][b] /= UNIT_P;
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



#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )
