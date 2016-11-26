#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined PARTICLE )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_ClusterMerger( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadProfile();
static bool CheckEmptyString( const char *InputString );
static int  CountRow( const char *FileName );


// global variables in the Cluster Merger test
// =======================================================================================
bool    ClusterMerger_Coll;                        // (true/false) --> test (cluster merger / single cluster)
char    ClusterMerger_File_Prof[1000];             // name of the input profile table
char    ClusterMerger_File_Par [1000];             // name of the input particle file

double *ClusterMerger_Prof[3] = {NULL,NULL,NULL};  // radial profiles of various quantities [radius/gas mass density/gas pressure]
int     ClusterMerger_NBin;                        // number of radial bins

/*
double ClusterMerger_Coll_D;           // distance between two ClusterMerger clouds for the ClusterMerger collision test
double ClusterMerger_Coll_ImpactPara;  // impact parameter
double ClusterMerger_Coll_BulkVel[3];  // bulk velocity
*/
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
   if ( OPT__INIT != INIT_RESTART )    LoadProfile();


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "   test mode     = %s\n", (ClusterMerger_Coll)? "merging cluster":"single cluster" );
      Aux_Message( stdout, "   profile file  = %s\n", ClusterMerger_File_Prof );
      Aux_Message( stdout, "   particle file = %s\n", ClusterMerger_File_Par );

      /*
      if ( ClusterMerger_Coll )
      Aux_Message( stdout, "       initial distance between two clouds          = %13.7e\n",  ClusterMerger_Coll_D );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "       bulk velocity [%d]                            = %14.7e\n", d, ClusterMerger_BulkVel[d] );
      */
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
      Aux_Error( ERROR_INFO, "NOT SUPPORETD YET !!\n" );
   }

   else
   {
      double r, Dens, Pres;

      r = sqrt( SQR(x-BoxCenter[0]) + SQR(y-BoxCenter[1]) + SQR(z-BoxCenter[2]) );

      Dens = Mis_InterpolateFromTable( ClusterMerger_NBin, ClusterMerger_Prof[0], ClusterMerger_Prof[1], r );
      Pres = Mis_InterpolateFromTable( ClusterMerger_NBin, ClusterMerger_Prof[0], ClusterMerger_Prof[2], r );

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
   sscanf( input_line, "%d%s",   &tmp_int,                  string );
   ClusterMerger_Coll = (bool)tmp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    ClusterMerger_File_Prof,  string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%s%s",    ClusterMerger_File_Par,   string );


// check

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
//-------------------------------------------------------------------------------------------------------
void LoadProfile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char  *Line       = new char [MaxLine];
   char  *FirstChar  = NULL;
   FILE  *File       = NULL;
   double tmp;
   int    NBin, NLoad;


// check the input file
   if (  ( File = fopen(ClusterMerger_File_Prof, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", ClusterMerger_File_Prof );


// get the total number of mass bins
   ClusterMerger_NBin = CountRow( ClusterMerger_File_Prof );
   NBin               = ClusterMerger_NBin;


// allocate data
   for (int v=0; v<3; v++)    ClusterMerger_Prof[v] = new double [NBin];


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
                    &tmp, &tmp, ClusterMerger_Prof[1]+NLoad, &tmp, &tmp, &tmp,
                    ClusterMerger_Prof[2]+NLoad, ClusterMerger_Prof[0]+NLoad );

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
      ClusterMerger_Prof[0][b] /= UNIT_L;
      ClusterMerger_Prof[1][b] /= UNIT_D;
      ClusterMerger_Prof[2][b] /= UNIT_P;
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
