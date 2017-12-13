#include "GAMER.h"

#if ( MODEL == ELBDM )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter( real &Scale_Temp );
static void LoadSolitonDensProfile();
static void ELBDM_TestProbSol_Soliton( real fluid[], const double x, const double y, const double z, const double Time );


// global variables in the ELBDM soliton test
// =======================================================================================
int    Soliton_N;                      // total number of solitons
int    Soliton_RSeed;                  // random seed for setting Soliton_Center (<=0 --> set manually)
int    Soliton_DensProfileN;           // array size of Soliton_DensProfileR and Soliton_DensProfile
real   Soliton_EmptyRegion;            // region where no soliton will reside at (useful only when Soliton_RSeed > 0)
real  *Soliton_DensProfileR = NULL;    // radial coord. of the soliton density profile
real  *Soliton_DensProfile  = NULL;    // soliton density profile
real  *Soliton_Scale        = NULL;    // scale factor of each soliton (<=0 --> set manually)
real (*Soliton_Center)[3]   = NULL;    // center coordinates of each soliton
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the ELBDM soliton test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{  

   const char *TestProb = "ELBDM soliton";

// check
#  if ( MODEL != ELBDM )
#  error : ERROR : "MODEL != ELBDM" in the ELBDM soliton test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the ELBDM soliton test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the ELBDM soliton test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = ELBDM_TestProbSol_Soliton;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   real Scale_Temp;
   LoadTestProbParameter( Scale_Temp );


// set the test problem parameters
   Soliton_Scale  = new real [Soliton_N];
   Soliton_Center = new real [Soliton_N][3];

   if ( Scale_Temp > (real)0.0 )
   {
      for (int t=0; t<Soliton_N; t++)  Soliton_Scale[t] = Scale_Temp;
   }

   else
   {
      if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : please set the scale factor of each soliton properly in \"file <%s>, function <%s>\"\n",
                   __FILE__, __FUNCTION__ );

      /*
      for (int t=0; t<Soliton_N; t++)  Soliton_Scale[t] = XXX;
      */
   } // if ( Scale_Temp > (real)0.0 ) ... else ...

   if ( Soliton_RSeed > 0 )
   {
      const real Coord_Min[3] = { Soliton_EmptyRegion, Soliton_EmptyRegion, Soliton_EmptyRegion };
      const real Coord_Max[3] = { amr->BoxSize[0] - Soliton_EmptyRegion,
                                  amr->BoxSize[1] - Soliton_EmptyRegion,
                                  amr->BoxSize[2] - Soliton_EmptyRegion };
      srand( Soliton_RSeed );

      for (int t=0; t<Soliton_N; t++)  
      for (int d=0; d<3; d++) 
         Soliton_Center[t][d] = ( (real)rand()/RAND_MAX )*(Coord_Max[d]-Coord_Min[d]) + Coord_Min[d];
   }

   else
   {
      if ( Soliton_N == 1 )
      {
         for (int d=0; d<3; d++)    Soliton_Center[0][d] = 0.5*amr->BoxSize[d];
      }

      else
      {
//       please comment out the following lines when the soliton positions are properly set
         if ( MPI_Rank == 0 )
            Aux_Message( stderr, "ERROR : please set the center coord. of each soliton in \"file <%s>, function <%s>\"\n",
                         __FILE__, __FUNCTION__ );

         End_GAMER();
      }

//    set the soliton positions here
      /*
      for (int t=0; t<Soliton_N; t++)  
      for (int d=0; d<3; d++)
         Soliton_Center[t][d] = XXX;
         */

   } // if ( Soliton_RSeed > 0 ) ... else ...


// load the soliton density profile
   LoadSolitonDensProfile();


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : total number of solitons                     = %d\n",     Soliton_N );
      Aux_Message( stdout, "       random seed for setting the center coord.    = %d\n",     Soliton_RSeed );
      Aux_Message( stdout, "       size of the empty zone                       = %13.7e\n", Soliton_EmptyRegion );
      Aux_Message( stdout, "       number of data points in the density profile = %d\n",     Soliton_DensProfileN );
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "       %7s  %13s  %13s  %13s  %13s\n", "ID", "Scale", "Center_X", "Center_Y", "Center_Z" );
      for (int t=0; t<Soliton_N; t++)
      Aux_Message( stdout, "       %7d  %13.6e  %13.6e  %13.6e  %13.6e\n", t, Soliton_Scale[t], Soliton_Center[t][0],
                                                                           Soliton_Center[t][1], Soliton_Center[t][2] ); 
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 1.0e3;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %ld in the %s test !!\n", "END_STEP", END_STEP, TestProb );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %13.7e in the %s test !!\n",
                      "END_T", END_T, TestProb );
   }

   if ( !OPT__INIT_RESTRICT )
   {
      OPT__INIT_RESTRICT = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__INIT_RESTRICT", OPT__INIT_RESTRICT, TestProb );
   }

   if ( OPT__OUTPUT_TEST_ERROR  &&  Output_TestProbErr_Ptr == NULL )
   {
      OPT__OUTPUT_TEST_ERROR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : option %s is not supported for the %s test and has been disabled !!\n",
                      "OPT__OUTPUT_TEST_ERROR", TestProb );
   }
   else if ( !OPT__OUTPUT_TEST_ERROR  &&  Output_TestProbErr_Ptr != NULL )
   {
      OPT__OUTPUT_TEST_ERROR = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : option %s is enabled for the %s test !!\n", 
                      "OPT__OUTPUT_TEST_ERROR", TestProb );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_TestProbSol_Soliton
// Description :  Initialize the ELBDM soliton test  
//
// Note        :  Invoked by "ELBDM_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void ELBDM_TestProbSol_Soliton( real *fluid, const double x, const double y, const double z, const double Time )
{

   const int MinIdx = 0;
   const int MaxIdx = Soliton_DensProfileN - 1;

   int    Idx;
   double ShellWidth, dr, r, Dens;


// initialize density as zero since there may be multiple solitons
   fluid[DENS] = 0.0;

// loop over all solitons to get the total density
   for (int t=0; t<Soliton_N; t++)
   {
      r   = sqrt( SQR(x-Soliton_Center[t][0]) + SQR(y-Soliton_Center[t][1]) + SQR(z-Soliton_Center[t][2]) );

//    rescale the radial coord. which ~ scale^-1 (we rescale "r" instead of "Soliton_DensProfile" just for convenience)
      r  *= Soliton_Scale[t];

//    find the proper array index in the density profile array
      Idx = Mis_BinarySearch_Real( Soliton_DensProfileR, MinIdx, MaxIdx, r );

//    linear interpolation
      if      ( Idx <  MinIdx )  Dens = Soliton_DensProfile[MinIdx];
      else if ( Idx >= MaxIdx )  Dens = Soliton_DensProfile[MaxIdx];
      else
      {
         ShellWidth = Soliton_DensProfileR[Idx+1] - Soliton_DensProfileR[Idx];
         dr         = r - Soliton_DensProfileR[Idx];
         Dens       = ( dr*Soliton_DensProfile[Idx+1] + (ShellWidth-dr)*Soliton_DensProfile[Idx] ) / ShellWidth;

#        ifdef GAMER_DEBUG
         if ( ShellWidth <= 0.0 )   Aux_Error( ERROR_INFO, "ShellWidth (%14.7e) <= 0.0 !!\n", ShellWidth );
         if ( dr < 0.0 )            Aux_Error( ERROR_INFO, "dr (%14.7e) < 0.0 !!\n", dr );
         if ( dr > ShellWidth )     Aux_Error( ERROR_INFO, "dr (%14.7e) > ShellWidth (%14.7e) !!\n", dr, ShellWidth );
#        endif
      }

#     ifdef GAMER_DEBUG
      if ( Dens < 0.0 )    
         Aux_Error( ERROR_INFO, "Dens (%14.7e) < 0.0 (dr %14.7e, ShellWidth %14.7e, DensL %14.7e, DensR %14.7e) !!\n",
                    Dens, dr, ShellWidth, Soliton_DensProfile[Idx], Soliton_DensProfile[Idx+1] );
#     endif

//    rescale the soliton density (~scale^4) and add to the fluid array
      fluid[DENS] += Dens*SQR( Soliton_Scale[t] )*SQR( Soliton_Scale[t] );

   } // for (int t=0; t<Soliton_N; t++)


// set the real and imaginary parts
   fluid[REAL] = SQRT( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : ELBDM_TestProbSol_Soliton



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter 
// Description :  Load parameters for the test problem 
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  Scale_Temp  : Temporary soliton scale factor
//
// Return      :  Scale_Temp
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter( real &Scale_Temp )
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;


   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Soliton_N,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Soliton_RSeed,            string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Scale_Temp,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Soliton_EmptyRegion,      string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Scale_Temp,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Soliton_EmptyRegion,      string );
#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check
   if ( Soliton_N <= 0 )   Aux_Error( ERROR_INFO, "Soliton_N (%d) <= 0 !!\n", Soliton_N );

   if ( Soliton_RSeed > 0  &&  Soliton_EmptyRegion < 0.0 )  Aux_Error( ERROR_INFO, "Soliton_EmptyRegion (%14.7e) < 0.0 !!\n",
                                                                       Soliton_EmptyRegion );

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadSolitonDensProfile
// Description :  Load the soliton density profile
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//
// Return      :  Soliton_DensProfileN, Soliton_DensProfileR, Soliton_DensProfile
//-------------------------------------------------------------------------------------------------------
void LoadSolitonDensProfile()
{

   const char FileName[] = "Input__SolitonDensProfile";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;


// get the number of lines
   Soliton_DensProfileN = 0;
   getline( &input_line, &len, File );    // skip the header
   while ( getline( &input_line, &len, File ) != -1 )    Soliton_DensProfileN ++;

   if ( Soliton_DensProfileN <= 0 )
         Aux_Error( ERROR_INFO, "incorrect reading of the file \"%s\" (number of data points == %d <= 0) !!\n", 
                    FileName, Soliton_DensProfileN );


// load the density profile
   Soliton_DensProfileR = new real [Soliton_DensProfileN];
   Soliton_DensProfile  = new real [Soliton_DensProfileN];

   fseek( File, 0, SEEK_SET );
   getline( &input_line, &len, File );    // skip the header

   for (int t=0; t<Soliton_DensProfileN; t++)
   {
      if (  getline( &input_line, &len, File ) == -1  )  
         Aux_Error( ERROR_INFO, "incorrect reading of the file \"%s\" !!\n", FileName );

#     ifdef FLOAT8
      sscanf( input_line, "%lf%lf", Soliton_DensProfileR+t, Soliton_DensProfile+t );
#     else
      sscanf( input_line, "%f%f",   Soliton_DensProfileR+t, Soliton_DensProfile+t );
#     endif
   }

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check the monotonicity of the radial coord.
   for (int t=1; t<Soliton_DensProfileN; t++)
   {
      if ( Soliton_DensProfileR[t] <= Soliton_DensProfileR[t-1] )
         Aux_Error( ERROR_INFO, "radial coord. of the density profile is not monotonic (R[%d] = %14.7e, R[%d] = %14.7e) !!\n",
                    t-1, Soliton_DensProfileR[t-1], t, Soliton_DensProfileR[t] );
   }

} // FUNCTION : LoadSolitonDensProfile



#endif // #if ( MODEL == ELBDM )
