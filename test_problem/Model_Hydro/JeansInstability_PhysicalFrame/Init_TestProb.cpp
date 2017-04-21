#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void HYDRO_TestProbSol_JeansInstability_Physical( real fluid[], const double x, const double y, const double z,
                                                         const double Time );
static void HYDRO_OutputError_JeansInstability_Physical( const bool BaseOnly );
static void WriteFile( FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk, double L1_Err[], const OptOutputPart_t Part  );
static void LoadTestProbParameter();


// global variables in the HYDRO Jeans instability test
// =======================================================================================
static real Jeans_WaveLength;    // wavelength
static real Jeans_WaveK;         // wavenumber
static real Jeans_WaveKj;        // critical wave number
static real Jeans_Omega;         // wave angular frequency
static real Jeans_RhoAmp;        // amplitude of the density perturbation (assuming background density = 1.0)
static real Jeans_Cs;            // sound speed
static real Jeans_v0;            // background velocity
static real Jeans_Sign;          // stable   : (+1/-1) --> (right/left-moving wave)
                                 // unstable : (+1/-1) --> (growing/decaying mode)
static real Jeans_Phase0;        // initial phase shift
static bool Jeans_Stable;        // true/false --> Jeans stable/unstable
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO Jeans instability test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_JeansInstability"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "HYDRO physical-frame Jeans instability";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO Jeans instability test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the HYDRO Jeans instability test !!
#  endif

#  ifndef FLOAT8
#  error : ERROR : "FLOAT8 must be ON" in the HYDRO Jeans instability test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO Jeans instability test !!
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in the %s test !!\n", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_JeansInstability_Physical;
   Output_TestProbErr_Ptr = HYDRO_OutputError_JeansInstability_Physical;


// load the test problem parameters
   LoadTestProbParameter();


// set global variables
   Jeans_WaveLength = amr->BoxSize[0]/sqrt(3.0);
   Jeans_WaveK      = 2.0*M_PI/Jeans_WaveLength;
   Jeans_WaveKj     = SQRT( (4.0*M_PI*NEWTON_G)/SQR(Jeans_Cs) );
   Jeans_Stable     = ( Jeans_WaveK > Jeans_WaveKj ) ? true : false;
   Jeans_Omega      = ( Jeans_Stable ) ? SQRT(  SQR(Jeans_Cs)*( SQR(Jeans_WaveK )-SQR(Jeans_WaveKj) )  )
                                       : SQRT(  SQR(Jeans_Cs)*( SQR(Jeans_WaveKj)-SQR(Jeans_WaveK ) )  );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : wavelength             = %14.7e\n", Jeans_WaveLength );
      Aux_Message( stdout, "       wavenumber             = %14.7e\n", Jeans_WaveK );
      Aux_Message( stdout, "       critial wavenumber     = %14.7e\n", Jeans_WaveKj );
      Aux_Message( stdout, "       stable                 = %s\n",     Jeans_Stable ? "YES" : "NO" );
      Aux_Message( stdout, "       wave angular frequency = %14.7e\n", Jeans_Omega );
      Aux_Message( stdout, "       density amplitude      = %14.7e\n", Jeans_RhoAmp );
      Aux_Message( stdout, "       sound speed            = %14.7e\n", Jeans_Cs );
      Aux_Message( stdout, "       background velocity    = %14.7e\n", Jeans_v0 );
      Aux_Message( stdout, "       sign (grow/decay;R/L)  = %14.7e\n", Jeans_Sign );
      Aux_Message( stdout, "       initial phase shift    = %14.7e\n", Jeans_Phase0 );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
// End_T : (stable/unstable) --> (1 period/grow by a factor of 50)
   const double End_T_Default    = ( Jeans_Stable) ? 2.0*M_PI/Jeans_Omega : LOG(50.0)/Jeans_Omega;
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
         Aux_Message( stdout, "NOTE : parameter %s is set to %13.7e in the %s test !!\n", "END_T", END_T, TestProb );
   }

   if ( !OPT__OUTPUT_TEST_ERROR )
   {
      OPT__OUTPUT_TEST_ERROR = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__OUTPUT_TEST_ERROR", OPT__OUTPUT_TEST_ERROR, TestProb );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_TestProbSol_JeansInstability_Physical
// Description :  Calculate the analytical solution in the HYDRO physical-frame Jeans instability test
//
// Note        :  1. Wave vector is along the diagonal direction
//                2. Background density is assumed to be ONE
//                3. This function is invoked by "HYDRO_Init_StartOver_AssignData" and "Output_TestProbErr"
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_JeansInstability_Physical( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double r         = 1.0/sqrt(3.0)*( x + y + z ) - Jeans_v0*Time;
   const double _Gamma_m1 = 1.0/(GAMMA-1.0);

   double v1, P0, P1, Phase;

   v1 = Jeans_Omega*Jeans_RhoAmp/Jeans_WaveK;
   P0 = SQR(Jeans_Cs)/GAMMA;
   P1 = SQR(Jeans_Cs)*Jeans_RhoAmp;

   if ( Jeans_Stable )
   {
      Phase       = Jeans_WaveK*r - Jeans_Sign*Jeans_Omega*Time + Jeans_Phase0;

      fluid[DENS] = 1.0 + Jeans_RhoAmp*COS( Phase );
      fluid[MOMX] = fluid[DENS]*v1    *COS( Phase )/sqrt(3.0) + fluid[DENS]*Jeans_v0/sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
      fluid[ENGY] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS]
                    + ( P0 + P1*COS(Phase) )*_Gamma_m1;
   }

   else
   {
      Phase       = Jeans_WaveK*r + Jeans_Phase0;

      fluid[DENS] =         1.0 + Jeans_RhoAmp*COS( Phase )*EXP( Jeans_Sign*Jeans_Omega*Time );
      fluid[MOMX] = -Jeans_Sign*fluid[DENS]*v1*SIN( Phase )*EXP( Jeans_Sign*Jeans_Omega*Time )/sqrt(3.0)
                    + fluid[DENS]*Jeans_v0/sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
      fluid[ENGY] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS]
                        + _Gamma_m1*(  P0 + P1*COS( Phase )*EXP( Jeans_Sign*Jeans_Omega*Time )  );
   }

} // FUNCTION : HYDRO_TestProbSol_JeansInstability_Physical



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_OutputError_JeansInstability_Physical
// Description :  Compare and output the numerical and analytical solutions in the HYDRO physical-frame
//                Jeans instability test problem
//
// Note        :  1. Invoked by "Output_TestProbErr"
//                2. This function has the similar form as Output_DumpData_Part, except that some code lines
//                   are modified to compute and output errors for this particular problem
//
// Parameter   :  BaseOnly :  Only output the base-level data
//-------------------------------------------------------------------------------------------------------
void HYDRO_OutputError_JeansInstability_Physical( const bool BaseOnly )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// set the default parameters (may need to be modified for different test problems)
// ===================================================================================================
   const OptOutputPart_t Part = OUTPUT_DIAG;
   const real            x    = NULL_REAL;
   const real            y    = NULL_REAL;
   const real            z    = NULL_REAL;


// output file name
   char FileName[NCOMP_FLUID][200];
   int  ID[6];

   ID[0] = DumpID/100000;
   ID[1] = DumpID%100000/10000;
   ID[2] = DumpID%10000/1000;
   ID[3] = DumpID%1000/100;
   ID[4] = DumpID%100/10;
   ID[5] = DumpID%10;

#  if   ( MODEL == HYDRO )
   sprintf( FileName[0], "HYDRO_JI_DENS_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[1], "HYDRO_JI_MOMX_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[2], "HYDRO_JI_MOMY_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[3], "HYDRO_JI_MOMZ_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[4], "HYDRO_JI_PRES_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   sprintf( FileName[0], "ELBDM_JI_DENS_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[1], "ELBDM_JI_REAL_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[2], "ELBDM_JI_IMAG_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
// ===================================================================================================


// check if the output files already exist
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<NCOMP_FLUID; v++)
      {
         FILE *File_Check = fopen( FileName[v], "r" );

         if ( File_Check != NULL )
         {
            Aux_Message( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n",
                         FileName[v] );
            fclose( File_Check );

            FILE *Temp = fopen( FileName[v], "w" );
            fclose( Temp );
         }
      }
   }


   const real dh_min = amr->dh[NLEVEL-1];
   const int  NLv    = ( BaseOnly ) ? 1 : NLEVEL;

   int  ii, jj, kk, scale;
   real dh, PW;
   real xx, yy, zz;        // grid physical coordinates
   int *Corner  = NULL;    // corner grid ID
   bool Check_x = false;
   bool Check_y = false;
   bool Check_z = false;

   double L1_Err[NCOMP_FLUID];
   static bool FirstTime = true;

   for (int v=0; v<NCOMP_FLUID; v++)   L1_Err[v] = 0.0;


   switch ( Part )
   {
      case OUTPUT_XY :                                      Check_z = true;   break;
      case OUTPUT_YZ :  Check_x = true;                                       break;
      case OUTPUT_XZ :                    Check_y = true;                     break;
      case OUTPUT_X  :                    Check_y = true;   Check_z = true;   break;
      case OUTPUT_Y  :  Check_x = true;                     Check_z = true;   break;
      case OUTPUT_Z  :  Check_x = true;   Check_y = true;                     break;
   }


   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File[NCOMP_FLUID];
         for (int v=0; v<NCOMP_FLUID; v++)   File[v] = fopen( FileName[v], "a" );

//       output header
         if ( TargetMPIRank == 0 )
         {
            for (int v=0; v<NCOMP_FLUID; v++)
               fprintf( File[v], "%9s %20s %20s %20s\n", "Coord.", "Numerical", "Analytical", "Error" );
         }


//       output data
         for (int lv=0; lv<NLv; lv++)
         {
            dh    = amr->dh   [lv];
            scale = amr->scale[lv];
            PW    = PATCH_SIZE*dh;

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
//             output the patch data only if it has no son (if the option "BaseOnly" is turned off)
               if ( amr->patch[0][lv][PID]->son == -1  ||  BaseOnly )
               {
                  Corner = amr->patch[0][lv][PID]->corner;

                  if ( Part == OUTPUT_DIAG ) // (+1,+1,+1) diagonal
                  {
                     if ( Corner[0] == Corner[1]  &&  Corner[0] == Corner[2] )
                     {
                        for (int k=0; k<PS1; k++)
                        {
                           kk = Corner[2] + k*scale;

                           WriteFile( File, lv, PID, k, k, k, kk, kk, kk, L1_Err, Part );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the targeted range
                     if (  !Check_x  ||  ( Corner[0]*dh_min <= x && Corner[0]*dh_min+PW > x )  )
                     if (  !Check_y  ||  ( Corner[1]*dh_min <= y && Corner[1]*dh_min+PW > y )  )
                     if (  !Check_z  ||  ( Corner[2]*dh_min <= z && Corner[2]*dh_min+PW > z )  )
                     {
//                      check whether the cell is within the targeted range
                        for (int k=0; k<PS1; k++)  {  kk = Corner[2] + k*scale;  zz = kk*dh_min;
                                                      if ( Check_z && ( zz>z || zz+dh<=z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  jj = Corner[1] + j*scale;  yy = jj*dh_min;
                                                      if ( Check_y && ( yy>y || yy+dh<=y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  ii = Corner[0] + i*scale;  xx = ii*dh_min;
                                                      if ( Check_x && ( xx>x || xx+dh<=x ) )    continue;

                           WriteFile( File, lv, PID, i, j, k, ii, jj, kk, L1_Err, Part );

                        }}}
                     } // if patch corner is within the targeted range

                  } // if ( Part == OUTPUT_DIAG ... else ... )
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

         } // for (int lv=0; lv<NLv; lv++)

         for (int v=0; v<NCOMP_FLUID; v++)   fclose( File[v] );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


// gather the L1 error from all ranks and output the results
   double L1_Err_Sum[NCOMP_FLUID];
   MPI_Reduce( L1_Err, L1_Err_Sum, NCOMP_FLUID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<NCOMP_FLUID; v++)   L1_Err_Sum[v] /= amr->BoxSize[0];

      FILE *File_L1 = fopen( "Record__L1Err", "a" );

//    output header
      if ( FirstTime )
      {
#        if   ( MODEL == HYDRO )
         fprintf( File_L1, "%5s %13s %19s %19s %19s %19s %19s\n",
                  "NGrid", "Time", "Error(DENS)", "Error(MOMX)", "Error(MOMY)", "Error(MOMZ)", "Error(PRES)" );

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )
         fprintf( File_L1, "%5s %13s %19s %19s %19s\n",
                  "NGrid", "Time", "Error(DENS)", "Error(REAL)", "Error(IMAG)" );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         FirstTime = false;
      } // if ( FirstTime )

//    output data
      fprintf( File_L1, "%5d %13.7e", NX0_TOT[0], Time[0] );

      for (int v=0; v<NCOMP_FLUID; v++)
      fprintf( File_L1, " %19.12e", L1_Err_Sum[v] );

      fprintf( File_L1, "\n" );


      fclose( File_L1 );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : HYDRO_OutputError_JeansInstability_Physical



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  Output data to file
//
// Parameter   :  File     : File pointer
//                lv       : Targeted refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                ii/jj/kk : Cell scale indices in the simulation domain
//                L1_Err   : Array to record the L1 errors of all variables
//                Part     : OUTPUT_XY   : xy plane
//                           OUTPUT_YZ   : yz plane
//                           OUTPUT_XZ   : xz plane
//                           OUTPUT_X    : x  line
//                           OUTPUT_Y    : y  line
//                           OUTPUT_Z    : z  line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//-------------------------------------------------------------------------------------------------------
void WriteFile( FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                const int ii, const int jj, const int kk, double L1_Err[], const OptOutputPart_t Part  )
{

   real fluid[NCOMP_FLUID], Anal[NCOMP_FLUID], Err[NCOMP_FLUID];

// get the numerical solution
   for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

// convert total energy to pressure
#  if   ( MODEL == HYDRO )
   fluid[ENGY] = (GAMMA-1.0)*(  fluid[ENGY] - 0.5*( fluid[MOMX]*fluid[MOMX] +
                                                    fluid[MOMY]*fluid[MOMY] +
                                                    fluid[MOMZ]*fluid[MOMZ] ) / fluid[DENS]  );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL


// get the analytical solution
   const real x = ( ii + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];
   const real y = ( jj + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];
   const real z = ( kk + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];

// ===================================================================================================
   HYDRO_TestProbSol_JeansInstability_Physical( Anal, x, y, z, Time[0] );
// ===================================================================================================

// convert total energy to pressure
#  if   ( MODEL == HYDRO )
   Anal[ENGY] = (GAMMA-1.0)*(  Anal[ENGY] - 0.5*( Anal[MOMX]*Anal[MOMX] +
                                                  Anal[MOMY]*Anal[MOMY] +
                                                  Anal[MOMZ]*Anal[MOMZ] ) / Anal[DENS]  );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // MODEL


// record the physical coordinate
   real r;

   switch ( Part )
   {
      case OUTPUT_X    :   r = x;            break;
      case OUTPUT_Y    :   r = y;            break;
      case OUTPUT_Z    :   r = z;            break;
      case OUTPUT_DIAG :   r = sqrt(3.0)*x;  break;
      default          :   Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [4/5/6/7] !!\n", Part );
   }


// estimate and output errors
   for (int v=0; v<NCOMP_FLUID; v++)
   {
      Err   [v]  = fabs( Anal[v] - fluid[v] );
      L1_Err[v] += Err[v]*amr->dh[lv];

      fprintf( File[v], "%9.7f %20.13e %20.13e %20.13e\n", r, fluid[v], Anal[v], Err[v] );
   }

} // FUNCTION : WriteFile



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter
// Description :  Load parameters for the test problem
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter()
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   int    temp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;

   getline( &input_line, &len, File );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jeans_RhoAmp,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jeans_Cs,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jeans_v0,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jeans_Sign,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jeans_Phase0,             string );

#  else

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Jeans_RhoAmp,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Jeans_Cs,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Jeans_v0,                 string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Jeans_Sign,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Jeans_Phase0,             string );
#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// force Jeans_Sign to be +1.0/-1.0
   if ( Jeans_Sign >= 0.0 )   Jeans_Sign = +1.0;
   else                       Jeans_Sign = -1.0;


} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
