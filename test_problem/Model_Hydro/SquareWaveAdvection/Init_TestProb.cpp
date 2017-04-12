#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void HYDRO_TestProbSol_SquareAdv( real fluid[], const double x, const double y, const double z,
                                                         const double Time );
static void HYDRO_OutputError_SquareAdv( const bool BaseOnly );
static void WriteFile( FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk, double L1_Err[], const OptOutputPart_t Part  );
static void LoadTestProbParameter();


// global variables in the HYDRO square wave advection test
// =======================================================================================
static double SquareAdv_RhoPeak;       // peak density (assuming background density = 1.0)
static double SquareAdv_Pres;          // background pressure
static double SquareAdv_Vel;           // background velocity
static double SquareAdv_Width;         // width of the square wave (default = 0.5*box_size)
static bool   SquareAdv_Smooth;        // smooth out the discontinuity in the square wave
static double SquareAdv_Smooth_Eps;    // width of the transition region used by SquareAdv_Smooth
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO square wave advection test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_SquareAdv"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "HYDRO square wave advection";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO square wave advection test !!
#  endif

#  ifdef GRAVITY
#  error : ERROR : "GRAVITY must be OFF" in the HYDRO square wave advection test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO square wave advection test !!
#  endif

#  ifdef PARTICLE
#  error : ERROR : "PARTICLE must be OFF" in the HYDRO square wave advection test !!
#  endif

   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "Please set \"OPT__BC_FLU = 1\" for the %s test!!\n", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_SquareAdv;
   Output_TestProbErr_Ptr = HYDRO_OutputError_SquareAdv;


// load the test problem parameters
   LoadTestProbParameter();


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : peak density          = % 14.7e\n", SquareAdv_RhoPeak    );
      Aux_Message( stdout, "       background pressure   = % 14.7e\n", SquareAdv_Pres       );
      Aux_Message( stdout, "       background velocity   = % 14.7e\n", SquareAdv_Vel        );
      Aux_Message( stdout, "       square wave width     = % 14.7e\n", SquareAdv_Width      );
      Aux_Message( stdout, "       smooth the transition = %d\n",      SquareAdv_Smooth     );
      if ( SquareAdv_Smooth )
      Aux_Message( stdout, "       transition widith     = % 14.7e\n", SquareAdv_Smooth_Eps );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 1.0;
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
// Function    :  HYDRO_TestProbSol_SquareAdv
// Description :  Calculate the analytical solution in the HYDRO square wave advection test
//
// Note        :  1. Only support 1D along the x direction
//                2. Background density is assumed to be ONE
//                3. This function is invoked by "HYDRO_Init_StartOver_AssignData" and "Output_TestProbErr"
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_SquareAdv( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double L         = amr->BoxSize[0];
   const double Cen       = fmod( 0.5*L + SquareAdv_Vel*Time, L );
   const double dx        = fmin(  fmin( fabs(x-Cen), fabs(x-Cen-L) ), fabs(x-Cen+L)  );
   const double _Gamma_m1 = 1.0/(GAMMA-1.0);
   const double RhoBg     = 1.0;
   const double HalfWidth = 0.5*SquareAdv_Width;

   if ( SquareAdv_Smooth )
   {
      fluid[DENS] = RhoBg + 0.5*( SquareAdv_RhoPeak-RhoBg )*
                                (  tanh( (dx+HalfWidth)/SquareAdv_Smooth_Eps ) + tanh( (-dx+HalfWidth)/SquareAdv_Smooth_Eps )  );
   }

   else
   {
      if ( dx < HalfWidth )   fluid[DENS] = SquareAdv_RhoPeak;
      else                    fluid[DENS] = RhoBg;
   }

   fluid[MOMX] = fluid[DENS]*SquareAdv_Vel;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = 0.5*fluid[DENS]*SQR(SquareAdv_Vel) + SquareAdv_Pres*_Gamma_m1;

} // FUNCTION : HYDRO_TestProbSol_SquareAdv



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_OutputError_SquareAdv
// Description :  Compare and output the numerical and analytical solutions in the HYDRO square wave advection test
//
// Note        :  1. Invoked by "Output_TestProbErr"
//                2. This function has the similar form as Output_DumpData_Part, except that some code lines
//                   are modified to compute and output errors for this particular problem
//
// Parameter   :  BaseOnly :  Only output the base-level data
//-------------------------------------------------------------------------------------------------------
void HYDRO_OutputError_SquareAdv( const bool BaseOnly )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// set the default parameters (may need to be modified for different test problems)
// ===================================================================================================
   const OptOutputPart_t Part = OUTPUT_X;
   const real            x    = NULL_REAL;
   const real            y    = 0.0;
   const real            z    = 0.0;


// output file name
   char FileName[NCOMP_FLUID][200];
   int  ID[6];

   ID[0] = DumpID/100000;
   ID[1] = DumpID%100000/10000;
   ID[2] = DumpID%10000/1000;
   ID[3] = DumpID%1000/100;
   ID[4] = DumpID%100/10;
   ID[5] = DumpID%10;

   sprintf( FileName[0], "HYDRO_SQUARE_DENS_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[1], "HYDRO_SQUARE_MOMX_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[2], "HYDRO_SQUARE_MOMY_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[3], "HYDRO_SQUARE_MOMZ_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[4], "HYDRO_SQUARE_PRES_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
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
         fprintf( File_L1, "%5s %13s %19s %19s %19s %19s %19s\n",
                  "NGrid", "Time", "Error(DENS)", "Error(MOMX)", "Error(MOMY)", "Error(MOMZ)", "Error(PRES)" );

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

} // FUNCTION : HYDRO_OutputError_SquareAdv



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
   const bool   CheckMinPres_No = false;
   const double Gamma_m1        = GAMMA - 1.0;

   fluid[ENGY] = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                  Gamma_m1, CheckMinPres_No, NULL_REAL );


// get the analytical solution
   const real x = ( ii + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];
   const real y = ( jj + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];
   const real z = ( kk + amr->scale[lv]/2 )*amr->dh[NLEVEL-1];

// ===================================================================================================
   HYDRO_TestProbSol_SquareAdv( Anal, x, y, z, Time[0] );
// ===================================================================================================

// convert total energy to pressure
   Anal[ENGY] = CPU_GetPressure( Anal[DENS], Anal[MOMX], Anal[MOMY], Anal[MOMZ], Anal[ENGY],
                                 Gamma_m1, CheckMinPres_No, NULL_REAL );


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

   int    tmp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;

   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SquareAdv_RhoPeak,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SquareAdv_Pres,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SquareAdv_Vel,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SquareAdv_Width,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",  &tmp_int,                string );
   SquareAdv_Smooth = (bool)tmp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SquareAdv_Smooth_Eps,  string );

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default values
   if ( SquareAdv_Width <= 0.0 )
   {
      SquareAdv_Width = 0.5*amr->BoxSize[0];

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                      "SquareAdv_Width", SquareAdv_Width );
   }


// check
   if ( SquareAdv_RhoPeak    < 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "SquareAdv_RhoPeak",    SquareAdv_RhoPeak );
   if ( SquareAdv_Pres       < 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "SquareAdv_Pres",       SquareAdv_Pres );
   if ( SquareAdv_Width      < 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "SquareAdv_Width",      SquareAdv_Width );
   if ( SquareAdv_Smooth_Eps < 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "SquareAdv_Smooth_Eps", SquareAdv_Smooth_Eps );

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
