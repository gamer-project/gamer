#include "GAMER.h"

#if ( MODEL == ELBDM )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void ELBDM_TestProbSol_JeansInstability_Physical( real fluid[], const double x, const double y, const double z,
                                                         const double Time );
static void ELBDM_OutputError_JeansInstability_Physical( const bool BaseOnly );
static void WriteFile( FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk, double L1_Err[], const OptOutputPart_t Part  );


// global variables in the ELBDM Jeans instability test
// =======================================================================================
static real Jeans_WaveLength;    // wave length
static real Jeans_WaveK;         // wave number
static real Jeans_WaveKj;        // critical wave number
static real Jeans_Omega;         // wave angular frequency
static real Jeans_RealAmp;       // wave real part amplitude
static real Jeans_ImagAmp;       // wave imaginary part amplitude
static real Jeans_Sign;          // stable   : (+1/-1) --> (right/left-moving wave)
                                 // unstable : (+1/-1) --> (growing/decaying mode)
static real Jeans_Phase0;        // initial phase shift
static bool Jeans_Stable;        // true/false --> Jeans stable/unstable
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the ELBDM Jeans instability test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "ELBDM_TestProbSol_JeansInstability"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "ELBDM physical-frame Jeans instability";

// check
#  if ( MODEL != ELBDM )
#  error : ERROR : "MODEL != ELBDM" in the ELBDM Jeans instability test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the ELBDM Jeans instability test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the ELBDM Jeans instability test !!
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in the %s test !!\n", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = ELBDM_TestProbSol_JeansInstability_Physical;
   Output_TestProbErr_Ptr = ELBDM_OutputError_JeansInstability_Physical;


// set global variables
   Jeans_WaveLength = amr->BoxSize[0]/sqrt(3.0);
   Jeans_WaveK      = 2.0*M_PI/Jeans_WaveLength;
   Jeans_WaveKj     = POW( 16.0*M_PI*NEWTON_G*SQR(ELBDM_ETA), 0.25 );
   Jeans_Stable     = ( Jeans_WaveK > Jeans_WaveKj ) ? true : false;
   Jeans_Omega      = ( Jeans_Stable ) ? 0.5/ELBDM_ETA*sqrt( POW(Jeans_WaveK,  4.0) - POW(Jeans_WaveKj, 4.0) )
                                       : 0.5/ELBDM_ETA*sqrt( POW(Jeans_WaveKj, 4.0) - POW(Jeans_WaveK,  4.0) );
   Jeans_Sign       = +1.0;
   Jeans_Phase0     = 0.0;
   Jeans_RealAmp    = 1.e-6;
   Jeans_ImagAmp    = 2.0*ELBDM_ETA*Jeans_Omega*Jeans_RealAmp/(Jeans_WaveK*Jeans_WaveK);


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : wave number            = %13.7e\n", Jeans_WaveK );
      Aux_Message( stdout, "       critial wave number    = %13.7e\n", Jeans_WaveKj );
      Aux_Message( stdout, "       stable                 = %s\n",     Jeans_Stable ? "YES" : "NO" );
      Aux_Message( stdout, "       wave angular frequency = %13.7e\n", Jeans_Omega );
      Aux_Message( stdout, "       real part amplitude    = %13.7e\n", Jeans_RealAmp );
      Aux_Message( stdout, "       imag part amplitude    = %13.7e\n", Jeans_ImagAmp );
      Aux_Message( stdout, "       sign (grow/decay;R/L)  = %13.7e\n", Jeans_Sign );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
// End_T : (stable/unstable) --> (1 period/grow by a factor of 50)
   const double End_T_Default    = ( Jeans_Stable) ? 2.0*M_PI/Jeans_Omega : log(50.0)/Jeans_Omega;
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
// Function    :  ELBDM_TestProbSol_JeansInstability_Physical
// Description :  Calculate the analytical solution in the ELBDM physical-frame Jeans instability test
//
// Note        :  1. Wave vector is along the diagonal direction
//                2. Background density is assumed to be ONE
//                3. This function is invoked by "ELBDM_Init_StartOver_AssignData" and "Output_TestProbErr"
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void ELBDM_TestProbSol_JeansInstability_Physical( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double r = 1.0/sqrt(3.0)*( x + y + z );
   double Phase;

   if ( Jeans_Stable )
   {
      Phase       = Jeans_WaveK*r - Jeans_Sign*Jeans_Omega*Time + Jeans_Phase0;
      fluid[REAL] = 1.0 + Jeans_RealAmp*cos( Phase );
      fluid[IMAG] =       Jeans_ImagAmp*sin( Phase )*Jeans_Sign;
   }

   else
   {
      Phase       = Jeans_WaveK*r + Jeans_Phase0;
      fluid[REAL] = 1.0 + Jeans_RealAmp*cos( Phase )*exp( Jeans_Sign*Jeans_Omega*Time );
      fluid[IMAG] =       Jeans_ImagAmp*cos( Phase )*exp( Jeans_Sign*Jeans_Omega*Time )*Jeans_Sign;
   }

   fluid[DENS] = fluid[REAL]*fluid[REAL] + fluid[IMAG]*fluid[IMAG];

} // FUNCTION : ELBDM_TestProbSol_JeansInstability_Physical



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_OutputError_JeansInstability_Physical
// Description :  Compare and output the numerical and analytical solutions in the ELBDM physical-frame
//                Jeans instability test problem
//
// Note        :  1. Invoked by "Output_TestProbErr"
//                2. This function has the similar form as Output_DumpData_Part, except that some code lines
//                   are modified to compute and output errors for this particular problem
//
// Parameter   :  BaseOnly :  Only output the base-level data
//-------------------------------------------------------------------------------------------------------
void ELBDM_OutputError_JeansInstability_Physical( const bool BaseOnly )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_Check_Synchronization( Time[0], Time[lv], __FUNCTION__, true );


// set the default parameters (may need to be modified for different test problems)
// ===================================================================================================
   const OptOutputPart_t Part = OUTPUT_DIAG;
   const real            x    = NULL_REAL;
   const real            y    = NULL_REAL;
   const real            z    = NULL_REAL;


// output file name
   char FileName[NCOMP][200];
   int  ID[6];

   ID[0] = DumpID/100000;
   ID[1] = DumpID%100000/10000;
   ID[2] = DumpID%10000/1000;
   ID[3] = DumpID%1000/100;
   ID[4] = DumpID%100/10;
   ID[5] = DumpID%10;

#  if   ( MODEL == HYDRO )
   sprintf( FileName[0], "HYDRO_XX_DENS_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[1], "HYDRO_XX_MOMX_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[2], "HYDRO_XX_MOMY_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[3], "HYDRO_XX_MOMZ_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );
   sprintf( FileName[4], "HYDRO_XX_PRES_%d%d%d%d%d%d", ID[0], ID[1], ID[2], ID[3], ID[4], ID[5] );

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
      for (int v=0; v<NCOMP; v++)
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

   double L1_Err[NCOMP];
   static bool FirstTime = true;

   for (int v=0; v<NCOMP; v++)   L1_Err[v] = 0.0;


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
         FILE *File[NCOMP];
         for (int v=0; v<NCOMP; v++)   File[v] = fopen( FileName[v], "a" );

//       output header
         if ( TargetMPIRank == 0 )
         {
            for (int v=0; v<NCOMP; v++)
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

         for (int v=0; v<NCOMP; v++)   fclose( File[v] );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


// gather the L1 error from all ranks and output the results
   double L1_Err_Sum[NCOMP];
   MPI_Reduce( L1_Err, L1_Err_Sum, NCOMP, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<NCOMP; v++)   L1_Err_Sum[v] /= amr->BoxSize[0];

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

      for (int v=0; v<NCOMP; v++)
      fprintf( File_L1, " %19.12e", L1_Err_Sum[v] );

      fprintf( File_L1, "\n" );


      fclose( File_L1 );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : ELBDM_OutputError_JeansInstability_Physical



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

   real fluid[NCOMP], Anal[NCOMP], Err[NCOMP];

// get the numerical solution
   for (int v=0; v<NCOMP; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

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
   ELBDM_TestProbSol_JeansInstability_Physical( Anal, x, y, z, Time[0] );
// ===================================================================================================


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
   for (int v=0; v<NCOMP; v++)
   {
      Err   [v]  = fabs( Anal[v] - fluid[v] );
      L1_Err[v] += Err[v]*amr->dh[lv];

      fprintf( File[v], "%9.7f %20.13e %20.13e %20.13e\n", r, fluid[v], Anal[v], Err[v] );
   }

} // FUNCTION : WriteFile



#endif // #if ( MODEL == ELBDM )
