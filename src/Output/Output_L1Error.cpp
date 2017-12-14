#include "GAMER.h"

static void WriteFile( void (*AnalFunc)( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] ),
                       FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                       double L1_Err[], const OptOutputPart_t Part );




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_L1Error
// Description :  Compare the numerical and analytical solutions and output the L1 errors
//
// Note        :  1. Mainly invoked by various test problems
//                2. Similar to Output_DumpData_Part()
//                3. L1 errors are recorded in "Record__L1Err"
//
// Parameter   :  AnalFunc : Function pointer to return the analytical solution
//                           --> Usually set to the same function pointer for initializing grids
//                               (e.g., SetGridIC() in various test problems)
//                Prefix   : Prefix of the output filename
//                Part     : OUTPUT_X    : x line
//                           OUTPUT_Y    : y line
//                           OUTPUT_Z    : z line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//                x/y/z    : spatial coordinates for Part
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_L1Error( void (*AnalFunc)( real fluid[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] ),
                     const char *Prefix, const OptOutputPart_t Part, const double x, const double y, const double z )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check
   if ( Part == OUTPUT_DIAG  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for \"OUTPUT_DIAG\" !!\n" );

   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// output filename
   char FileName[NCOMP_FLUID][200];

#  if   ( MODEL == HYDRO )
   sprintf( FileName[0], "%s_Dens_%06d", Prefix, DumpID );
   sprintf( FileName[1], "%s_MomX_%06d", Prefix, DumpID );
   sprintf( FileName[2], "%s_MomY_%06d", Prefix, DumpID );
   sprintf( FileName[3], "%s_MomZ_%06d", Prefix, DumpID );
   sprintf( FileName[4], "%s_Pres_%06d", Prefix, DumpID );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   sprintf( FileName[0], "%s_Dens_%06d", Prefix, DumpID );
   sprintf( FileName[1], "%s_Real_%06d", Prefix, DumpID );
   sprintf( FileName[2], "%s_Imag_%06d", Prefix, DumpID );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL


// check if the output files already exist
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<NCOMP_FLUID; v++)
      {
         FILE *File_Check = fopen( FileName[v], "r" );

         if ( File_Check != NULL )
         {
            Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n",
                         FileName[v] );
            fclose( File_Check );

            FILE *Temp = fopen( FileName[v], "w" );
            fclose( Temp );
         }
      }
   }


// prepare to output errors
   double  dh, xx, yy, zz;
   int    *Corner  = NULL;
   double *EdgeL   = NULL;
   double *EdgeR   = NULL;
   bool    Check_x = false;
   bool    Check_y = false;
   bool    Check_z = false;

   double L1_Err[NCOMP_FLUID];
   static bool FirstTime = true;

   for (int v=0; v<NCOMP_FLUID; v++)   L1_Err[v] = 0.0;

   switch ( Part )
   {
      case OUTPUT_X    :                    Check_y = true;   Check_z = true;   break;
      case OUTPUT_Y    :  Check_x = true;                     Check_z = true;   break;
      case OUTPUT_Z    :  Check_x = true;   Check_y = true;                     break;
      case OUTPUT_DIAG :  Check_x = false;  Check_y = false;  Check_z = false;  break;
      default          :  Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [4/5/6/7] !!\n", Part );
   }


// output one MPI rank at a time
   for (int TRank=0; TRank<MPI_NRank; TRank++)
   {
      if ( MPI_Rank == TRank )
      {
         FILE *File[NCOMP_FLUID];
         for (int v=0; v<NCOMP_FLUID; v++)   File[v] = fopen( FileName[v], "a" );

//       output header
         if ( TRank == 0 )
         {
            for (int v=0; v<NCOMP_FLUID; v++)
               fprintf( File[v], "#%20s %20s %20s %20s\n", "Coord.", "Numerical", "Analytical", "Error" );
         }


//       output data
         for (int lv=0; lv<NLEVEL; lv++)
         {
            dh = amr->dh[lv];

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
//             only check the leaf patches
               if ( amr->patch[0][lv][PID]->son == -1 )
               {
                  Corner = amr->patch[0][lv][PID]->corner;
                  EdgeL  = amr->patch[0][lv][PID]->EdgeL;
                  EdgeR  = amr->patch[0][lv][PID]->EdgeR;

                  if ( Part == OUTPUT_DIAG ) // (+1,+1,+1) diagonal
                  {
                     if ( Corner[0] == Corner[1]  &&  Corner[0] == Corner[2] )
                     {
                        for (int k=0; k<PS1; k++)
                        {
                           WriteFile( AnalFunc, File, lv, PID, k, k, k, L1_Err, Part );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the target range
                     if (  !Check_x  ||  ( EdgeL[0] <= x && EdgeR[0] > x )  )
                     if (  !Check_y  ||  ( EdgeL[1] <= y && EdgeR[1] > y )  )
                     if (  !Check_z  ||  ( EdgeL[2] <= z && EdgeR[2] > z )  )
                     {
//                      check whether the cell is within the target range
                        for (int k=0; k<PS1; k++)  {  zz = EdgeL[2] + k*dh;
                                                      if ( Check_z && ( zz>z || zz+dh<=z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  yy = EdgeL[1] + j*dh;
                                                      if ( Check_y && ( yy>y || yy+dh<=y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  xx = EdgeL[0] + i*dh;
                                                      if ( Check_x && ( xx>x || xx+dh<=x ) )    continue;

                           WriteFile( AnalFunc, File, lv, PID, i, j, k, L1_Err, Part );

                        }}}
                     }
                  } // if ( Part == OUTPUT_DIAG ... else ... )
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int lv=0; lv<NLEVEL; lv++)

         for (int v=0; v<NCOMP_FLUID; v++)   fclose( File[v] );

      } // if ( MPI_Rank == TRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TRank=0; TRank<MPI_NRank; TRank++)


// gather the L1 error from all ranks and output the results
   double L1_Err_Sum[NCOMP_FLUID], Norm;
   MPI_Reduce( L1_Err, L1_Err_Sum, NCOMP_FLUID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      switch ( Part )
      {
         case OUTPUT_X    :  Norm = amr->BoxSize[0];  break;
         case OUTPUT_Y    :  Norm = amr->BoxSize[1];  break;
         case OUTPUT_Z    :  Norm = amr->BoxSize[2];  break;
         case OUTPUT_DIAG :  Norm = amr->BoxSize[0];  break;
         default          :  Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [4/5/6/7] !!\n", Part );
      }

      for (int v=0; v<NCOMP_FLUID; v++)   L1_Err_Sum[v] /= Norm;

      FILE *File_L1 = fopen( "Record__L1Err", "a" );

//    output header
      if ( FirstTime )
      {
#        if   ( MODEL == HYDRO )
         fprintf( File_L1, "#%5s %13s %19s %19s %19s %19s %19s\n",
                  "NGrid", "Time", "Error(Dens)", "Error(MomX)", "Error(MomY)", "Error(MomZ)", "Error(Pres)" );

#        elif ( MODEL == MHD )
#        warning : WAIT MHD !!!

#        elif ( MODEL == ELBDM )
         fprintf( File_L1, "#%5s %13s %19s %19s %19s\n",
                  "NGrid", "Time", "Error(Dens)", "Error(Real)", "Error(Imag)" );

#        else
#        error : unsupported MODEL !!
#        endif // MODEL

         FirstTime = false;
      } // if ( FirstTime )

//    output data
      fprintf( File_L1, "%6d %13.7e", NX0_TOT[0], Time[0] );

      for (int v=0; v<NCOMP_FLUID; v++)
      fprintf( File_L1, " %19.12e", L1_Err_Sum[v] );

      fprintf( File_L1, "\n" );

      fclose( File_L1 );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_L1Error



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  WRite the data of a single cell
//
// Note        :  1. Invoked by Output_L1Error()
//
// Parameter   :  AnalFunc : Function pointer to return the analytical solution
//                File     : File pointer
//                lv       : Target refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                L1_Err   : Array to record the L1 errors of all variables
//                Part     : OUTPUT_X    : x line
//                           OUTPUT_Y    : y line
//                           OUTPUT_Z    : z line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//
// Return      :  L1_Err
//-------------------------------------------------------------------------------------------------------
void WriteFile( void (*AnalFunc)( real fluid[], const double x, const double y, const double z, const double Time,
                                  const int lv, double AuxArray[] ),
                FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                double L1_Err[], const OptOutputPart_t Part )
{

   real fluid[NCOMP_FLUID], Anal[NCOMP_FLUID], Err[NCOMP_FLUID];


// get the numerical solution
   for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];


// convert total energy to pressure
#  if ( MODEL == HYDRO )
   const bool   CheckMinPres_No = false;
   const double Gamma_m1        = GAMMA - 1.0;

   fluid[ENGY] = CPU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                  Gamma_m1, CheckMinPres_No, NULL_REAL );
#  endif


// get the analytical solution
   const double dh = amr->dh[lv];
   const double x  = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;
   const double y  = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
   const double z  = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;

   AnalFunc( Anal, x, y, z, Time[0], lv, NULL );


// convert total energy to pressure
#  if ( MODEL == HYDRO )
   Anal[ENGY] = CPU_GetPressure( Anal[DENS], Anal[MOMX], Anal[MOMY], Anal[MOMZ], Anal[ENGY],
                                 Gamma_m1, CheckMinPres_No, NULL_REAL );
#  endif


// record the physical coordinate
   double r;

   switch ( Part )
   {
      case OUTPUT_X    : r = x;              break;
      case OUTPUT_Y    : r = y;              break;
      case OUTPUT_Z    : r = z;              break;
      case OUTPUT_DIAG : r = sqrt(3.0)*x;    break;
      default          : Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [4/5/6/7] !!\n", Part );
   }


// estimate and output errors
   for (int v=0; v<NCOMP_FLUID; v++)
   {
      Err   [v]  = FABS( Anal[v] - fluid[v] );
      L1_Err[v] += Err[v]*dh;

      fprintf( File[v], " %20.13e %20.13e %20.13e %20.13e\n", r, fluid[v], Anal[v], Err[v] );
   }

} // FUNCTION : WriteFile
