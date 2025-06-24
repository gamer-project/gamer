#include "GAMER.h"

static void WriteFile( void (*AnalFunc_Flu)( real fluid[], const double x, const double y, const double z, const double Time,
                                             const int lv, double AuxArray[] ),
                       void (*AnalFunc_Mag)( real magnetic[], const double x, const double y, const double z, const double Time,
                                             const int lv, double AuxArray[] ),
                       FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                       double L1_Err[], const OptOutputPart_t Part );


// NEXTRA: number of extra fields for computing errors
#if ( MODEL == HYDRO )
#define NEXTRA       1     // temperature
#else
#define NEXTRA       0     // none
#endif

#define NBASIC    ( NCOMP_TOTAL + NCOMP_MAG )
#define NERR      ( NBASIC + NEXTRA )




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_L1Error
// Description :  Compare the numerical and analytical solutions and output the L1 errors
//
// Note        :  1. Mainly invoked by various test problems
//                2. Similar to Output_DumpData_Part()
//                3. L1 errors are recorded in "Record__L1Err"
//                4. For MHD, this function uses the average **cell-centered** magnetic field to compute errors
//
// Parameter   :  AnalFunc_Flu : Function pointer to return the analytical solution of the fluid variables
//                               --> Usually set to the same function pointer for initializing grids
//                                   (e.g., SetGridIC() in various test problems)
//                               --> For MHD, the total energy set by this function must NOT include magnetic energy
//                AnalFunc_Mag : Function pointer to return the analytical solution of the magnetic field (MHD only)
//                               --> Usually set to the same function pointer for initializing B field
//                                   (e.g., SetBFieldIC() in various test problems)
//                Prefix       : Prefix of the output filename
//                Part         : OUTPUT_XY   : x-y plane
//                               OUTPUT_XZ   : x-z plane
//                               OUTPUT_YZ   : y-z plane
//                               OUTPUT_X    : x line
//                               OUTPUT_Y    : y line
//                               OUTPUT_Z    : z line
//                               OUTPUT_DIAG : diagonal along (+1,+1,+1)
//                               OUTPUT_BOX  : entire box
//                x/y/z        : spatial coordinates for Part
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_L1Error( void (*AnalFunc_Flu)( real fluid[], const double x, const double y, const double z, const double Time,
                                           const int lv, double AuxArray[] ),
                     void (*AnalFunc_Mag)( real magnetic[], const double x, const double y, const double z, const double Time,
                                           const int lv, double AuxArray[] ),
                     const char *Prefix, const OptOutputPart_t Part, const double x, const double y, const double z )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check
   if ( Part == OUTPUT_DIAG  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for \"OUTPUT_DIAG\" !!\n" );

   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );

   if ( AnalFunc_Flu == NULL )   Aux_Error( ERROR_INFO, "AnalyFunc_Flu == NULL !!\n" );

#  ifdef MHD
   if ( AnalFunc_Mag == NULL )   Aux_Error( ERROR_INFO, "AnalyFunc_Mag == NULL !!\n" );
#  endif


// output filename
   char FileName[NERR][2*MAX_STRING];

#  if   ( MODEL == HYDRO )
   sprintf( FileName[            0], "%s/%s_Dens_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            1], "%s/%s_MomX_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            2], "%s/%s_MomY_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            3], "%s/%s_MomZ_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            4], "%s/%s_Pres_%06d", OUTPUT_DIR, Prefix, DumpID );

   for (int v=0; v<NCOMP_PASSIVE; v++)
   sprintf( FileName[NCOMP_FLUID+v], "%s/%s_Passive%02d_%06d", OUTPUT_DIR, Prefix, v, DumpID );

#  ifdef MHD
   sprintf( FileName[NCOMP_TOTAL+0], "%s/%s_MagX_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[NCOMP_TOTAL+1], "%s/%s_MagY_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[NCOMP_TOTAL+2], "%s/%s_MagZ_%06d", OUTPUT_DIR, Prefix, DumpID );
#  endif

   sprintf( FileName[     NBASIC+0], "%s/%s_Temp_%06d", OUTPUT_DIR, Prefix, DumpID );

#  elif ( MODEL == ELBDM )
#  if   ( ELBDM_SCHEME == ELBDM_WAVE )
   sprintf( FileName[            0], "%s/%s_Dens_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            1], "%s/%s_Real_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            2], "%s/%s_Imag_%06d", OUTPUT_DIR, Prefix, DumpID );
#  elif ( ELBDM_SCHEME == ELBDM_HYBRID )
   sprintf( FileName[            0], "%s/%s_Dens_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            1], "%s/%s_Phas_%06d", OUTPUT_DIR, Prefix, DumpID );
   sprintf( FileName[            2], "%s/%s_Stub_%06d", OUTPUT_DIR, Prefix, DumpID );
#  else
#  error : ERROR : unsupported ELBDM_SCHEME !!
#  endif // ELBDM_SCHEME

   for (int v=0; v<NCOMP_PASSIVE; v++)
   sprintf( FileName[NCOMP_FLUID+v], "%s/%s_Passive%02d_%06d", OUTPUT_DIR, Prefix, v, DumpID );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


// check if the output files already exist
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<NERR; v++)
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
   int     dim;
   char    coord1, coord2, coord3;
   double  dh, xx, yy, zz;
   int    *Corner  = NULL;
   double *EdgeL   = NULL;
   double *EdgeR   = NULL;
   bool    Check_x = false;
   bool    Check_y = false;
   bool    Check_z = false;
   long    NGrid   = 0L;

   double L1_Err[NERR];
   static bool FirstTime = true;

   for (int v=0; v<NERR; v++)    L1_Err[v] = 0.0;

   switch ( Part )
   {
      case OUTPUT_XY   :                                      Check_z = true;   dim = 2;  coord1 = 'X';  coord2 = 'Y';                 break;
      case OUTPUT_YZ   :  Check_x = true;                                       dim = 2;  coord1 = 'Y';  coord2 = 'Z';                 break;
      case OUTPUT_XZ   :                    Check_y = true;                     dim = 2;  coord1 = 'X';  coord2 = 'Z';                 break;
      case OUTPUT_X    :                    Check_y = true;   Check_z = true;   dim = 1;  coord1 = 'X';                                break;
      case OUTPUT_Y    :  Check_x = true;                     Check_z = true;   dim = 1;  coord1 = 'Y';                                break;
      case OUTPUT_Z    :  Check_x = true;   Check_y = true;                     dim = 1;  coord1 = 'Z';                                break;
      case OUTPUT_DIAG :  Check_x = false;  Check_y = false;  Check_z = false;  dim = 1;  coord1 = 'X';                                break;
      case OUTPUT_BOX  :  Check_x = false;  Check_y = false;  Check_z = false;  dim = 3;  coord1 = 'X';  coord2 = 'Y';  coord3 = 'Z';  break;
      default          :  Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [1/2/3/4/5/6/7/8] !!\n", Part );
   }

// output one MPI rank at a time
   for (int TRank=0; TRank<MPI_NRank; TRank++)
   {
      if ( MPI_Rank == TRank )
      {
         FILE *File[NERR];
         for (int v=0; v<NERR; v++)    File[v] = fopen( FileName[v], "a" );

//       output header
         if ( TRank == 0 )
         {
            for (int v=0; v<NERR; v++)
            {
               fprintf( File[v], "#%*s %c", StrLen_Flt, "Coord.", coord1 );
               if ( dim >= 2 )
               fprintf( File[v], " %*s %c", StrLen_Flt, "Coord.", coord2 );
               if ( dim == 3 )
               fprintf( File[v], " %*s %c", StrLen_Flt, "Coord.", coord3 );

               fprintf( File[v], " %*s %*s %*s\n",
                        StrLen_Flt, "Numerical", StrLen_Flt, "Analytical", StrLen_Flt, "Error" );
            }
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
                           WriteFile( AnalFunc_Flu, AnalFunc_Mag, File, lv, PID, k, k, k, L1_Err, Part );
                           NGrid++;
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

                           WriteFile( AnalFunc_Flu, AnalFunc_Mag, File, lv, PID, i, j, k, L1_Err, Part );
                           NGrid++;

                        }}}
                     }
                  } // if ( Part == OUTPUT_DIAG ) ... else ...
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int lv=0; lv<NLEVEL; lv++)

         for (int v=0; v<NERR; v++)    fclose( File[v] );

      } // if ( MPI_Rank == TRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TRank=0; TRank<MPI_NRank; TRank++)


// gather the L1 error from all ranks and output the results
   double L1_Err_Sum[NERR], Norm;
   long   NGrid_Sum;
   MPI_Reduce(  L1_Err,  L1_Err_Sum, NERR, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &NGrid,  &NGrid_Sum,  1,    MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      switch ( Part )
      {
         case OUTPUT_XY   :  Norm = amr->BoxSize[0]*amr->BoxSize[1];                  break;
         case OUTPUT_XZ   :  Norm = amr->BoxSize[0]*amr->BoxSize[2];                  break;
         case OUTPUT_YZ   :  Norm =                 amr->BoxSize[1]*amr->BoxSize[2];  break;
         case OUTPUT_X    :  Norm = amr->BoxSize[0];                                  break;
         case OUTPUT_Y    :  Norm =                 amr->BoxSize[1];                  break;
         case OUTPUT_Z    :  Norm =                                 amr->BoxSize[2];  break;
         case OUTPUT_DIAG :  Norm = amr->BoxSize[0];                                  break;
         case OUTPUT_BOX  :  Norm = amr->BoxSize[0]*amr->BoxSize[1]*amr->BoxSize[2];  break;
         default          :  Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [1/2/3/4/5/6/7/8] !!\n", Part );
      }

      for (int v=0; v<NERR; v++)    L1_Err_Sum[v] /= Norm;

      char FileName_L1[2*MAX_STRING];
      sprintf( FileName_L1, "%s/Record__L1Err", OUTPUT_DIR );
      FILE *File_L1 = fopen( FileName_L1, "a" );

//    output header
      if ( FirstTime )
      {
#        if   ( MODEL == HYDRO )
         fprintf( File_L1, "#%11s %13s %*s %*s %*s %*s %*s", "NGrid", "Time", StrLen_Flt, "Error(Dens)",
                  StrLen_Flt, "Error(MomX)", StrLen_Flt, "Error(MomY)", StrLen_Flt, "Error(MomZ)", StrLen_Flt, "Error(Pres)" );

         for (int v=0; v<NCOMP_PASSIVE; v++) {
            char tmp_str[MAX_STRING];
            sprintf( tmp_str, "Error(Passive%02d)", v );
            fprintf( File_L1, " %*s", StrLen_Flt, tmp_str );
         }

#        ifdef MHD
         fprintf( File_L1, " %*s %*s %*s",
                  StrLen_Flt, "Error(MagX)", StrLen_Flt, "Error(MagY)", StrLen_Flt, "Error(MagZ)" );
#        endif

         fprintf( File_L1, " %*s", StrLen_Flt, "Error(Temp)" );

         fprintf( File_L1, "\n" );

#        elif ( MODEL == ELBDM )
#        if   ( ELBDM_SCHEME == ELBDM_WAVE )

         fprintf( File_L1, "#%11s %13s %*s %*s %*s", "NGrid", "Time", StrLen_Flt, "Error(Dens)",
                  StrLen_Flt, "Error(Real)", StrLen_Flt, "Error(Imag)" );
#        elif ( ELBDM_SCHEME == ELBDM_HYBRID )
         fprintf( File_L1, "#%11s %13s %*s %*s %*s", "NGrid", "Time", StrLen_Flt, "Error(Dens)",
                  StrLen_Flt, "Error(Phas)", StrLen_Flt, "Stub" );
#        else
#        error : ERROR : unsupported ELBDM_SCHEME !!
#        endif // ELBDM_SCHEME

         for (int v=0; v<NCOMP_PASSIVE; v++) {
            char tmp_str[MAX_STRING];
            sprintf( tmp_str, "Error(Passive%02d)", v );
            fprintf( File_L1, " %*s", StrLen_Flt, tmp_str );
         }

         fprintf( File_L1, "\n" );

#        else
#        error : ERROR : unsupported MODEL !!
#        endif // MODEL

         FirstTime = false;
      } // if ( FirstTime )

//    output data
      fprintf( File_L1, "%12ld %13.7e", NGrid_Sum, Time[0] );

      for (int v=0; v<NERR; v++)   fprintf( File_L1, BlankPlusFormat_Flt, L1_Err_Sum[v] );

      fprintf( File_L1, "\n" );

      fclose( File_L1 );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_L1Error



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  Write the data of a single cell
//
// Note        :  1. Invoked by Output_L1Error()
//
// Parameter   :  AnalFunc_Flu : Function pointer to return the analytical solution of the fluid variables
//                               --> For MHD, the total energy set by this function must NOT include magnetic energy
//                AnalFunc_Mag : Function pointer to return the analytical solution of the magnetic field (MHD only)
//                File         : File pointer
//                lv           : Target refinement level
//                PID          : Patch ID
//                i/j/k        : Cell indices within the patch
//                L1_Err       : Array to record the L1 errors of all variables
//                Part         : OUTPUT_XY   : x-y plane
//                               OUTPUT_XZ   : x-z plane
//                               OUTPUT_YZ   : y-z plane
//                               OUTPUT_X    : x line
//                               OUTPUT_Y    : y line
//                               OUTPUT_Z    : z line
//                               OUTPUT_DIAG : diagonal along (+1,+1,+1)
//                               OUTPUT_BOX  : entire box
//
// Return      :  L1_Err
//-------------------------------------------------------------------------------------------------------
void WriteFile( void (*AnalFunc_Flu)( real fluid[], const double x, const double y, const double z, const double Time,
                                      const int lv, double AuxArray[] ),
                void (*AnalFunc_Mag)( real magnetic[], const double x, const double y, const double z, const double Time,
                                      const int lv, double AuxArray[] ),
                FILE *File[], const int lv, const int PID, const int i, const int j, const int k,
                double L1_Err[], const OptOutputPart_t Part )
{

   real Nume[NERR], Anal[NERR], Err[NERR];

// get the numerical solution
   for (int v=0; v<NCOMP_TOTAL; v++)
      Nume[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];


// note that we use the cell-centered B field to compute errors
#  ifdef MHD
   MHD_GetCellCenteredBFieldInPatch( Nume+NCOMP_TOTAL, lv, PID, i, j, k, amr->MagSg[lv] );
#  endif

// get pressure and temperature
#  if ( MODEL == HYDRO )
   const bool  CheckMinPres_No = false;
   const bool  CheckMinTemp_No = false;
#  ifdef MHD
   const real *B_Nume          = Nume + NCOMP_TOTAL;
   const real  Emag_Nume       = (real)0.5*( SQR(B_Nume[MAGX]) + SQR(B_Nume[MAGY]) + SQR(B_Nume[MAGZ]) );
#  else
   const real  Emag_Nume       = NULL_REAL;
#  endif

   const real  Pres_Nume       = Hydro_Con2Pres( Nume[DENS], Nume[MOMX], Nume[MOMY], Nume[MOMZ],
                                                 Nume[ENGY], Nume+NCOMP_FLUID,
                                                 CheckMinPres_No, NULL_REAL, Emag_Nume,
                                                 EoS_DensEint2Pres_CPUPtr,
                                                 EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
   const real  Temp_Nume       = Hydro_Con2Temp( Nume[DENS], Nume[MOMX], Nume[MOMY], Nume[MOMZ],
                                                 Nume[ENGY], Nume+NCOMP_FLUID,
                                                 CheckMinTemp_No, NULL_REAL, Emag_Nume,
                                                 EoS_DensEint2Temp_CPUPtr,
                                                 EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Nume[ENGY    ] = Pres_Nume;
   Nume[NBASIC+0] = Temp_Nume;
#  endif // #if ( MODEL == HYDRO )


// get the analytical solution
   const double dh = amr->dh[lv];
   const double x  = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh;
   const double y  = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh;
   const double z  = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh;

   AnalFunc_Flu( Anal,             x, y, z, Time[0], lv, NULL );
#  ifdef MHD
   AnalFunc_Mag( Anal+NCOMP_TOTAL, x, y, z, Time[0], lv, NULL );
#  endif

// get pressure and temperature
#  if ( MODEL == HYDRO )
   const real Emag_Zero = 0.0;   // Anal[ENGY] set by AnalFunc_Flu() does NOT include magnetic energy
   const real Pres_Anal = Hydro_Con2Pres( Anal[DENS], Anal[MOMX], Anal[MOMY], Anal[MOMZ], Anal[ENGY],
                                          Anal+NCOMP_FLUID, CheckMinPres_No, NULL_REAL, Emag_Zero,
                                          EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

   const real Temp_Anal = Hydro_Con2Temp( Anal[DENS], Anal[MOMX], Anal[MOMY], Anal[MOMZ], Anal[ENGY],
                                          Anal+NCOMP_FLUID, CheckMinTemp_No, NULL_REAL, Emag_Zero,
                                          EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

   Anal[ENGY    ] = Pres_Anal;
   Anal[NBASIC+0] = Temp_Anal;
#  endif

// convert real and imaginary part to phase for wave patches in hybrid scheme
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
      Anal[PHAS] = SATAN2(Anal[IMAG], Anal[REAL]);
      Nume[PHAS] = SATAN2(Nume[IMAG], Nume[REAL]);
      Anal[STUB] = 0;
      Nume[STUB] = 0;
   }
#  endif


// record the physical coordinate(s)
   double r1, r2, r3;
   int dim;

   switch ( Part )
   {
      case OUTPUT_XY   : r1 = x;            r2 = y;           dim = 2;  break;
      case OUTPUT_XZ   : r1 = x;            r2 = z;           dim = 2;  break;
      case OUTPUT_YZ   : r1 = y;            r2 = z;           dim = 2;  break;
      case OUTPUT_X    : r1 = x;                              dim = 1;  break;
      case OUTPUT_Y    : r1 = y;                              dim = 1;  break;
      case OUTPUT_Z    : r1 = z;                              dim = 1;  break;
      case OUTPUT_DIAG : r1 = sqrt(3.0)*x;                    dim = 1;  break;
      case OUTPUT_BOX  : r1 = x;            r2 = y;  r3 = z;  dim = 3;  break;
      default          : Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [1/2/3/4/5/6/7/8] !!\n", Part );
   }

   const double dv = (dim==1) ? dh : (dim==2) ? SQR(dh) : CUBE(dh);

// estimate and output errors
   for (int v=0; v<NERR; v++)
   {
      Err   [v]  = FABS( Anal[v] - Nume[v] );
      L1_Err[v] += Err[v]*dv;

      fprintf( File[v], "  "                         );
      fprintf( File[v], BlankPlusFormat_Flt, r1      );
      if ( dim >= 2 ) {
      fprintf( File[v], "  "                         );
      fprintf( File[v], BlankPlusFormat_Flt, r2      ); }
      if ( dim == 3 ) {
      fprintf( File[v], "  "                         );
      fprintf( File[v], BlankPlusFormat_Flt, r3      ); }
      fprintf( File[v], BlankPlusFormat_Flt, Nume[v] );
      fprintf( File[v], BlankPlusFormat_Flt, Anal[v] );
      fprintf( File[v], BlankPlusFormat_Flt, Err[v]  );
      fprintf( File[v], "\n");
   }

} // FUNCTION : WriteFile
