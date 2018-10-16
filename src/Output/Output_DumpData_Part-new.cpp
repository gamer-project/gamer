#include "GAMER.h"

static void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                       const int ii, const int jj, const int kk );

#if (MODEL == SR_HYDRO)
void CPU_4Velto3Vel( const real In[], real Out[] );
void CPU_Con2Pri( const real In[], real Out[], const real Gamma );
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpData_Part
// Description :  Output part of data in the ASCII form
//
// Parameter   :  Part     : OUTPUT_XY   : xy plane
//                           OUTPUT_YZ   : yz plane
//                           OUTPUT_XZ   : xz plane
//                           OUTPUT_X    : x  line
//                           OUTPUT_Y    : y  line
//                           OUTPUT_Z    : z  line
//                           OUTPUT_DIAG : diagonal along (+1,+1,+1)
//
//                BaseOnly : Only output the base-level data
//
//                x        : x coordinate
//                y        : y coordinate
//                z        : z coordinate
//
//                FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const double x, const double y,
                           const double z, const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check the input parameters
   if ( Part != OUTPUT_XY  &&  Part != OUTPUT_YZ  &&  Part != OUTPUT_XZ  &&
        Part != OUTPUT_X   &&  Part != OUTPUT_Y   &&  Part != OUTPUT_Z   &&  Part != OUTPUT_DIAG )
      Aux_Error( ERROR_INFO, "unsupported option \"Part = %d\" [0 ~ 6] !!\n", Part );

   if (  ( Part == OUTPUT_YZ  ||  Part == OUTPUT_Y  ||  Part == OUTPUT_Z )  &&
         ( x < 0.0  ||  x >= amr->BoxSize[0] )  )
      Aux_Error( ERROR_INFO, "incorrect x (out of range [0<=X<%lf]) !!\n", amr->BoxSize[0] );

   if (  ( Part == OUTPUT_XZ  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Z )  &&
         ( y < 0.0  ||  y >= amr->BoxSize[1] )  )
      Aux_Error( ERROR_INFO, "incorrect y (out of range [0<=Y<%lf]) !!\n", amr->BoxSize[1] );

   if (  ( Part == OUTPUT_XY  ||  Part == OUTPUT_X  ||  Part == OUTPUT_Y )  &&
         ( z < 0.0  ||  z >= amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "incorrect z (out of range [0<=Z<%lf]) !!\n", amr->BoxSize[2] );

   if ( Part == OUTPUT_DIAG  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for \"OUTPUT_DIAG\" !!\n" );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// check if the file already exists
   if ( MPI_Rank == 0 )
   {
      if ( Aux_CheckFileExist(FileName) )
      {
         Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

         FILE *Temp = fopen( FileName, "w" );
         fclose( Temp );
      }
   }


   const double dh_min = amr->dh[NLEVEL-1];
   const int    NLv    = ( BaseOnly ) ? 1 : NLEVEL;

   int     ii, jj, kk, scale;
   double  dh, xx, yy, zz;    // xx,yy,zz => physical coordinates of cell left edge
   int    *Corner  = NULL;    // patch corner in scale
   double *EdgeL   = NULL;    // patch corner in physical coord.
   double *EdgeR   = NULL;
   bool    Check_x = false;
   bool    Check_y = false;
   bool    Check_z = false;

   switch ( Part )
   {
      case OUTPUT_XY :                                      Check_z = true;   break;
      case OUTPUT_YZ :  Check_x = true;                                       break;
      case OUTPUT_XZ :                    Check_y = true;                     break;
      case OUTPUT_X  :                    Check_y = true;   Check_z = true;   break;
      case OUTPUT_Y  :  Check_x = true;                     Check_z = true;   break;
      case OUTPUT_Z  :  Check_x = true;   Check_y = true;                     break;

      case OUTPUT_DIAG :
      case OUTPUT_PART_NONE : break; // do nothing
   }


   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File = fopen( FileName, "a" );

//       output header
         if ( TargetMPIRank == 0 )
         {
            fprintf( File, "#%10s %10s %10s %20s %20s %20s", "i[1]", "j[2]", "k[3]", "x[4]", "y[5]", "z[6]" );


#           if ( MODEL == SR_HYDRO )
            fprintf( File, "%14s", "Dens[7]" );
            fprintf( File, "%14s", "MomX[8]" );
            fprintf( File, "%14s", "MomY[9]" );
            fprintf( File, "%14s", "MomZ[10]" );
            fprintf( File, "%14s", "Engy[11]" );
#           else
            for (int v=0; v<NCOMP_TOTAL; v++)
            fprintf( File, "%14s", FieldLabel[v] );
#           endif

#           ifdef GRAVITY
            if ( OPT__OUTPUT_POT )
            fprintf( File, "%14s", PotLabel );
#           endif

//          other derived fields
#           if ( MODEL == HYDRO  ||  MODEL == MHD )
            fprintf( File, "%14s", "Pressure[12]" );
#           elif ( MODEL == SR_HYDRO )
            fprintf( File, "%14s", "PrimDens[12]" );
            fprintf( File, "%14s", "Ux[13]" ); // 4-velocity
            fprintf( File, "%14s", "Uy[14]" );
            fprintf( File, "%14s", "Uz[15]" );
            fprintf( File, "%14s", "Pressure[16]" );
            fprintf( File, "%19s", "Lorentz Fac[17]" );
            fprintf( File, "%14s", "Vx[18]" ); // 3-velocity
            fprintf( File, "%14s", "Vy[19]" );
            fprintf( File, "%14s", "Vz[20]" );
#           endif
            fprintf( File, "\n" );
         } // if ( TargetMPIRank == 0 )


//       output data
         for (int lv=0; lv<NLv; lv++)
         {
            dh    = amr->dh   [lv];
            scale = amr->scale[lv];

            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
//             output the patch data only if it has no son (if the option "BaseOnly" is turned off)
               if ( amr->patch[0][lv][PID]->son == -1  ||  BaseOnly )
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
                           kk = Corner[2] + k*scale;

                           WriteFile( File, lv, PID, k, k, k, kk, kk, kk );
                        }
                     }
                  } // if ( Part == OUTPUT_DIAG )


                  else // x/y/z lines || xy/yz/xz slices
                  {
//                   check whether the patch corner is within the target range
                     if (  !Check_x  ||  ( EdgeL[0]<=x && EdgeR[0]>x )  )
                     if (  !Check_y  ||  ( EdgeL[1]<=y && EdgeR[1]>y )  )
                     if (  !Check_z  ||  ( EdgeL[2]<=z && EdgeR[2]>z )  )
                     {
//                      check whether the cell is within the target range
                        for (int k=0; k<PS1; k++)  {  kk = Corner[2] + k*scale;  zz = kk*dh_min;
                                                      if ( Check_z && ( zz>z || zz+dh<=z ) )    continue;

                        for (int j=0; j<PS1; j++)  {  jj = Corner[1] + j*scale;  yy = jj*dh_min;
                                                      if ( Check_y && ( yy>y || yy+dh<=y ) )    continue;

                        for (int i=0; i<PS1; i++)  {  ii = Corner[0] + i*scale;  xx = ii*dh_min;
                                                      if ( Check_x && ( xx>x || xx+dh<=x ) )    continue;

                           WriteFile( File, lv, PID, i, j, k, ii, jj, kk );

                        }}}
                     } // if patch corner is within the target range

                  } // if ( Part == OUTPUT_DIAG ... else ... )
               } // if ( amr->patch[0][lv][PID]->son == -1 )
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // for (int lv=0; lv<NLv; lv++)

         fclose( File );

      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Output_DumpData_Part



//-------------------------------------------------------------------------------------------------------
// Function    :  WriteFile
// Description :  Output data to file
//
// Parameter   :  File     : File pointer
//                lv       : Target refinement level
//                PID      : Patch ID
//                i/j/k    : Cell indices within the patch
//                ii/jj/kk : Cell scale indices in the simulation domain
//-------------------------------------------------------------------------------------------------------
void WriteFile( FILE *File, const int lv, const int PID, const int i, const int j, const int k,
                const int ii, const int jj, const int kk )
{

   const double dh_min  = amr->dh[TOP_LEVEL];
   const double scale_2 = 0.5*amr->scale[lv];
   real u[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)   u[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

// output cell indices and coordinates
   fprintf( File, " %10d %10d %10d %20.14e %20.14e %20.14e",
            ii, jj, kk, (ii+scale_2)*dh_min, (jj+scale_2)*dh_min, (kk+scale_2)*dh_min );

// output all conserved variables in the fluid or sr-fluid array
   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13.6e", u[v] );

#  if (MODEL == SR_HYDRO)
   real Pri4Vel[NCOMP_FLUID];
   CPU_Con2Pri(u,Pri4Vel,GAMMA);

// output all primitive variables in the sr-fluid array
   for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13.6e", Pri4Vel[v] );
#  endif

// output potential
#  ifdef GRAVITY
   if ( OPT__OUTPUT_POT )
   fprintf( File, " %13.6e", amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i] );
#  endif

// output other derived fields
#  if   ( MODEL == HYDRO )
   const bool CheckMinPres_Yes = true;
   fprintf( File, " %13.6e", CPU_GetPressure(u[DENS], u[MOMX], u[MOMY], u[MOMZ], u[ENGY], GAMMA-1.0, CheckMinPres_Yes, MIN_PRES) );
#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  elif ( MODEL == SR_HYDRO )
// output Lorentz factor
   fprintf( File, " %17.6e", SQRT(1+SQR(Pri4Vel[1])+SQR(Pri4Vel[2])+SQR(Pri4Vel[3])) );

   real Pri3Vel[NCOMP_FLUID];

   CPU_4Velto3Vel(Pri4Vel,Pri3Vel);

// output 3-velocity
   fprintf( File, " %13.6e", Pri3Vel[1] );
   fprintf( File, " %13.6e", Pri3Vel[2] );
   fprintf( File, " %13.6e", Pri3Vel[3] );

#  endif // MODEL

   fprintf( File, "\n" );

} // FUNCTION : WriteFile
