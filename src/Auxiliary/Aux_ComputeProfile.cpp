#include "GAMER.h"


// indices of fields not defined in Macro.h
static const int INTERNAL_ENGY  = 97;
static const int VRAD           = 98;
static const int PRESSURE       = 99;


//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeProfile
// Description :  Compute the average radial profile of a target field
//
// Note        :  1. Results will be stored in the input "Prof" object
//                   --> Prof->Radius[]: Radial coordinate at each bin
//                       Prof->Data  []: Profile data at each bin
//                       Prof->Weight[]: Total weighting at each bin
//                       Prof->NCell []: Number of cells at each bin
//                       Prof->NBin    : Total number of bins
//                   --> See the "Profile_t" structure defined in "include/Profile.h" for details
//                   --> These arrays will be free'd when deleting "Prof"
//                2. Maximum radius adopted when actually computing the profile may be larger than the input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//                3. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same profile data after invoking this function
//
// Parameter   :  Prof        : Profile_t object to store the results
//                Center      : Target center coordinates
//                r_max_input : Maximum radius for computing the profile
//                              --> See also "Note-2" above
//                dr_min      : Minimum bin size
//                              --> For linear bin, this is the size of all bins
//                                  For log    bin, this is the size of the 0th bin
//                LogBin      : true/false --> log/linear bins
//                LogBinRatio : Ratio of adjacent log bins
//                              --> Right edge of log bin n = dr_min*LogBinRatio^n
//                RemoveEmpty : true  --> remove empty bins from the data
//                              false --> these empty bins will still be in the profile arrays with
//                                        Data[empty_bin]=Weight[empty_bin]=NCell[empty_bin]=0
//                Quantity    : Quantity to be averaged spherically
//                              Support field indicies defined in Macro.h, and VRAD (98) and PRESSURE (99)
//                              The weight function is the cell volume.
//
// Example     :  Profile_t Prof;
//
//                const double Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double MaxRadius      = 0.5*amr->BoxSize[0];
//                const double MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool   LogBin         = true;
//                const double LogBinRatio    = 1.25;
//                const bool   RemoveEmptyBin = true;
//                const int    Quantity       = DENS;
//
//                Aux_ComputeProfile( &Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin, Quantity );
//
//                if ( MPI_Rank == 0 )
//                {
//                   FILE *File = fopen( "Profile.txt", "w" );
//                   fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
//                   for (int b=0; b<Prof.NBin; b++)
//                      fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
//                               Prof.Radius[b], Prof.Data[b], Prof.Weight[b], Prof.NCell[b] );
//                   fclose( File );
//                }
//
// Return      :  Prof
//-------------------------------------------------------------------------------------------------------
void Aux_ComputeProfile( Profile_t *Prof, const double Center[], const double r_max_input, const double dr_min,
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const int Quantity )
{

// check
#  ifdef GAMER_DEBUG
   if ( r_max_input <= 0.0 )
      Aux_Error( ERROR_INFO, "r_max_input (%14.7e) <= 0.0 !!\n", r_max_input );

   if ( dr_min <= 0.0 )
      Aux_Error( ERROR_INFO, "dr_min (%14.7e) <= 0.0 !!\n", dr_min );

   if ( LogBin  &&  LogBinRatio <= 1.0 )
      Aux_Error( ERROR_INFO, "LogBinRatio (%14.7e) <= 1.0 !!\n", LogBinRatio );
#  endif


// get the total number of radial bins and the corresponding maximum radius
   if ( LogBin )
   {
      Prof->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
      Prof->MaxRadius = dr_min*pow( LogBinRatio, Prof->NBin-1 );
   }

   else // linear bin
   {
      Prof->NBin      = (int)ceil( r_max_input / dr_min );
      Prof->MaxRadius = dr_min*Prof->NBin;
   }


// record profile parameters
   for (int d=0; d<3; d++)    Prof->Center[d] = Center[d];

   Prof->LogBin = LogBin;

   if ( LogBin )  Prof->LogBinRatio = LogBinRatio;


// allocate all member arrays of Prof
   Prof->AllocateMemory();


// record radial coordinates
   if ( LogBin )
      for (int b=0; b<Prof->NBin; b++)    Prof->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
   else
      for (int b=0; b<Prof->NBin; b++)    Prof->Radius[b] = (b+0.5)*dr_min;


// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double **OMP_Data=NULL, **OMP_Weight=NULL;
   long   **OMP_NCell=NULL;

   Aux_AllocateArray2D( OMP_Data,   NT, Prof->NBin );
   Aux_AllocateArray2D( OMP_Weight, NT, Prof->NBin );
   Aux_AllocateArray2D( OMP_NCell,  NT, Prof->NBin );


// collect profile dat in this rank
   const double r_max2 = SQR( Prof->MaxRadius );

#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int b=0; b<Prof->NBin; b++)
      {
         OMP_Data  [TID][b] = 0.0;
         OMP_Weight[TID][b] = 0.0;
         OMP_NCell [TID][b] = 0;
      }

      for (int lv=0; lv<NLEVEL; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );

#        pragma omp for schedule( runtime )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            if ( amr->patch[0][lv][PID]->son != -1 )  continue;

            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

            for (int k=0; k<PS1; k++)  {  const double dz = z0 + k*dh;
            for (int j=0; j<PS1; j++)  {  const double dy = y0 + j*dh;
            for (int i=0; i<PS1; i++)  {  const double dx = x0 + i*dh;

               const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

               if ( r2 < r_max2 )
               {
                  const double r   = sqrt( r2 );
                  const int    bin = ( LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(LogBinRatio) ) + 1  )
                                                : int( r/dr_min );
//                prevent from round-off errors
                  if ( bin >= Prof->NBin )   continue;

//                check
#                 ifdef GAMER_DEBUG
                  if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                 endif

//                user-specified quantity; for case of MODEL = HYDRO
                  switch ( Quantity )
                  {
                     case DENS:
                     case ENGY:
                     case MOMX:
                     case MOMY:
                     case MOMZ:
                     {
                        OMP_Data  [TID][bin] += amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[Quantity][k][j][i]*dv;
                        OMP_Weight[TID][bin] += dv;
                     }
                     break;

                     case VRAD:
                     {
                        const double Phi      = ATAN2(dy, dx);
                        const double cosPhi   = COS(Phi);
                        const double sinPhi   = SIN(Phi);
                        const double cosTheta = dz / r;
                        const double sinTheta = SQRT(1. - SQR(cosTheta));

                        const double MomRad = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i]*sinTheta*cosPhi
                                            + amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i]*sinTheta*sinPhi
                                            + amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i]*cosTheta;

                        OMP_Data  [TID][bin] += MomRad*dv;
                        OMP_Weight[TID][bin] += amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i]*dv;
                     }
                     break;

                     case PRESSURE:
                     {
                        const double Pres = Hydro_GetPressure(amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i],
                                                              amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i],
                                                              amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i],
                                                              amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i],
                                                              amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i],
                                                              GAMMA - (real)1.0, false, NULL_REAL);
                        OMP_Data  [TID][bin] += Pres*dv;
                        OMP_Weight[TID][bin] += dv;
                     }
                     break;

                     case INTERNAL_ENGY:
                     {
                        const double intengy = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY][k][j][i]
                                             - 0.5 * ( SQR( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX][k][j][i] )
                                                     + SQR( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY][k][j][i] )
                                                     + SQR( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ][k][j][i] ) )
                                             / amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

                        OMP_Data  [TID][bin] += intengy*dv;
                        OMP_Weight[TID][bin] += dv;
                     }
                     break;

                     default:
                        Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Quantity", Quantity );
                  } // switch ( Quantity )

                  OMP_NCell [TID][bin] ++;
               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=0; lv<NLEVEL; lv++)
   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int b=0; b<Prof->NBin; b++)
   {
      Prof->Data  [b]  = OMP_Data  [0][b];
      Prof->Weight[b]  = OMP_Weight[0][b];
      Prof->NCell [b]  = OMP_NCell [0][b];
   }

   for (int t=1; t<NT; t++)
   for (int b=0; b<Prof->NBin; b++)
   {
      Prof->Data  [b] += OMP_Data  [t][b];
      Prof->Weight[b] += OMP_Weight[t][b];
      Prof->NCell [b] += OMP_NCell [t][b];
   }


// free per-thread arrays
   Aux_DeallocateArray2D( OMP_Data );
   Aux_DeallocateArray2D( OMP_Weight );
   Aux_DeallocateArray2D( OMP_NCell );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   if ( MPI_Rank == 0 )
   {
      MPI_Reduce( MPI_IN_PLACE, Prof->Data,   Prof->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( MPI_IN_PLACE, Prof->Weight, Prof->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( MPI_IN_PLACE, Prof->NCell , Prof->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
   }

   else
   {
      MPI_Reduce( Prof->Data,   NULL,         Prof->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( Prof->Weight, NULL,         Prof->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( Prof->NCell,  NULL,         Prof->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
   }
#  endif


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   for (int b=0; b<Prof->NBin; b++)
   {
//    skip empty bins since both their data and weight are zero
      if ( Prof->NCell[b] > 0L )
         switch ( Quantity )
         {
            case DENS    :
            case ENGY    :
            case MOMX    :
            case MOMY    :
            case MOMZ    :
            case PRESSURE:
               Prof->Data[b] /= Prof->Weight[b];
            break;

            case VRAD:
//             Avoid division by zero when denisty is zero
               if ( Prof->Weight[b] > 0.0 )   Prof->Data[b] /= Prof->Weight[b];
            break;
         } // switch ( Quantity )
   }


// broadcast data to all ranks
   MPI_Bcast( Prof->Data,   Prof->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( Prof->Weight, Prof->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   MPI_Bcast( Prof->NCell,  Prof->NBin, MPI_LONG,   0, MPI_COMM_WORLD );


// remove the empty bins
// --> all ranks do the same work so that no data broadcast is required
   if ( RemoveEmpty )
   for (int b=0; b<Prof->NBin; b++)
   {
      if ( Prof->NCell[b] != 0L )   continue;

      for (int b_up=b+1; b_up<Prof->NBin; b_up++)
      {
         const int b_up_m1 = b_up - 1;

         Prof->Radius[b_up_m1] = Prof->Radius[b_up];
         Prof->Data  [b_up_m1] = Prof->Data  [b_up];
         Prof->Weight[b_up_m1] = Prof->Weight[b_up];
         Prof->NCell [b_up_m1] = Prof->NCell [b_up];
      }

//    reset the total number of bins
      Prof->NBin --;

//    reset the maximum radius if we are removing the last bin
      if ( b == Prof->NBin )
      {
         if ( LogBin )
            Prof->MaxRadius = dr_min*pow( LogBinRatio, Prof->NBin-1 );
         else
            Prof->MaxRadius = dr_min*Prof->NBin;
      }

//    reduce counter since all bins above b have been shifted downward
      b --;
   } // for (int b=0; b<Prof->NBin; b++)

} // FUNCTION : Aux_ComputeProfile
