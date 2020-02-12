#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeProfile
// Description :  Compute the average radial profile of target field(s)
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
//                4. Use cell volume as the weighting of each cell
//                   --> Will support other weighting functions in the future
//
// Parameter   :  Prof        : Profile_t object array to store the results
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
//                TVar        : Target variables for computing the profiles
//                              --> Supported field indicies (defined in Macro.h):
//                                     HYDRO : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELR, _PRES, _EINT [, _POTE]
//                                     ELBDM : _DENS, _REAL, _IMAG [, _POTE]
//                              --> For a passive scalar with an integer field index FieldIdx returned by AddField(),
//                                  one can convert it to a bitwise field index by BIDX(FieldIdx)
//                NProf       : Number of Profile_t objects in Prof
//                SingleLv    : Only consider patches on the specified level
//                              --> If SingleLv<0, loop over all levels
//
// Example     :  const double Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double MaxRadius      = 0.5*amr->BoxSize[0];
//                const double MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool   LogBin         = true;
//                const double LogBinRatio    = 1.25;
//                const bool   RemoveEmptyBin = true;
//                const long   TVar           = _DENS | _PRES;
//                const int    NProf          = 2;
//                const int    SingleLv       = -1;
//
//                Profile_t Prof_Dens, Prof_Pres;
//                Profile_t *Prof[] = { &Prof_Dens, &Prof_Pres };
//
//                Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
//                                    TVar, NProf, SingleLv );
//
//                if ( MPI_Rank == 0 )
//                {
//                   for (int p=0; p<NProf; p++)
//                   {
//                      char Filename[MAX_STRING];
//                      sprintf( Filename, "Profile%d.txt", p );
//                      FILE *File = fopen( Filename, "w" );
//                      fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
//                      for (int b=0; b<Prof[p]->NBin; b++)
//                         fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
//                                  Prof[p]->Radius[b], Prof[p]->Data[b], Prof[p]->Weight[b], Prof[p]->NCell[b] );
//                      fclose( File );
//                   }
//                }
//
// Return      :  Prof
//-------------------------------------------------------------------------------------------------------
void Aux_ComputeProfile( Profile_t *Prof[], const double Center[], const double r_max_input, const double dr_min,
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const long TVar,
                         const int NProf, const int SingleLv )
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


// extract the target field indices
// --> treat intrinsic and derived fields separately
   int NFlu, NFound, TVarIdxList_Flu[NCOMP_TOTAL];

// intrinsic fields
   NFlu = 0;
   for (int v=0; v<NCOMP_TOTAL; v++)
      if ( TVar & (1L<<v) )   TVarIdxList_Flu[ NFlu ++ ] = v;

   NFound = NFlu;

// gravitational potential
#  ifdef GRAVITY
   bool PrepPote=false;
   if ( TVar & _POTE )  {  PrepPote=true;  NFound++; };
#  endif

// derived fields
#  if ( MODEL == HYDRO )
   bool PrepVelR=false, PrepPres=false, PrepEint=false;
   if ( TVar & _VELR )  {  PrepVelR=true;  NFound++; };
   if ( TVar & _PRES )  {  PrepPres=true;  NFound++; };
   if ( TVar & _EINT )  {  PrepEint=true;  NFound++; };

#  else
#  error : ERROR : unsupported MODEL !!
#  endif

   if ( NFound != NProf )
      Aux_Error( ERROR_INFO, "inconsistent number of fields (found %d != expected %d) !!\n", NFound, NProf );


// initialize the profile objects
   for (int p=0; p<NProf; p++)
   {
//    get the total number of radial bins and the corresponding maximum radius
      if ( LogBin )
      {
         Prof[p]->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
         Prof[p]->MaxRadius = dr_min*pow( LogBinRatio, Prof[p]->NBin-1 );
      }

      else // linear bin
      {
         Prof[p]->NBin      = (int)ceil( r_max_input / dr_min );
         Prof[p]->MaxRadius = dr_min*Prof[p]->NBin;
      }


//    record profile parameters
      for (int d=0; d<3; d++)    Prof[p]->Center[d] = Center[d];

      Prof[p]->LogBin = LogBin;

      if ( LogBin )  Prof[p]->LogBinRatio = LogBinRatio;


//    allocate all member arrays of Prof
      Prof[p]->AllocateMemory();


//    record radial coordinates
      if ( LogBin )
         for (int b=0; b<Prof[0]->NBin; b++)    Prof[p]->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
      else
         for (int b=0; b<Prof[0]->NBin; b++)    Prof[p]->Radius[b] = (b+0.5)*dr_min;

   } // for (int p=0; p<NProf; p++)


// allocate memory for the per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double ***OMP_Data=NULL, ***OMP_Weight=NULL;
   long   ***OMP_NCell=NULL;

   Aux_AllocateArray3D( OMP_Data,   NProf, NT, Prof[0]->NBin );
   Aux_AllocateArray3D( OMP_Weight, NProf, NT, Prof[0]->NBin );
   Aux_AllocateArray3D( OMP_NCell,  NProf, NT, Prof[0]->NBin );


// collect profile data in this rank
#  if ( MODEL == HYDRO )
   const real   Gamma_m1    = GAMMA - (real)1.0;
#  endif
   const double r_max2      = SQR( Prof[0]->MaxRadius );
   const double HalfBox[3]  = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[4] == BC_FLU_PERIODIC };

#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize arrays
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         OMP_Data  [p][TID][b] = 0.0;
         OMP_Weight[p][TID][b] = 0.0;
         OMP_NCell [p][TID][b] = 0;
      }

//    determine which levels to be considered
      const int lv_min = ( SingleLv < 0 ) ? 0         : SingleLv;
      const int lv_max = ( SingleLv < 0 ) ? TOP_LEVEL : SingleLv;

      for (int lv=lv_min; lv<=lv_max; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );

#        pragma omp for schedule( runtime )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            if ( amr->patch[0][lv][PID]->son != -1 )  continue;

            const real (*FluidPtr)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
#           ifdef GRAVITY
            const real (*PotPtr  )[PS1][PS1]      = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
#           endif

            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

            for (int k=0; k<PS1; k++)  {  double dz = z0 + k*dh;
                                          if ( Periodic[2] ) {
                                             if      ( dz > +HalfBox[2] )  {  dz -= amr->BoxSize[2];  }
                                             else if ( dz < -HalfBox[2] )  {  dz += amr->BoxSize[2];  }
                                          }
            for (int j=0; j<PS1; j++)  {  double dy = y0 + j*dh;
                                          if ( Periodic[1] ) {
                                             if      ( dy > +HalfBox[1] )  {  dy -= amr->BoxSize[1];  }
                                             else if ( dy < -HalfBox[1] )  {  dy += amr->BoxSize[1];  }
                                          }
            for (int i=0; i<PS1; i++)  {  double dx = x0 + i*dh;
                                          if ( Periodic[0] ) {
                                             if      ( dx > +HalfBox[0] )  {  dx -= amr->BoxSize[0];  }
                                             else if ( dx < -HalfBox[0] )  {  dx += amr->BoxSize[0];  }
                                          }

               const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

               if ( r2 < r_max2 )
               {
                  const double r   = sqrt( r2 );
                  const int    bin = ( LogBin ) ? (  (r<dr_min) ? 0 : int( log(r/dr_min)/log(LogBinRatio) ) + 1  )
                                                : int( r/dr_min );
//                prevent from round-off errors
                  if ( bin >= Prof[0]->NBin )   continue;

//                check
#                 ifdef GAMER_DEBUG
                  if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                 endif

                  int ProfIdx = 0;

//                intrinsic fields
                  for (int v=0; v<NFlu; v++)
                  {
                     const real Weight = dv;
                     const int  FluIdx = TVarIdxList_Flu[v];

                     OMP_Data  [ProfIdx][TID][bin] += FluidPtr[FluIdx][k][j][i]*Weight;
                     OMP_Weight[ProfIdx][TID][bin] += Weight;
                     OMP_NCell [ProfIdx][TID][bin] ++;
                     ProfIdx                       ++;
                  }

//                gravitational potential
#                 ifdef GRAVITY
                  if ( PrepPote )
                  {
                     const real Weight = FluidPtr[DENS][k][j][i]*dv;    // weighted by cell mass

                     OMP_Data  [ProfIdx][TID][bin] += PotPtr[k][j][i]*Weight;
                     OMP_Weight[ProfIdx][TID][bin] += Weight;
                     OMP_NCell [ProfIdx][TID][bin] ++;
                     ProfIdx                       ++;
                  }
#                 endif

//                derived fields
#                 if ( MODEL == HYDRO )
                  if ( PrepVelR )
                  {
                     const real Weight = FluidPtr[DENS][k][j][i]*dv;   // weighted by cell mass
                     const real MomR   = ( FluidPtr[MOMX][k][j][i]*dx +
                                           FluidPtr[MOMY][k][j][i]*dy +
                                           FluidPtr[MOMZ][k][j][i]*dz ) / r;

                     OMP_Data  [ProfIdx][TID][bin] += MomR*dv;    // vr*(rho*dv)
                     OMP_Weight[ProfIdx][TID][bin] += Weight;
                     OMP_NCell [ProfIdx][TID][bin] ++;
                     ProfIdx                       ++;
                  }

                  if ( PrepPres )
                  {
                     const real Weight = dv;
#                    ifdef MHD
                     const real EngyB  = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#                    else
                     const real EngyB  = NULL_REAL;
#                    endif
                     const real Pres   = Hydro_GetPressure( FluidPtr[DENS][k][j][i], FluidPtr[MOMX][k][j][i], FluidPtr[MOMY][k][j][i],
                                                            FluidPtr[MOMZ][k][j][i], FluidPtr[ENGY][k][j][i],
                                                            Gamma_m1, false, NULL_REAL, EngyB );

                     OMP_Data  [ProfIdx][TID][bin] += Pres*Weight;
                     OMP_Weight[ProfIdx][TID][bin] += Weight;
                     OMP_NCell [ProfIdx][TID][bin] ++;
                     ProfIdx                       ++;
                  }

                  if ( PrepEint )
                  {
                     const real Weight = dv;
#                    ifdef MHD
                     const real EngyB  = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#                    else
                     const real EngyB  = NULL_REAL;
#                    endif
                     const real Pres   = Hydro_GetPressure( FluidPtr[DENS][k][j][i], FluidPtr[MOMX][k][j][i], FluidPtr[MOMY][k][j][i],
                                                            FluidPtr[MOMZ][k][j][i], FluidPtr[ENGY][k][j][i],
                                                            Gamma_m1, false, NULL_REAL, EngyB );
                     const real Eint   = Pres/Gamma_m1;  // assuming gamma law for now; will be replaced by a general EoS

                     OMP_Data  [ProfIdx][TID][bin] += Eint*Weight;
                     OMP_Weight[ProfIdx][TID][bin] += Weight;
                     OMP_NCell [ProfIdx][TID][bin] ++;
                     ProfIdx                       ++;
                  }
#                 endif // HYDRO
               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=lv_min; lv<=lv_max; lv++)
   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int p=0; p<NProf; p++)
   {
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         Prof[p]->Data  [b]  = OMP_Data  [p][0][b];
         Prof[p]->Weight[b]  = OMP_Weight[p][0][b];
         Prof[p]->NCell [b]  = OMP_NCell [p][0][b];
      }

      for (int t=1; t<NT; t++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         Prof[p]->Data  [b] += OMP_Data  [p][t][b];
         Prof[p]->Weight[b] += OMP_Weight[p][t][b];
         Prof[p]->NCell [b] += OMP_NCell [p][t][b];
      }
   }

// free per-thread arrays
   Aux_DeallocateArray3D( OMP_Data );
   Aux_DeallocateArray3D( OMP_Weight );
   Aux_DeallocateArray3D( OMP_NCell );


// collect data from all ranks (in-place reduction)
#  ifndef SERIAL
   for (int p=0; p<NProf; p++)
   {
      if ( MPI_Rank == 0 )
      {
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->Data,   Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->Weight, Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Prof[p]->NCell , Prof[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }

      else
      {
         MPI_Reduce( Prof[p]->Data,   NULL,            Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Prof[p]->Weight, NULL,            Prof[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Prof[p]->NCell,  NULL,            Prof[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }
   }
#  endif


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   {
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Prof[0]->NBin; b++)
      {
//       skip empty bins since both their data and weight are zero
         if ( Prof[p]->NCell[b] > 0L )    Prof[p]->Data[b] /= Prof[p]->Weight[b];
      }
   }


// broadcast data to all ranks
   for (int p=0; p<NProf; p++)
   {
      MPI_Bcast( Prof[p]->Data,   Prof[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Prof[p]->Weight, Prof[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Prof[p]->NCell,  Prof[p]->NBin, MPI_LONG,   0, MPI_COMM_WORLD );
   }


// remove the empty bins
// --> all ranks do the same work so that no data broadcast is required
   if ( RemoveEmpty )
   {
      for (int b=0; b<Prof[0]->NBin; b++)
      {
         if ( Prof[0]->NCell[b] != 0L )   continue;

//       remove consecutive empty bins at the same time for better performance
         int b_up;
         for (b_up=b+1; b_up<Prof[0]->NBin; b_up++)
            if ( Prof[0]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (b_up=b+stride; b_up<Prof[0]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            for (int p=0; p<NProf; p++)
            {
               Prof[p]->Radius[b_up_ms] = Prof[p]->Radius[b_up];
               Prof[p]->Data  [b_up_ms] = Prof[p]->Data  [b_up];
               Prof[p]->Weight[b_up_ms] = Prof[p]->Weight[b_up];
               Prof[p]->NCell [b_up_ms] = Prof[p]->NCell [b_up];
            }
         }

//       reset the total number of bins
         for (int p=0; p<NProf; p++)
            Prof[p]->NBin -= stride;
      } // for (int b=0; b<Prof->NBin; b++)

//    update the maximum radius since the last bin may have been removed
      for (int p=0; p<NProf; p++)
      {
         const int LastBin = Prof[p]->NBin-1;

         Prof[p]->MaxRadius = ( LogBin ) ? Prof[p]->Radius[LastBin]*sqrt( LogBinRatio )
                                         : Prof[p]->Radius[LastBin] + 0.5*dr_min;
      }
   } // if ( RemoveEmpty )

} // FUNCTION : Aux_ComputeProfile
