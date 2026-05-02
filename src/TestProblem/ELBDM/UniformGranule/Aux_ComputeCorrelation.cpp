#include "GAMER.h"

extern void SetTempIntPara( const int lv, const int Sg0, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );


//-------------------------------------------------------------------------------------------------------
// Function    :  InterpolateMeanAndStd
// Description :  Interpolate the initial density for a given radius
//
// Note        :  1. Use only linear interpolation
//                2. When invoking in Aux_ComputeCorrelation function, the prof_init is always assumed to be in linear bin, so the bin_index
//                   passed-in is also estimated based on linear bin; to use this function for non-linear bin, user need to compute bin_index
//                   accordingly
//
// Parameter   :  mean_inter : interpolated mean value for all the target fields at r
//                std_inter  : interpolated standard deviation for all the target fields at r
//                prof_init  : Profile_t object array for storing the mean and standard deviation quantities used for interpolation
//                NProf      : Number of Profile_t objects in prof_init
//                bin_index  : estimated bin index for a given target radius r
//                r          : target radius
//
// Return      :  mean_inter, std_inter
//-------------------------------------------------------------------------------------------------------
void InterpolateMeanAndStd(real *mean_inter, real *std_inter, const Profile_t *prof_init[], const int NProf, const int bin_index, const double r)
{
   double delta_r, x;
   if ( r > prof_init[0]->Radius[bin_index] )
   {
      int bin_index_right = bin_index+1;
      delta_r             = prof_init[0]->Radius[bin_index_right] - prof_init[0]->Radius[bin_index];
      x                   = (r - prof_init[0]->Radius[bin_index]) / delta_r;
//    check x
      if (x<(real)(0.0))
         Aux_Error( ERROR_INFO, "x (%14.7e) < 0.0 !! index = %d ; r = %14.7e ; left-hand point = %14.7e ; right-hand point = %14.7e \n", x , bin_index, r, prof_init[0]->Radius[bin_index], prof_init[0]->Radius[bin_index_right] );
      else if (x>(real)(1.0))
         Aux_Error( ERROR_INFO, "x (%14.7e) > 1.0 !! index = %d ; r = %14.7e ; left-hand point = %14.7e ; right-hand point = %14.7e \n", x , bin_index, r, prof_init[0]->Radius[bin_index], prof_init[0]->Radius[bin_index_right] );
//    interpolate
      for  (int i=0; i<NProf; i++)
      {
         mean_inter[i] = prof_init[i]->Data      [bin_index]*(real)(1.-x) + prof_init[i]->Data      [bin_index_right]*(real)x;
         std_inter[i]  = prof_init[i]->Data_Sigma[bin_index]*(real)(1.-x) + prof_init[i]->Data_Sigma[bin_index_right]*(real)x;
      }
   }
   else
   {
//    no left hand side bin, no interpolation
      if (bin_index==0)
      {
         for (int i=0; i<NProf; i++)
         {
            mean_inter[i] = prof_init[i]->Data      [bin_index];
            std_inter[i]  = prof_init[i]->Data_Sigma[bin_index];
         }
      }
      else
      {
         int bin_index_left  = bin_index-1;
         delta_r             = prof_init[0]->Radius[bin_index] - prof_init[0]->Radius[bin_index_left];
         x                   = (r - prof_init[0]->Radius[bin_index_left]) / delta_r;
//       check x
         if (x<(real)(0.0))
            Aux_Error( ERROR_INFO, "x (%14.7e) < 0.0 !! index = %d ; r = %14.7e ; left-hand point = %14.7e ; right-hand point = %14.7e \n", x , bin_index, r, prof_init[0]->Radius[bin_index_left], prof_init[0]->Radius[bin_index] );
         else if (x>(real)(1.0))
            Aux_Error( ERROR_INFO, "x (%14.7e) > 1.0 !! index = %d ; r = %14.7e ; left-hand point = %14.7e ; right-hand point = %14.7e \n", x , bin_index, r, prof_init[0]->Radius[bin_index_left], prof_init[0]->Radius[bin_index] );
//       interpolate
         for  (int i=0; i<NProf; i++)
         {
            mean_inter[i] = prof_init[i]->Data      [bin_index_left]*(real)(1.-x) + prof_init[i]->Data      [bin_index]*(real)x;
            std_inter[i]  = prof_init[i]->Data_Sigma[bin_index_left]*(real)(1.-x) + prof_init[i]->Data_Sigma[bin_index]*(real)x;
         }
      }
   }
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_ComputeCorrelation
// Description :  Compute the average radial correlation profile of target field(s)
//
// Note        :  1. Results will be stored in the input "Correlation" object
//                   --> Correlation->Radius[]: Radial coordinate at each bin
//                       Correlation->Data  []: Correlation profile data at each bin
//                       Correlation->Weight[]: Total weighting at each bin
//                       Correlation->NCell []: Number of cells at each bin
//                       Correlation->NBin    : Total number of bins
//                   --> See the "Profile_t" structure defined in "include/Profile.h" for details
//                   --> These arrays will be free'd when deleting "Correlation"
//                2. Maximum radius adopted when actually computing the profile may be larger than the input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//                3. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same profile data after invoking this function
//                4. Use cell volume as the weighting of each cell
//                   --> Will support other weighting functions in the future
//                5. Currently only support computing _DENS field
//                   --> If multiple fields are supported in the future, the order of fields to be returned follows TVarBitIdx[]
//
// Parameter   :  Correlation : Profile_t object array to store the correlation function
//                prof_init   : Profile_t object array for storing the mean and standard deviation quantities for calculating correlation function
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
//                TVarBitIdx  : Bitwise indices of target variables for computing the profiles
//                              --> Supported indices (defined in Macro.h):
//                                     ELBDM : _DENS
//                              --> For a passive scalar with an integer field index FieldIdx returned by AddField(),
//                                  one can convert it to a bitwise field index by BIDX(FieldIdx)
//                NProf       : Number of Profile_t objects in Correlation
//                Min/MaxLv   : Consider patches on levels from MinLv to MaxLv
//                PatchType   : Only consider patches of the specified type
//                              --> Supported types: PATCH_LEAF, PATCH_NONLEAF, PATCH_BOTH, PATCH_LEAF_PLUS_MAXNONLEAF
//                              --> PATCH_LEAF_PLUS_MAXNONLEAF includes leaf patches on all target levels
//                                  (i.e., MinLv ~ MaxLv) and non-leaf patches only on MaxLv
//                PrepTime    : Target physical time to prepare data
//                              --> If PrepTime<0, turn off temporal interpolation and always use the most recent data
//
// Example     :  const double      Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double      MaxRadius      = 0.5*amr->BoxSize[0];
//                const double      MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool        LogBin         = true;
//                const double      LogBinRatio    = 1.25;
//                const bool        RemoveEmptyBin = true;
//                const long        TVarBitIdx[]   = { _DENS };
//                const int         NProf          = 1;
//                const int         MinLv          = 0;
//                const int         MaxLv          = MAX_LEVEL;
//                const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
//                const double      PrepTime       = -1.0;
//
//                Profile_t Correlation_Dens;
//                Profile_t *Correlation[] = { &Correlation_Dens };
//
//                Aux_ComputeCorrelation( Correlation, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
//                                        TVarBitIdx, NProf, MinLv, MaxLv, PatchType, PrepTime );
//
//                if ( MPI_Rank == 0 )
//                {
//                   for (int p=0; p<NProf; p++)
//                   {
//                      char Filename[MAX_STRING];
//                      sprintf( Filename, "Correlation_function%d.txt", p );
//                      FILE *File = fopen( Filename, "w" );
//                      fprintf( File, "#%19s  %21s  %21s  %10s\n", "Radius", "Data", "Weight", "Cells" );
//                      for (int b=0; b<Correlation[p]->NBin; b++)
//                         fprintf( File, "%20.14e  %21.14e  %21.14e  %10ld\n",
//                                  Correlation[p]->Radius[b], Correlation[p]->Data[b], Correlation[p]->Weight[b], Correlation[p]->NCell[b] );
//                      fclose( File );
//                   }
//                }
//
// Return      :  Correlation
//-------------------------------------------------------------------------------------------------------
void Aux_ComputeCorrelation( Profile_t *Correlation[], const Profile_t *prof_init[], const double Center[],
                             const double r_max_input, const double dr_min, const bool LogBin, const double LogBinRatio,
                             const bool RemoveEmpty, const long TVarBitIdx[], const int NProf, const int MinLv, const int MaxLv,
                             const PatchType_t PatchType, const double PrepTime, const double dr_min_prof)
{

// check
#  ifdef GAMER_DEBUG
   if ( r_max_input <= 0.0 )
      Aux_Error( ERROR_INFO, "r_max_input (%14.7e) <= 0.0 !!\n", r_max_input );

   if ( dr_min <= 0.0 )
      Aux_Error( ERROR_INFO, "dr_min (%14.7e) <= 0.0 !!\n", dr_min );

   if ( LogBin  &&  LogBinRatio <= 1.0 )
      Aux_Error( ERROR_INFO, "LogBinRatio (%14.7e) <= 1.0 !!\n", LogBinRatio );

   if ( MinLv < 0  ||  MinLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MinLv (%d) !!\n", MinLv );

   if ( MaxLv < 0  ||  MaxLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MaxLv (%d) !!\n", MaxLv );

   if ( MinLv > MaxLv )
      Aux_Error( ERROR_INFO, "MinLv (%d) > MaxLv (%d) !!\n", MinLv, MaxLv );

   if ( PatchType != PATCH_LEAF  &&  PatchType != PATCH_NONLEAF  &&
        PatchType != PATCH_BOTH  &&  PatchType != PATCH_LEAF_PLUS_MAXNONLEAF )
      Aux_Error( ERROR_INFO, "incorrect PatchType (%d) !!\n", PatchType );
#  endif
   if ( NProf != NCOMP_PASSIVE )
      Aux_Error( ERROR_INFO, "NProf(%d) != NCOMP_PASSIVE(%d) !! Currently only support NProf = NCOMP_PASSIVE for computing correlation !!\n", NProf, NCOMP_PASSIVE );
   if ( TVarBitIdx[0] != _DENS )
      Aux_Error( ERROR_INFO, "TVarBitIdx[0](%ld) != _DENS(%ld) !! Currently only support TVarBitIdx[0] = _DENS for computing correlation !!\n", TVarBitIdx[0], _DENS );
#  if ( MODEL == HYDRO )
      Aux_Error( ERROR_INFO, "Does not support HDRDO for computing correlation function yet!!\n" );
#  endif


// precompute the integer indices of intrinsic fluid fields for better performance
   const int IdxUndef = -1;
   int TFluIntIdx[NProf];

   for (int p=0; p<NProf; p++)
   {
      TFluIntIdx[p] = IdxUndef;

      for (int v=0; v<NCOMP_TOTAL; v++)
         if ( TVarBitIdx[p] & (1L<<v) )   TFluIntIdx[p] = v;
   }


// check whether _POTE is in TVarBitIdx since the potential array may have not been computed during initialization
#  ifdef GRAVITY
   bool InclPot = false;

   for (int p=0; p<NProf; p++)
      if ( TVarBitIdx[p] & _POTE )   InclPot = true;
#  endif

// check whether phase field is accessed in hybrid scheme
// currently computing the profile of the phase field is not supported
#  if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID)
   bool UsePhaseStub = false;
   for (int p=0; p<NProf; p++)
      if ( TVarBitIdx[p] & _PHAS ||  TVarBitIdx[p] & _STUB ) UsePhaseStub = true;

   if ( UsePhaseStub )
      for (int lv=MinLv; lv<=MaxLv; lv++)
         if ( !amr->use_wave_flag[lv] )
            Aux_Error( ERROR_INFO, "Retrieving PHAS and STUB to compute profile in hybrid scheme is not supported !!\n" );
#  endif // #  if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID)

// initialize the profile objects
   for (int p=0; p<NProf; p++)
   {
//    get the total number of radial bins and the corresponding maximum radius
      if ( LogBin )
      {
//       MaxRadius will be smaller than r_max_input if Correlation[p]->NBin = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 1;
//                 will be greater than r_max_input if Correlation[p]->NBin = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 2;
         Correlation[p]->NBin      = int( log(r_max_input/dr_min)/log(LogBinRatio) ) + 1;
         Correlation[p]->MaxRadius = dr_min*pow( LogBinRatio, Correlation[p]->NBin-1 );
      }

      else // linear bin
      {
         Correlation[p]->NBin      = (int)ceil( r_max_input / dr_min );
         Correlation[p]->MaxRadius = dr_min*Correlation[p]->NBin;
      }

      if ( Correlation[p]->MaxRadius > prof_init[p]->MaxRadius )  Aux_Error( ERROR_INFO, "Correlation[%d]->MaxRadius ( %14.7e ) > Profile[%d]->MaxRadius ( %14.7e ) !!\n", p, Correlation[p]->MaxRadius, p, prof_init[p]->MaxRadius );


//    record profile parameters
      for (int d=0; d<3; d++)    Correlation[p]->Center[d] = Center[d];

      Correlation[p]->LogBin = LogBin;

      if ( LogBin )  Correlation[p]->LogBinRatio = LogBinRatio;


//    allocate all member arrays of Correlation
      Correlation[p]->AllocateMemory();


//    record radial coordinates
      if ( LogBin )
         for (int b=0; b<Correlation[0]->NBin; b++)    Correlation[p]->Radius[b] = dr_min*pow( LogBinRatio, b-0.5 );
      else
         for (int b=0; b<Correlation[0]->NBin; b++)    Correlation[p]->Radius[b] = (b+0.5)*dr_min;

   } // for (int p=0; p<NProf; p++)


// allocate memory for the per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double ***OMP_Data=NULL, ***OMP_Weight=NULL;
   long   ***OMP_NCell=NULL;

   Aux_AllocateArray3D( OMP_Data,   NProf, NT, Correlation[0]->NBin );
   Aux_AllocateArray3D( OMP_Weight, NProf, NT, Correlation[0]->NBin );
   Aux_AllocateArray3D( OMP_NCell,  NProf, NT, Correlation[0]->NBin );


// collect profile data in this rank
   const double r_max2      = SQR( Correlation[0]->MaxRadius );
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
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         OMP_Data  [p][TID][b] = 0.0;
         OMP_Weight[p][TID][b] = 0.0;
         OMP_NCell [p][TID][b] = 0;
      }

//    allocate passive scalar arrays
      real *Passive      = new real [NCOMP_PASSIVE];

//    loop over all target levels
      for (int lv=MinLv; lv<=MaxLv; lv++)
      {
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );


//       determine temporal interpolation parameters
         bool FluIntTime = false;
         int  FluSg      = amr->FluSg[lv];
         int  FluSg_IntT;
         real FluWeighting, FluWeighting_IntT;

         if ( PrepTime >= 0.0 )
         {
//          fluid
            const int FluSg0 = amr->FluSg[lv];
            SetTempIntPara( lv, FluSg0, PrepTime, amr->FluSgTime[lv][FluSg0], amr->FluSgTime[lv][1-FluSg0],
                            FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );
         }


//       use the "static" schedule for reproducibility
#        pragma omp for schedule( static )
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
//          skip untargeted patches
            if ( amr->patch[0][lv][PID]->son != -1 )
            {
               if ( PatchType == PATCH_LEAF )                                    continue;
               if ( PatchType == PATCH_LEAF_PLUS_MAXNONLEAF  &&  lv != MaxLv )   continue;
            }

            else
            {
               if ( PatchType == PATCH_NONLEAF )                                 continue;
            }

            const real (*FluidPtr)[PS1][PS1][PS1] = amr->patch[ FluSg ][lv][PID]->fluid;

//          pointer for temporal interpolation
            const real (*FluidPtr_IntT)[PS1][PS1][PS1] = ( FluIntTime ) ? amr->patch[ FluSg_IntT ][lv][PID]->fluid : NULL;
            const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
            const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
            const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

            for (int k=0; k<PS1; k++)  {  double dz = z0 + k*dh;
                                          if ( Periodic[2] )
                                          {
                                             if      ( dz > +HalfBox[2] )  {  dz -= amr->BoxSize[2];  }
                                             else if ( dz < -HalfBox[2] )  {  dz += amr->BoxSize[2];  }
                                          }
            for (int j=0; j<PS1; j++)  {  double dy = y0 + j*dh;
                                          if ( Periodic[1] )
                                          {
                                             if      ( dy > +HalfBox[1] )  {  dy -= amr->BoxSize[1];  }
                                             else if ( dy < -HalfBox[1] )  {  dy += amr->BoxSize[1];  }
                                          }
            for (int i=0; i<PS1; i++)  {  double dx = x0 + i*dh;
                                          if ( Periodic[0] )
                                          {
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
                  if ( bin >= Correlation[0]->NBin )   continue;

//                check
#                 ifdef GAMER_DEBUG
                  if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                 endif

//                interpolate to get mean value at r
                  real mean_value[NProf], std_value[NProf];
//                ****find corresponding bin index in density profile, which always uses linear bin!!*****
                  const int    bin_prof = int (r/dr_min_prof);
                  InterpolateMeanAndStd( mean_value, std_value, prof_init, NProf, bin_prof, r );

//                prepare passive scalars (for better sustainability, always do it even when unnecessary)
                  for (int v_out=0; v_out<NCOMP_PASSIVE; v_out++)
                  {
                     const int v_in = v_out + NCOMP_FLUID;

                     Passive[v_out] = FluidPtr[v_in][k][j][i];
                  }

                  for (int p=0; p<NProf; p++)
                  {
//                   intrinsic fluid fields
                     if ( TFluIntIdx[p] != IdxUndef )
                     {
                        const real Weight = dv;
                        real delta  = ( FluIntTime )
                                          ? ( FluWeighting     *FluidPtr     [ TFluIntIdx[p] ][k][j][i]
                                            + FluWeighting_IntT*FluidPtr_IntT[ TFluIntIdx[p] ][k][j][i] )
                                          :                     FluidPtr     [ TFluIntIdx[p] ][k][j][i]  ;
                        real delta_passive = Passive[p];
                        delta         -= mean_value[p];
                        delta_passive -= mean_value[p];

                        if (std_value[p]>(real)0.)
                           OMP_Data[p][TID][bin] += delta*delta_passive/std_value[p]/std_value[p]*Weight;
                        else
                           OMP_Data[p][TID][bin] += delta*delta_passive/mean_value[p]/mean_value[p]*Weight;

                        OMP_Weight[p][TID][bin] += Weight;
                        OMP_NCell [p][TID][bin] ++;
                     }
                  } // for (int p=0; p<NProf; p++)
               } // if ( r2 < r_max2 )
            }}} // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=MinLv; lv<=MaxLv; lv++)

      delete [] Passive;
      Passive = NULL;

   } // OpenMP parallel region


// sum over all OpenMP threads
   for (int p=0; p<NProf; p++)
   {
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         Correlation[p]->Data  [b]  = OMP_Data  [p][0][b];
         Correlation[p]->Weight[b]  = OMP_Weight[p][0][b];
         Correlation[p]->NCell [b]  = OMP_NCell [p][0][b];
      }

      for (int t=1; t<NT; t++)
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         Correlation[p]->Data  [b] += OMP_Data  [p][t][b];
         Correlation[p]->Weight[b] += OMP_Weight[p][t][b];
         Correlation[p]->NCell [b] += OMP_NCell [p][t][b];
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
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->Data,   Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->Weight, Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( MPI_IN_PLACE,    Correlation[p]->NCell , Correlation[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }

      else
      {
         MPI_Reduce( Correlation[p]->Data,   NULL,            Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Correlation[p]->Weight, NULL,            Correlation[p]->NBin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
         MPI_Reduce( Correlation[p]->NCell,  NULL,            Correlation[p]->NBin, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
      }
   }
#  endif


// compute profile by the root rank
   if ( MPI_Rank == 0 )
   {
      for (int p=0; p<NProf; p++)
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
//       skip empty bins since both their data and weight are zero
         if ( Correlation[p]->NCell[b] > 0L )    Correlation[p]->Data[b] /= Correlation[p]->Weight[b];
      }
   }


// broadcast data to all ranks
   for (int p=0; p<NProf; p++)
   {
      MPI_Bcast( Correlation[p]->Data,   Correlation[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Correlation[p]->Weight, Correlation[p]->NBin, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Bcast( Correlation[p]->NCell,  Correlation[p]->NBin, MPI_LONG,   0, MPI_COMM_WORLD );
   }


// remove the empty bins
// --> all ranks do the same work so that no data broadcast is required
   if ( RemoveEmpty )
   {
      for (int b=0; b<Correlation[0]->NBin; b++)
      {
         if ( Correlation[0]->NCell[b] != 0L )   continue;

//       remove consecutive empty bins at the same time for better performance
         int b_up;
         for (b_up=b+1; b_up<Correlation[0]->NBin; b_up++)
            if ( Correlation[0]->NCell[b_up] != 0L )   break;

         const int stride = b_up - b;

         for (b_up=b+stride; b_up<Correlation[0]->NBin; b_up++)
         {
            const int b_up_ms = b_up - stride;

            for (int p=0; p<NProf; p++)
            {
               Correlation[p]->Radius[b_up_ms] = Correlation[p]->Radius[b_up];
               Correlation[p]->Data  [b_up_ms] = Correlation[p]->Data  [b_up];
               Correlation[p]->Weight[b_up_ms] = Correlation[p]->Weight[b_up];
               Correlation[p]->NCell [b_up_ms] = Correlation[p]->NCell [b_up];
            }
         }

//       reset the total number of bins
         for (int p=0; p<NProf; p++)
            Correlation[p]->NBin -= stride;
      } // for (int b=0; b<Correlation->NBin; b++)

//    update the maximum radius since the last bin may have been removed
      for (int p=0; p<NProf; p++)
      {
         const int LastBin = Correlation[p]->NBin-1;

         Correlation[p]->MaxRadius = ( LogBin ) ? Correlation[p]->Radius[LastBin]*sqrt( LogBinRatio )
                                                : Correlation[p]->Radius[LastBin] + 0.5*dr_min;
      }
   } // if ( RemoveEmpty )

} // FUNCTION : Aux_ComputeCorrelation
