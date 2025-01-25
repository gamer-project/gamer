#include "GAMER.h"

extern void SetTempIntPara( const int lv, const int Sg0, const double PrepTime, const double Time0, const double Time1,
                            bool &IntTime, int &Sg, int &Sg_IntT, real &Weighting, real &Weighting_IntT );




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
//                2. The exact maximum radius adopted may be slightly larger than the input "r_max"
//                   --> Because "r_max" in general does not coincide with the right edge of the maximum bin
//                3. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same profile data after invoking this function
//                4. Weighting of each cell:
//                      Cell mass  : gas velocity, gravitational potential
//                      Cell volume: other fields
//                   --> Will support more weighting fields in the future
//                5. Support computing multiple fields
//                   --> The order of fields to be returned follows TVarBitIdx[]
//                6. This routine is thread-unsafe when the temporal interpolation set by PrepTime and OPT__INT_TIME
//                   are inconsistent with each other
//                   --> But it shouldn't be a big issue since this routine itself has been parallelized with OpenMP
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
//                TVarBitIdx  : Bitwise indices of target variables for computing the profiles
//                              --> Supported indices (defined in Macro.h):
//                                     HYDRO        : _DENS, _MOMX, _MOMY, _MOMZ, _ENGY, _VELX, _VELY, _VELZ, _VELR,
//                                                    _PRES, _TEMP, _ENTR, _EINT
//                                                    [, _DUAL, _CRAY, _POTE, __MAGX_CC, _MAGY_CC, _MAGZ_CC, _MAGE_CC]
//                                     ELBDM_WAVE   : _DENS, _REAL, _IMAG [, _POTE]
//                                     ELBDM_HYBRID : _DENS, _PHAS [, _POTE]
//                              --> All fields supported by Prepare_PatchData() are also supported here
//                              --> For a passive scalar with an integer field index FieldIdx returned by AddField(),
//                                  one can convert it to a bitwise field index by BIDX(FieldIdx)
//                NProf       : Number of Profile_t objects in Prof
//                Min/MaxLv   : Consider patches on levels from MinLv to MaxLv
//                PatchType   : Only consider patches of the specified type
//                              --> Supported types: PATCH_LEAF, PATCH_NONLEAF, PATCH_BOTH, PATCH_LEAF_PLUS_MAXNONLEAF
//                              --> PATCH_LEAF_PLUS_MAXNONLEAF includes leaf patches on all target levels
//                                  (i.e., MinLv ~ MaxLv) and non-leaf patches only on MaxLv
//                PrepTimeIn  : Target physical time to prepare data
//                              --> If PrepTimeIn<0, turn off temporal interpolation and always use the most recent data
//
// Example     :  const double      Center[3]      = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double      MaxRadius      = 0.5*amr->BoxSize[0];
//                const double      MinBinSize     = amr->dh[MAX_LEVEL];
//                const bool        LogBin         = true;
//                const double      LogBinRatio    = 1.25;
//                const bool        RemoveEmptyBin = true;
//                const long        TVar[]         = { _DENS, _PRES };
//                const int         NProf          = 2;
//                const int         MinLv          = 0;
//                const int         MaxLv          = MAX_LEVEL;
//                const PatchType_t PatchType      = PATCH_LEAF_PLUS_MAXNONLEAF;
//                const double      PrepTime       = -1.0;
//
//                Profile_t Prof_Dens, Prof_Pres;
//                Profile_t *Prof[] = { &Prof_Dens, &Prof_Pres };
//
//                Aux_ComputeProfile( Prof, Center, MaxRadius, MinBinSize, LogBin, LogBinRatio, RemoveEmptyBin,
//                                    TVar, NProf, MinLv, MaxLv, PatchType, PrepTime );
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
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const long TVarBitIdx[],
                         const int NProf, const int MinLv, const int MaxLv, const PatchType_t PatchType,
                         const double PrepTimeIn )
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

#  ifdef OPENMP
   if (  omp_in_parallel()  &&  ( (PrepTimeIn >= 0.0 && !OPT__INT_TIME) ||
                                  (PrepTimeIn <  0.0 &&  OPT__INT_TIME) )  )
      Aux_Error( ERROR_INFO, "this routine is thread-unsafe when the temporal interpolation set "
                             "by PrepTimeIn and OPT__INT_TIME are inconsistent !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG


// list all supported fields
// --> all fields supported by Prepare_PatchData() should be supported here
   long SupportedFields = ( _TOTAL | _DERIVED );
#  ifdef GRAVITY
   SupportedFields |= _POTE;
#  endif
#  ifdef MASSIVE_PARTICLES
   SupportedFields |= _PAR_DENS;
   SupportedFields |= _TOTAL_DENS;
#  endif

   for (int p=0; p<NProf; p++) {
      if ( TVarBitIdx[p] & ~SupportedFields )
         Aux_Error( ERROR_INFO, "unsupported field (TVarBitIdx[%d] = %ld) !!\n", p, TVarBitIdx[p] );
   }


// record whether particle density is requested
#  ifdef MASSIVE_PARTICLES
   bool NeedPar = false;
   for (int p=0; p<NProf; p++) {
      if ( TVarBitIdx[p] == _PAR_DENS  ||  TVarBitIdx[p] == _TOTAL_DENS ) {
         NeedPar = true;
         break;
      }
   }
#  endif


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

// initialize profile arrays
   for (int p=0; p<NProf; p++)
   for (int t=0; t<NT; t++)
   for (int b=0; b<Prof[0]->NBin; b++)
   {
      OMP_Data  [p][t][b] = 0.0;
      OMP_Weight[p][t][b] = 0.0;
      OMP_NCell [p][t][b] = 0;
   }

   real (*Patch_Data)[8][PS1][PS1][PS1] = new real [NT][8][PS1][PS1][PS1];  // field data of each cell
   int  (*Patch_Bin )[8][PS1][PS1][PS1] = new int  [NT][8][PS1][PS1][PS1];  // radial bin of each cell


// set global constants
   const double r_max2         = SQR( Prof[0]->MaxRadius );
   const double HalfBox[3]     = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3]    = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                   OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                   OPT__BC_FLU[4] == BC_FLU_PERIODIC };
   const int    CellSkip       = -1;
   const int    WeightByVolume = 1;
   const int    WeightByMass   = 2;


// temporarily overwrite OPT__INT_TIME
// --> necessary because SetTempIntPara() called by Prepare_PatchData() relies on OPT__INT_TIME
// --> must restore it before exiting this routine
// --> note that modifying OPT__INT_TIME renders this routine thread-unsafe
//###REVISE: make temporal interpolation a function parameter in Prepare_PatchData() to solve this thread-safety issue
   const bool IntTimeBackup = OPT__INT_TIME;
   OPT__INT_TIME = ( PrepTimeIn >= 0.0 ) ? true : false;


// loop over all target levels
   for (int lv=MinLv; lv<=MaxLv; lv++)
   {
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );


//    determine the temporal interpolation parameters
//    --> mainly for computing cell mass for weighting; Prepare_PatchData() needs PrepTime
      const int    FluSg0 = amr->FluSg[lv];
      const double PrepTime = ( PrepTimeIn >= 0.0 ) ? PrepTimeIn : amr->FluSgTime[lv][FluSg0];

      bool FluIntTime;
      int  FluSg, FluSg_IntT;
      real FluWeighting, FluWeighting_IntT;

      SetTempIntPara( lv, FluSg0, PrepTime, amr->FluSgTime[lv][FluSg0], amr->FluSgTime[lv][1-FluSg0],
                      FluIntTime, FluSg, FluSg_IntT, FluWeighting, FluWeighting_IntT );


//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef MASSIVE_PARTICLES
      const bool TimingSendPar_No = false;
      const bool JustCountNPar_No = false;
#     ifdef LOAD_BALANCE
      const bool PredictPos       = amr->Par->PredictPos;
      const bool SibBufPatch      = true;
      const bool FaSibBufPatch    = true;
#     else
      const bool PredictPos       = false;
      const bool SibBufPatch      = NULL_BOOL;
      const bool FaSibBufPatch    = NULL_BOOL;
#     endif

      if ( NeedPar )
      {
//       these two routines should NOT be put inside an OpenMP parallel region
         Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictPos,
                                       PrepTime, SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

         Prepare_PatchData_InitParticleDensityArray( lv, PrepTime );
      } // if ( NeedPar )
#     endif // #ifdef MASSIVE_PARTICLES


//    different OpenMP threads and MPI processes first compute profiles independently
//    --> their data will be combined later
#     pragma omp parallel
      {
#        ifdef OPENMP
         const int TID = omp_get_thread_num();
#        else
         const int TID = 0;
#        endif

//       use the "static" schedule for reproducibility
#        pragma omp for schedule( static )
         for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
         {
//          skip untargeted patches
            bool SkipPatch[8], SkipPatchGroup=true;

            for (int LocalID=0; LocalID<8; LocalID++)
            {
               const int PID = PID0 + LocalID;
               SkipPatch[LocalID] = false;

               if ( amr->patch[0][lv][PID]->son != -1 )
               {
                  if ( PatchType == PATCH_LEAF )                                    SkipPatch[LocalID] = true;
                  if ( PatchType == PATCH_LEAF_PLUS_MAXNONLEAF  &&  lv != MaxLv )   SkipPatch[LocalID] = true;
               }

               else
               {
                  if ( PatchType == PATCH_NONLEAF )                                 SkipPatch[LocalID] = true;
               }

               if ( ! SkipPatch[LocalID] )   SkipPatchGroup = false;
            } // for (int LocalID=0; LocalID<8; LocalID++)

            if ( SkipPatchGroup )   continue;


//          store the radial bin associated with each cell
//          --> do it before looping over all target fields to avoid redundant calculations
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               if ( SkipPatch[LocalID] )  continue;

               const int    PID = PID0 + LocalID;
               const double x0  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh - Center[0];
               const double y0  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh - Center[1];
               const double z0  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh - Center[2];

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
//                   prevent from round-off errors
                     if ( bin >= Prof[0]->NBin )   Patch_Bin[TID][LocalID][k][j][i] = CellSkip;
                     else                          Patch_Bin[TID][LocalID][k][j][i] = bin;

//                   check
#                    ifdef GAMER_DEBUG
                     if ( bin < 0 )    Aux_Error( ERROR_INFO, "bin (%d) < 0 !!\n", bin );
#                    endif
                  }

                  else
                     Patch_Bin[TID][LocalID][k][j][i] = CellSkip;
               }}} // i,j,k
            } // for (int LocalID=0; LocalID<8; LocalID++)


//          compute one field at a time
            for (int p=0; p<NProf; p++)
            {
//             collect the data of the target field
               switch ( TVarBitIdx[p] )
               {
//                _VELR is currently not supported by Prepare_PatchData()
#                 ifdef _VELR
                  case _VELR:
                     for (int LocalID=0; LocalID<8; LocalID++)
                     {
                        if ( SkipPatch[LocalID] )  continue;

                        const int PID = PID0 + LocalID;
                        const real (*FluidPtr     )[PS1][PS1][PS1] =                  amr->patch[ FluSg      ][lv][PID]->fluid;
                        const real (*FluidPtr_IntT)[PS1][PS1][PS1] = ( FluIntTime ) ? amr->patch[ FluSg_IntT ][lv][PID]->fluid : NULL;

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

                           if ( Patch_Bin[TID][LocalID][k][j][i] == CellSkip )   continue;

                           const double r        = sqrt( SQR(dx) + SQR(dy) + SQR(dz) );
                           const real _Dens      =                  (real)1.0 / FluidPtr     [DENS][k][j][i];
                           const real _Dens_IntT = ( FluIntTime ) ? (real)1.0 / FluidPtr_IntT[DENS][k][j][i] : NULL_REAL;

                           real VelR;
                           if ( r == 0.0 ) {
                              VelR = (real)0.0;    // take care of the corner case where the profile center coincides with a cell center
                           }

                           else {
                              VelR = ( FluIntTime )
                                   ? ( FluWeighting     *( FluidPtr     [MOMX][k][j][i]*dx +
                                                           FluidPtr     [MOMY][k][j][i]*dy +
                                                           FluidPtr     [MOMZ][k][j][i]*dz )*_Dens
                                     + FluWeighting_IntT*( FluidPtr_IntT[MOMX][k][j][i]*dx +
                                                           FluidPtr_IntT[MOMY][k][j][i]*dy +
                                                           FluidPtr_IntT[MOMZ][k][j][i]*dz )*_Dens_IntT ) / r
                                   :                     ( FluidPtr     [MOMX][k][j][i]*dx +
                                                           FluidPtr     [MOMY][k][j][i]*dy +
                                                           FluidPtr     [MOMZ][k][j][i]*dz )*_Dens / r;
                           }

                           Patch_Data[TID][LocalID][k][j][i] = VelR;
                        }}} // i,j,k
                     } // for (int LocalID=0; LocalID<8; LocalID++)
                  break; // _VELR
#                 endif // #ifdef _VELR

                  default:
                     const int  NGhost             = 0;
                     const int  NPG                = 1;
                     const bool IntPhase_No        = false;
                     const real MinDens_No         = -1.0;
                     const real MinPres_No         = -1.0;
                     const real MinTemp_No         = -1.0;
                     const real MinEntr_No         = -1.0;
                     const bool DE_Consistency_Yes = true;

                     Prepare_PatchData( lv, PrepTime, &Patch_Data[TID][0][0][0][0], NULL, NGhost, NPG, &PID0,
                                        TVarBitIdx[p], _NONE, INT_NONE, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No,
                                        OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No,
                                        DE_Consistency_Yes );
                  break; // default
               } // switch ( TVarBitIdx[p] )


//             set the weight field
//###REVISE: allow users to choose the weight field
               int WeightField=-1;

               switch ( TVarBitIdx[p] )
               {
#                 ifdef _VELX
                  case _VELX : WeightField = WeightByMass;     break;
#                 endif
#                 ifdef _VELY
                  case _VELY : WeightField = WeightByMass;     break;
#                 endif
#                 ifdef _VELZ
                  case _VELZ : WeightField = WeightByMass;     break;
#                 endif
#                 ifdef _VELR
                  case _VELR : WeightField = WeightByMass;     break;
#                 endif
#                 ifdef _POTE
                  case _POTE : WeightField = WeightByMass;     break;
#                 endif
                  default    : WeightField = WeightByVolume;   break;
               } // switch ( TVarBitIdx[p] )


//             compute the radial profile
               for (int LocalID=0; LocalID<8; LocalID++)
               {
                  if ( SkipPatch[LocalID] )  continue;

                  const int PID = PID0 + LocalID;
                  const real (*DensPtr     )[PS1][PS1] =                  amr->patch[ FluSg      ][lv][PID]->fluid[DENS];
                  const real (*DensPtr_IntT)[PS1][PS1] = ( FluIntTime ) ? amr->patch[ FluSg_IntT ][lv][PID]->fluid[DENS] : NULL;

                  for (int k=0; k<PS1; k++)
                  for (int j=0; j<PS1; j++)
                  for (int i=0; i<PS1; i++)
                  {
                     if ( Patch_Bin[TID][LocalID][k][j][i] == CellSkip )   continue;

//                   compute the weight
                     real Weight;
                     switch ( WeightField )
                     {
                        case WeightByMass   :   Weight = ( FluIntTime )
                                                       ? ( FluWeighting     *DensPtr     [k][j][i]
                                                         + FluWeighting_IntT*DensPtr_IntT[k][j][i] )*dv
                                                       :                     DensPtr     [k][j][i]  *dv;
                        break;

                        case WeightByVolume :   Weight = dv;
                        break;

                        default:
                           Aux_Error( ERROR_INFO, "unsupported weight field (%d) !!\n", WeightField );
                           exit( 1 );
                     }


//                   update the profile
                     const int bin = Patch_Bin[TID][LocalID][k][j][i];

                     OMP_Data  [p][TID][bin] += Patch_Data[TID][LocalID][k][j][i]*Weight;
                     OMP_Weight[p][TID][bin] += Weight;
                     OMP_NCell [p][TID][bin] ++;
                  } // i,j,k
               } // for (int LocalID=0; LocalID<8; LocalID++)
            } // for (int p=0; p<NProf; p++)
         } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      } // OpenMP parallel region


//    free particle resources
//    --> these two routines should NOT be put inside an OpenMP parallel region
#     ifdef MASSIVE_PARTICLES
      if ( NeedPar )
      {
         Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

         Prepare_PatchData_FreeParticleDensityArray( lv );
      }
#     endif
   } // for (int lv=MinLv; lv<=MaxLv; lv++)


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

   delete [] Patch_Data;
   delete [] Patch_Bin;


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


// restore the original temporal interpolation setup
   OPT__INT_TIME = IntTimeBackup;

} // FUNCTION : Aux_ComputeProfile
