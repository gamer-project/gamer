#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_FindExtrema
// Description :  Find the value and location of the extreme value of a given field within a target
//                spherical region
//
// Note        :  1. Must set "Field, Radius, Center" in the input "Extrema" object in advance
//                   --> Extrema->Field  : Target field (e.g., _DENS or BIDX(DENS))
//                       Extrema->Radius : Target radius
//                       Extrema->Center : Target center
//                   --> "Extrema_t" structure is defined in "include/Extrema.h"
//                2. Results will be stored in "Extrema"
//                   --> Extrema->Value  : Output extreme value
//                       Extrema->Coord  : Coordinates of the extrema
//                       Extrema->Rank   : MPI rank of the extrema
//                       Extrema->Level  : AMR level of the extrema
//                       Extrema->PID    : Local Patch index of the extrema
//                       Extrema->Cell   : Cell indices within a patch of the extrema
//                3. Support fields that are supported by Prepare_PatchData()
//                4. Does not support computing multiple fields at once
//                5. Support periodic BC
//                6. Support finding either the maximum or minimum value
//                   --> Controlled by the parameter "Mode" (see below)
//                7. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same results after invoking this function
//                8. The parameters Min/MaxLv and PatchType are mainly for performance optimization
//                   --> Simply set MinLv=0, MaxLv=TOP_LEVEL, PatchType=PATCH_LEAF should be OK for most cases
//                9. In case different AMR levels are not synchronized, this function currently only checks
//                   the most recent data on each level (i.e., data associated with FluSg[lv]/PotSg[lv]) without temporal interpolation
//
// Parameter   :  Extrema   : Extrema_t object storing the input and output information of the extrema
//                Mode      : EXTREMA_MIN/MAX --> find the minimum/maximum value
//                Min/MaxLv : Consider patches on levels from MinLv to MaxLv
//                PatchType : Only consider patches of the specified type
//                            --> Supported types: PATCH_LEAF, PATCH_NONLEAF, PATCH_BOTH, PATCH_LEAF_PLUS_MAXNONLEAF
//                            --> PATCH_LEAF_PLUS_MAXNONLEAF includes leaf patches on all target levels
//                                (i.e., MinLv ~ MaxLv) and non-leaf patches only on MaxLv
//
// Example     :  Extrema_t Extrema;
//                Extrema.Field     = _DENS;
//                Extrema.Radius    = __FLT_MAX__; // entire domain
//                Extrema.Center[0] = amr->BoxCenter[0];
//                Extrema.Center[1] = amr->BoxCenter[1];
//                Extrema.Center[2] = amr->BoxCenter[2];
//
//                Aux_FindExtrema( &Extrema, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );
//
//                if ( MPI_Rank == 0 )
//                {
//                   char Filename[MAX_STRING];
//                   sprintf( Filename, "%s", "Extrema.txt" );
//                   FILE *File = fopen( Filename, "w" );
//                   fprintf( File, "#%13s%14s%3s%14s %10s %13s %13s %13s %13s %14s %13s %13s %13s %5s %5s %5s %4s %4s %4s\n",
//                            "Time", "Step", "", "dt", "Field", "Radius", "Center[x]", "Center[y]", "Center[z]",
//                            "Value", "Coord[x]", "Coord[y]", "Coord[z]", "Rank", "Level", "PID", "i", "j", "k" );
//                   fprintf( File, "%14.7e%14ld%3s%14.7e", Time[0], Step, "", dTime_Base );
//                   fprintf( File, " %10ld %13.7e %13.7e %13.7e %13.7e %14.7e %13.7e %13.7e %13.7e %5d %5d %5d %4d %4d %4d\n",
//                            Extrema.Field, Extrema.Radius, Extrema.Center[0], Extrema.Center[1], Extrema.Center[2],
//                            Extrema.Value, Extrema.Coord[0], Extrema.Coord[1], Extrema.Coord[2],
//                            Extrema.Rank, Extrema.Level, Extrema.PID, Extrema.Cell[0], Extrema.Cell[1], Extrema.Cell[2] );
//                   fclose( File );
//                }
//
// Return      :  Extrema
//-------------------------------------------------------------------------------------------------------
void Aux_FindExtrema( Extrema_t *Extrema, const ExtremaMode_t Mode, const int MinLv, const int MaxLv,
                      const PatchType_t PatchType )
{

// check
#  ifdef GAMER_DEBUG
   if ( Extrema == NULL )
      Aux_Error( ERROR_INFO, "Extrema == NULL !!\n" );

   if ( Extrema->Field == _NONE )
      Aux_Error( ERROR_INFO, "Field == _NONE !!\n" );

   if ( Extrema->Field & (Extrema->Field-1) )
      Aux_Error( ERROR_INFO, "not support computing multiple fields at once (Extrema->Field = %ld) !!\n", Extrema->Field );

   if ( Extrema->Radius <= 0.0 )
      Aux_Error( ERROR_INFO, "Radius (%14.7e) <= 0.0 !!\n", Extrema->Radius );

   if ( !Aux_IsFinite( SQR(Extrema->Radius) ) )
      Aux_Error( ERROR_INFO, "SQR(Extrema->Radius) (%14.7e) overflow !!\n", SQR(Extrema->Radius) );

   for (int d=0; d<3; d++) {
      if ( Extrema->Center[d] < amr->BoxEdgeL[d]  ||  Extrema->Center[d] > amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "Center[%d] (%14.7e) lies outside the simulation box !!\n", d, Extrema->Center[d] );
   }

   if ( Mode != EXTREMA_MIN  &&  Mode != EXTREMA_MAX )
      Aux_Error( ERROR_INFO, "incorrect Mode (%d) !!\n", Mode );

   if ( MinLv < 0  ||  MinLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MinLv (%d) !!\n", MinLv );

   if ( MaxLv < 0  ||  MaxLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect MaxLv (%d) !!\n", MaxLv );

   if ( MinLv > MaxLv )
      Aux_Error( ERROR_INFO, "MinLv (%d) > MaxLv (%d) !!\n", MinLv, MaxLv );

   if ( PatchType != PATCH_LEAF  &&  PatchType != PATCH_NONLEAF  &&
        PatchType != PATCH_BOTH  &&  PatchType != PATCH_LEAF_PLUS_MAXNONLEAF )
      Aux_Error( ERROR_INFO, "incorrect PatchType (%d) !!\n", PatchType );
#  endif // #ifdef GAMER_DEBUG


// get the integer index of the target intrinsic fluid field
   const int FluIdxUndef = -1;
   int TFluIntIdx = FluIdxUndef;
   bool UsePrepare = true;

   for (int v=0; v<NCOMP_TOTAL; v++) {
      if ( Extrema->Field & BIDX(v) ) {
         TFluIntIdx = v;
         UsePrepare = false;
         break;
      }
   }

#  ifdef GRAVITY
   if ( Extrema->Field & _POTE ) UsePrepare = false;
#  endif


   const long   Field       = Extrema->Field;
   const double MaxR        = Extrema->Radius;
   const double MaxR2       = SQR( MaxR );
   const double HalfBox[3]  = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3] = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                OPT__BC_FLU[4] == BC_FLU_PERIODIC };
#  ifdef OPENMP
   const int    NT          = OMP_NTHREAD;
#  else
   const int    NT          = 1;
#  endif

   Extrema_t OMP_Extrema[NT];


// set the target region that takes into account periodic BC for excluding distant patches
   double Center[3], Center_Img[3], RegMax[3], RegMin[3], RegMax_Img[3], RegMin_Img[3];

   for (int d=0; d<3; d++)
   {
      Center    [d] = Extrema->Center[d];
      if ( Periodic[d] )
      Center_Img[d] = ( Center[d] > amr->BoxCenter[d] ) ? Center[d]-amr->BoxSize[d] : Center[d]+amr->BoxSize[d];
      else
      Center_Img[d] = Center[d];

      RegMax    [d] = Center    [d] + MaxR;
      RegMin    [d] = Center    [d] - MaxR;
      RegMax_Img[d] = Center_Img[d] + MaxR;
      RegMin_Img[d] = Center_Img[d] - MaxR;
   }

   const int    NPG_Max           = FLU_GPU_NPGROUP;
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef MASSIVE_PARTICLES
   const bool   TimingSendPar_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   PredictPos        = amr->Par->PredictPos;
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   PredictPos        = false;
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef MASSIVE_PARTICLES


// initialize the extrema
   for (int TID=0; TID<NT; TID++) OMP_Extrema[TID].Value = ( Mode == EXTREMA_MIN ) ? HUGE_NUMBER : -HUGE_NUMBER;

// allocate memory for Prepare_PatchData
   real (*FieldPtr)[PS1][PS1][PS1] = NULL;
   if ( UsePrepare )
   {
      FieldPtr = new real [8*NPG_Max][PS1][PS1][PS1]; // 8: number of local patches
   }


// loop over all target levels
   for (int lv=MinLv; lv<=MaxLv; lv++)
   {
      const double dh  = amr->dh[lv];
      const int NTotal = amr->NPatchComma[lv][1] / 8;

      int   *PID0_List = new int [NTotal];
      for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;

//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef MASSIVE_PARTICLES
      if ( Extrema->Field & _PAR_DENS  ||  Extrema->Field & _TOTAL_DENS )
      {
         Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictPos, Time[lv],
                                       SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

         Prepare_PatchData_InitParticleDensityArray( lv, Time[lv] );
      }
#     endif

      for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)
      {
         const int NPG = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp;

         if ( UsePrepare )
         {
            Prepare_PatchData( lv, Time[lv], FieldPtr[0][0][0], NULL, 0, NPG, PID0_List+Disp, Field, _NONE,
                               INT_NONE, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                               MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
         }

#        pragma omp parallel
         {
#           ifdef OPENMP
            const int TID = omp_get_thread_num();
#           else
            const int TID = 0;
#           endif

#           pragma omp for schedule( runtime )
            for (int t=0; t<8*NPG; t++)
            {
               const int PID = 8*Disp + t;

//             skip untargeted patches
               if ( amr->patch[0][lv][PID]->son != -1 )
               {
                  if ( PatchType == PATCH_LEAF )                                    continue;
                  if ( PatchType == PATCH_LEAF_PLUS_MAXNONLEAF  &&  lv != MaxLv )   continue;
               }

               else
               {
                  if ( PatchType == PATCH_NONLEAF )                                 continue;
               }


//             skip distant patches
               const double *EdgeL = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR = amr->patch[0][lv][PID]->EdgeR;

               if (   (  ( EdgeL[0]>RegMax[0] || EdgeR[0]<RegMin[0] )  &&  ( EdgeL[0]>RegMax_Img[0] || EdgeR[0]<RegMin_Img[0] )  )   ||
                      (  ( EdgeL[1]>RegMax[1] || EdgeR[1]<RegMin[1] )  &&  ( EdgeL[1]>RegMax_Img[1] || EdgeR[1]<RegMin_Img[1] )  )   ||
                      (  ( EdgeL[2]>RegMax[2] || EdgeR[2]<RegMin[2] )  &&  ( EdgeL[2]>RegMax_Img[2] || EdgeR[2]<RegMin_Img[2] )  )    )
                  continue;


//             loop over all cells
               const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
               const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
               const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

               for (int k=0; k<PS1; k++)  {  const double z = z0 + k*dh;
                                             double dz = z - Center[2];
                                             if ( Periodic[2] ) {
                                                if      ( dz > +HalfBox[2] )  {  dz -= amr->BoxSize[2];  }
                                                else if ( dz < -HalfBox[2] )  {  dz += amr->BoxSize[2];  }
                                             }
               for (int j=0; j<PS1; j++)  {  const double y = y0 + j*dh;
                                             double dy = y - Center[1];
                                             if ( Periodic[1] ) {
                                                if      ( dy > +HalfBox[1] )  {  dy -= amr->BoxSize[1];  }
                                                else if ( dy < -HalfBox[1] )  {  dy += amr->BoxSize[1];  }
                                             }
               for (int i=0; i<PS1; i++)  {  const double x = x0 + i*dh;
                                             double dx = x - Center[0];
                                             if ( Periodic[0] ) {
                                                if      ( dx > +HalfBox[0] )  {  dx -= amr->BoxSize[0];  }
                                                else if ( dx < -HalfBox[0] )  {  dx += amr->BoxSize[0];  }
                                             }

                  const double r2 = SQR(dx) + SQR(dy) + SQR(dz);

//                only include cells within the target sphere
                  if ( r2 < MaxR2 )
                  {
                     real Value;
                     if ( TFluIntIdx != FluIdxUndef )
                        Value = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[TFluIntIdx][k][j][i];
                     else if ( UsePrepare )
                        Value = FieldPtr[t][k][j][i];
#                    ifdef GRAVITY
                     else if ( Extrema->Field & _POTE )
                        Value = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];
#                    endif
                     else
                        Aux_Error( ERROR_INFO, "unsupported field (%ld) !!\n", Extrema->Field );

                     if (  ( Mode == EXTREMA_MAX && Value > OMP_Extrema[TID].Value )  ||
                           ( Mode == EXTREMA_MIN && Value < OMP_Extrema[TID].Value )   )
                     {
                        OMP_Extrema[TID].Value    = Value;
                        OMP_Extrema[TID].Coord[0] = x;
                        OMP_Extrema[TID].Coord[1] = y;
                        OMP_Extrema[TID].Coord[2] = z;
                        OMP_Extrema[TID].Rank     = MPI_Rank;
                        OMP_Extrema[TID].Level    = lv;
                        OMP_Extrema[TID].PID      = PID;
                        OMP_Extrema[TID].Cell[0]  = i;
                        OMP_Extrema[TID].Cell[1]  = j;
                        OMP_Extrema[TID].Cell[2]  = k;
                     }
                  } // if ( r2 < MaxR2 )
               }}} // i,j,k
            } // for (int t=0; t<8*NPG; t++)
         } // OpenMP parallel region
      } // for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)

//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef MASSIVE_PARTICLES
      if ( Extrema->Field & _PAR_DENS  ||  Extrema->Field & _TOTAL_DENS )
      {
         Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

         Prepare_PatchData_FreeParticleDensityArray( lv );
      }
#     endif

      delete [] PID0_List;

   } // for (int lv=MinLv; lv<=MaxLv; lv++)


// free memory for Prepare_PatchData
   if ( UsePrepare )
   {
      delete [] FieldPtr;
   }

// find the extrema over all OpenMP threads
   Extrema->Value = ( Mode == EXTREMA_MIN ) ? HUGE_NUMBER : -HUGE_NUMBER;

   for (int TID=0; TID<NT; TID++)
   {
      if (  ( Mode == EXTREMA_MAX && OMP_Extrema[TID].Value > Extrema->Value )  ||
            ( Mode == EXTREMA_MIN && OMP_Extrema[TID].Value < Extrema->Value )   )
         *Extrema = OMP_Extrema[TID];
   }

// restore the input fields
   Extrema->Field     = Field;
   Extrema->Radius    = MaxR;
   Extrema->Center[0] = Center[0];
   Extrema->Center[1] = Center[1];
   Extrema->Center[2] = Center[2];


// for MPI only
#  ifndef SERIAL
// define an MPI derived datatype for Extrema_t
   MPI_Datatype MPI_Extrema_t;

   Extrema->CreateMPIType( &MPI_Extrema_t );


// get the extrema among all ranks
   const int RootRank = 0;
   Extrema_t Extrema_AllRank[MPI_NRank];

   MPI_Gather( Extrema, 1, MPI_Extrema_t, Extrema_AllRank, 1, MPI_Extrema_t, RootRank, MPI_COMM_WORLD );

   if ( MPI_Rank == RootRank )
   for (int r=0; r<MPI_NRank; r++)
   {
      if (  ( Mode == EXTREMA_MAX && Extrema_AllRank[r].Value > Extrema->Value )  ||
            ( Mode == EXTREMA_MIN && Extrema_AllRank[r].Value < Extrema->Value )   )
         *Extrema = Extrema_AllRank[r];
   }


// broadcast the extrema to all ranks
   MPI_Bcast( Extrema, 1, MPI_Extrema_t, RootRank, MPI_COMM_WORLD );


// free the MPI derived datatype
   MPI_Type_free( &MPI_Extrema_t );
#  endif // #ifndef SERIAL

} // FUNCTION : Aux_FindExtrema
