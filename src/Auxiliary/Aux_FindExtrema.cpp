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
//                3. Only support intrinsic fluid fields and gravitational potential for now
//                   --> Does not support magnetic field, derived fields, and deposited particle fields
//                4. Does not support computing multiple fields at once
//                5. Support periodic BC
//                6. Support finding either the maximum or minimum value
//                   --> Controlled by the parameter "Mode" (see below)
//                7. Support hybrid OpenMP/MPI parallelization
//                   --> All ranks will share the same results after invoking this function
//                8. The parameters Min/MaxLv and PatchType are mainly for performance optimization
//                   --> Simply set MinLv=0, MaxLv=TOP_LEVEL, PatchType=PATCH_LEAF should be OK for most cases
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
//                Extrema.Radius    = HUGE_NUMBER; // entire domain
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

   long SupportedField = _TOTAL;
#  ifdef GRAVITY
   SupportedField |= _POTE;
#  endif
   if ( Extrema->Field == _NONE )
      Aux_Error( ERROR_INFO, "Field == _NONE !!\n" );

   if ( Extrema->Field & ~SupportedField )
      Aux_Error( ERROR_INFO, "unsupported field (%ld) !!\n", Extrema->Field );

   if ( Extrema->Radius <= 0.0 )
      Aux_Error( ERROR_INFO, "Radius (%14.7e) <= 0.0 !!\n", Extrema->Radius );

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

   for (int v=0; v<NCOMP_TOTAL; v++) {
      if ( Extrema->Field & BIDX(v) ) {
         TFluIntIdx = v;
         break;
      }
   }


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


#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

//    initialize the extrema
      OMP_Extrema[TID].Value = ( Mode == EXTREMA_MIN ) ? HUGE_NUMBER : -HUGE_NUMBER;

//    loop over all target levels
      for (int lv=MinLv; lv<=MaxLv; lv++)
      {
         const double dh = amr->dh[lv];

#        pragma omp for schedule( runtime )
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


//          skip distant patches
            const double *EdgeL = amr->patch[0][lv][PID]->EdgeL;
            const double *EdgeR = amr->patch[0][lv][PID]->EdgeR;

            if (   (  ( EdgeL[0]>RegMax[0] || EdgeR[0]<RegMin[0] )  &&  ( EdgeL[0]>RegMax_Img[0] || EdgeR[0]<RegMin_Img[0] )  )   ||
                   (  ( EdgeL[1]>RegMax[1] || EdgeR[1]<RegMin[1] )  &&  ( EdgeL[1]>RegMax_Img[1] || EdgeR[1]<RegMin_Img[1] )  )   ||
                   (  ( EdgeL[2]>RegMax[2] || EdgeR[2]<RegMin[2] )  &&  ( EdgeL[2]>RegMax_Img[2] || EdgeR[2]<RegMin_Img[2] )  )    )
               continue;


//          loop over all cells
            const real (*FluidPtr)[PS1][PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid;
#           ifdef GRAVITY
            const real (*PotPtr  )[PS1][PS1]      = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;
#           endif
            const double x0  = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
            const double y0  = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
            const double z0  = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

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

//             only include cells within the target sphere
               if ( r2 < MaxR2 )
               {
                  real Value;
                  if ( TFluIntIdx != FluIdxUndef )
                     Value = FluidPtr[TFluIntIdx][k][j][i];
#                 ifdef GRAVITY
                  else if ( Extrema->Field & _POTE )
                     Value = PotPtr[k][j][i];
#                 endif
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
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // for (int lv=MinLv; lv<=MaxLv; lv++)
   } // OpenMP parallel region


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
