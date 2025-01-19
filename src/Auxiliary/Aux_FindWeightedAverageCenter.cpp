#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_FindWeightedAverageCenter
// Description :  Find the location of the weighted-by-field average center within a target region
//
// Note        :  The weighted average center is defined as
//                  ( \int w(x,y,z)*r(x,y,z) dxdydz ) / ( \int w(x,y,z) dxdydz ),
//                where w is the weighting density field and r is a vector of the position coordinate.
//                If the weighting density field is the mass density, the weighted average center is the center of mass.
//
// Parameter   :  WeightedAverageCenter  : Coordinate of the weighted average center to be returned
//                Center_ref             : Center coordinates for reference (i.e., an initial guess of the center coordinates)
//                MaxR                   : Maximum radius to specify the region to compute the weighted average center
//                MinWD                  : Minimum weighting density to specify the region to compute the weighted average center
//                WeightingDensityField  : Weighting density field w(x,y,z) in the above Note
//                TolErrR                : Maximum tolerated error of distance between the center coordinates during iterations
//                MaxIter                : Maximum number of iterations to find the weighted average center
//                FinaldR                : Record of the distance of the center coordinates update in the last iteration
//                FinalNIter             : Record of the total number of iterations
//
// Example     :  double CoM_Coord[3];   // To find the location of center of mass
//                double CoM_FinaldR;
//                int    CoM_FinalNIter;
//
//                const double CoM_ref[3]    = { amr->BoxCenter[0], amr->BoxCenter[1], amr->BoxCenter[2] };
//                const double CoM_MaxR      = __FLT_MAX__; // entire domain
//                const double CoM_MinRho    = 0.0;
//                const long   CoM_Field     = _TOTAL_DENS;
//                const double CoM_TolErrR   = __FLT_MAX__; // disable it since CoM_MaxR includes the entire domain
//                const int    CoM_MaxIter   = 1;           // no iteration is required since CoM_MaxR includes the entire domain
//
//                Aux_FindWeightedAverageCenter( CoM_Coord, CoM_ref, CoM_MaxR, CoM_MinRho,
//                                               CoM_Field, CoM_TolErrR, CoM_MaxIter, &CoM_FinaldR, &CoM_FinalNIter );
//
//                if ( MPI_Rank == 0 )
//                {
//                   char Filename[MAX_STRING];
//                   sprintf( Filename, "%s", "CoM.txt" );
//                   FILE *File = fopen( Filename, "w" );
//                   fprintf( File, "#%13s%14s%3s%14s %15s %14s %14s %14s %14s %14s %14s %14s %15s %14s %14s %14s %14s\n",
//                                    "Time", "Step", "", "dt",
//                                    "CoM_Field", "CoM_MaxR", "CoM_MinRho", "CoM_TolErrR", "CoM_MaxIter",
//                                    "CoM_ref[x]", "CoM_ref[y]", "CoM_ref[z]",
//                                    "CoM_FinalNIter", "CoM_FinaldR", "CoM_Coord[x]", "CoM_Coord[y]", "CoM_Coord[z]" );
//                   fprintf( File, "%14.7e%14ld%3s%14.7e %15ld %14.7e %14.7e %14.7e %14d %14.7e %14.7e %14.7e %15d %14.7e %14.7e %14.7e %14.7e\n",
//                                    Time[0], Step, "", dTime_Base,
//                                    CoM_Field, CoM_MaxR, CoM_MinRho, CoM_TolErrR, CoM_MaxIter,
//                                    CoM_ref[0], CoM_ref[1], CoM_ref[2],
//                                    CoM_FinalNIter, CoM_FinaldR, CoM_Coord[0], CoM_Coord[1], CoM_Coord[2] );
//                   fclose( File );
//                }
//
// Return      :  WeightedAverageCenter[], FinaldR, FinalNIter
//-------------------------------------------------------------------------------------------------------
void Aux_FindWeightedAverageCenter( double WeightedAverageCenter[], const double Center_ref[], const double MaxR, const double MinWD,
                                    const long WeightingDensityField, const double TolErrR, const int MaxIter, double *FinaldR, int *FinalNIter )
{

// check
#  ifdef GAMER_DEBUG
   if ( WeightingDensityField == _NONE )
      Aux_Error( ERROR_INFO, "WeightingDensityField == _NONE !!\n" );

   if ( WeightingDensityField & (WeightingDensityField-1) )
      Aux_Error( ERROR_INFO, "not support computing multiple fields at once (WeightingDensityField = %ld) !!\n", WeightingDensityField );

   if ( MaxR <= 0.0 )
      Aux_Error( ERROR_INFO, "MaxR (%14.7e) <= 0.0 !!\n", MaxR );

   if ( !Aux_IsFinite( SQR(MaxR) ) )
      Aux_Error( ERROR_INFO, "SQR(MaxR) (%14.7e) overflow !!\n", SQR(MaxR) );

   if ( !Aux_IsFinite( SQR(TolErrR) ) )
      Aux_Error( ERROR_INFO, "SQR(TolErrR) (%14.7e) overflow !!\n", SQR(TolErrR) );

   for (int d=0; d<3; d++) {
      if ( Center_ref[d] < amr->BoxEdgeL[d]  ||  Center_ref[d] > amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "Center_ref[%d] (%14.7e) lies outside the simulation box (BoxEdgeL = %14.7e, BoxEdgeR = %14.7e) !!\n",
                    d, Center_ref[d], amr->BoxEdgeL[d], amr->BoxEdgeR[d] );
   }
#  endif // #ifdef GAMER_DEBUG


// get the integer index of the target intrinsic fluid field
   const int FluIdxUndef = -1;
   int TFluIntIdx = FluIdxUndef;
   bool UsePrepare = true;

   for (int v=0; v<NCOMP_TOTAL; v++) {
      if ( WeightingDensityField & BIDX(v) ) {
         TFluIntIdx = v;
         UsePrepare = false;
         break;
      }
   }

#  ifdef GRAVITY
   if ( WeightingDensityField & _POTE ) UsePrepare = false;
#  endif


   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic[3]       = { OPT__BC_FLU[0] == BC_FLU_PERIODIC,
                                      OPT__BC_FLU[2] == BC_FLU_PERIODIC,
                                      OPT__BC_FLU[4] == BC_FLU_PERIODIC };

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

// initialize the referenced center in the first iteration as the input Center_ref
   const double MaxR2             = SQR( MaxR );
   const double TolErrR2          = SQR( TolErrR );
   double dR2, Center_ref_Iter[3];
   for (int d=0; d<3; d++) Center_ref_Iter[d] = Center_ref[d];
   int NIter = 0;

// if the target region is the entire domain, then there is no need to iterate
   const bool   isWholeBox        = ( MaxR2 >= (SQR(amr->BoxSize[0]) + SQR(amr->BoxSize[1]) + SQR(amr->BoxSize[2])) );

// allocate memory for Prepare_PatchData
   real (*WeightingDensity)[PS1][PS1][PS1] = NULL;
   if ( UsePrepare )
   {
      WeightingDensity = new real [8*NPG_Max][PS1][PS1][PS1]; // 8: number of local patches
   }


// start the iteration to find the center until convergence
   while ( true )
   {
//    set the target region that takes into account periodic BC for excluding distant patches
      double Center_Img[3], RegMax[3], RegMin[3], RegMax_Img[3], RegMin_Img[3];

      for (int d=0; d<3; d++)
      {
         if ( Periodic[d] )
         Center_Img[d] = ( Center_ref_Iter[d] > amr->BoxCenter[d] ) ? Center_ref_Iter[d]-amr->BoxSize[d] : Center_ref_Iter[d]+amr->BoxSize[d];
         else
         Center_Img[d] = Center_ref_Iter[d];

         RegMax    [d] = Center_ref_Iter[d] + MaxR;
         RegMin    [d] = Center_ref_Iter[d] - MaxR;
         RegMax_Img[d] = Center_Img     [d] + MaxR;
         RegMin_Img[d] = Center_Img     [d] - MaxR;
      }

      double W_ThisRank, WR_ThisRank[3], W_AllRank, WR_AllRank[3];
      double WX_ThisRank, WY_ThisRank, WZ_ThisRank;

      W_ThisRank  = 0.0;
      WX_ThisRank = 0.0;
      WY_ThisRank = 0.0;
      WZ_ThisRank = 0.0;


//    loop over all levels
      for (int lv=0; lv<NLEVEL; lv++)
      {
         const int NTotal = amr->NPatchComma[lv][1] / 8;
         int   *PID0_List = new int [NTotal];
         for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;

//       initialize the particle density array (rho_ext) and collect particles to the target level
#        ifdef MASSIVE_PARTICLES
         if ( WeightingDensityField & _PAR_DENS  ||  WeightingDensityField & _TOTAL_DENS )
         {
            Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictPos, Time[lv],
                                          SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );

            Prepare_PatchData_InitParticleDensityArray( lv, Time[lv] );
         }
#        endif

//       calculate the weighted average center
         const double dh = amr->dh[lv];
         const double dv = CUBE( dh );

         for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)
         {
            const int NPG = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp;

//          get the weighting density on grids
            if ( UsePrepare )
            {
               Prepare_PatchData( lv, Time[lv], WeightingDensity[0][0][0], NULL, 0, NPG, PID0_List+Disp, WeightingDensityField, _NONE,
                                  INT_NONE, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                                  MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
            }

//          use the "static" schedule for reproducibility
#           pragma omp parallel for reduction ( +:W_ThisRank, WX_ThisRank, WY_ThisRank, WZ_ThisRank ) schedule( static )
            for (int t=0; t<8*NPG; t++)
            {
               const int PID = 8*Disp + t;

//             skip non-leaf patches
               if ( amr->patch[0][lv][PID]->son != -1 )  continue;

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

               for (int k=0; k<PS1; k++)  {  double z = z0 + k*dh;
                                             double dz = z - Center_ref_Iter[2];
                                             if ( Periodic[2] ) {
                                                if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                                else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                             }
               for (int j=0; j<PS1; j++)  {  double y = y0 + j*dh;
                                             double dy = y - Center_ref_Iter[1];
                                             if ( Periodic[1] ) {
                                                if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                                else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                             }
               for (int i=0; i<PS1; i++)  {  double x = x0 + i*dh;
                                             double dx = x - Center_ref_Iter[0];
                                             if ( Periodic[0] ) {
                                                if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                                else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                             }

                  const double R2 = SQR(dx) + SQR(dy) + SQR(dz);

//                only include cells that are within a sphere with radius MaxR and with the weighting density larger than MinWD
                  if ( R2 < MaxR2 )
                  {
                     double WD;
                     if ( TFluIntIdx != FluIdxUndef )
                        WD = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[TFluIntIdx][k][j][i];
                     else if ( UsePrepare )
                        WD =  WeightingDensity[t][k][j][i];
#                    ifdef GRAVITY
                     else if ( WeightingDensityField & _POTE )
                        WD = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];
#                    endif
                     else
                        Aux_Error( ERROR_INFO, "unsupported field (%ld) !!\n", WeightingDensityField );

                     if ( WD > MinWD )
                     {
                        const double dw = WD*dv; // weighting

                        W_ThisRank  += dw;
                        WX_ThisRank += dw*x;
                        WY_ThisRank += dw*y;
                        WZ_ThisRank += dw*z;
                     }
                  }
               }}} // i,j,k
            } // for (int t=0; t<8*NPG; t++)
         } // for (int Disp=0; Disp<NTotal; Disp+=NPG_Max)

//       free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#        ifdef MASSIVE_PARTICLES
         if ( WeightingDensityField & _PAR_DENS  ||  WeightingDensityField & _TOTAL_DENS )
         {
            Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

            Prepare_PatchData_FreeParticleDensityArray( lv );
         }
#        endif

         delete [] PID0_List;

      } // for (int lv=0; lv<NLEVEL; lv++)

      WR_ThisRank[0] = WX_ThisRank;
      WR_ThisRank[1] = WY_ThisRank;
      WR_ThisRank[2] = WZ_ThisRank;


//    collect data from all ranks to calculate the weighted average center
//    --> note that all ranks will get WeightedAverageCenter[]
      MPI_Allreduce( &W_ThisRank, &W_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( WR_ThisRank, WR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      if ( W_AllRank == 0.0 )
      {
         if ( MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : Weighted average center cannot be found because the total weighting (W_AllRank) = %14.7e !!\n", W_AllRank );

//       return a huge number to indicate this routine fails to find the center
         for (int d=0; d<3; d++) WeightedAverageCenter[d] = HUGE_NUMBER;
         dR2 = HUGE_NUMBER;
         NIter++;

         break;
      }

      for (int d=0; d<3; d++) WeightedAverageCenter[d] = WR_AllRank[d] / W_AllRank;

//    map the new center back to the simulation domain
      for (int d=0; d<3; d++)
      {
         if ( Periodic[d] )
         {
            if      ( WeightedAverageCenter[d] >= amr->BoxEdgeR[d] ) WeightedAverageCenter[d] -= amr->BoxSize[d];
            else if ( WeightedAverageCenter[d] <  amr->BoxEdgeL[d] ) WeightedAverageCenter[d] += amr->BoxSize[d];
         }
      }

      for (int d=0; d<3; d++)
         if ( WeightedAverageCenter[d] >= amr->BoxEdgeR[d]  ||  WeightedAverageCenter[d] < amr->BoxEdgeL[d] )
            Aux_Error( ERROR_INFO, "WeightedAverageCenter[%d] = %14.7e lies outside the domain (BoxEdgeL = %14.7e, BoxEdgeR = %14.7e) !!\n",
                       d, WeightedAverageCenter[d], amr->BoxEdgeL[d], amr->BoxEdgeR[d] );


      dR2 = SQR( Center_ref_Iter[0] - WeightedAverageCenter[0] )
          + SQR( Center_ref_Iter[1] - WeightedAverageCenter[1] )
          + SQR( Center_ref_Iter[2] - WeightedAverageCenter[2] );
      NIter++;

//    check the convergence and number of iterations to decide whether to end the iteration
      if ( dR2 <= TolErrR2  ||  NIter >= MaxIter  ||  isWholeBox )
         break;
      else
//       use the weighted average center as the referenced center in the next iteration
         memcpy( Center_ref_Iter, WeightedAverageCenter, sizeof(double)*3 );

   } // while ( true )


// free memory for Prepare_PatchData
   if ( UsePrepare )
   {
      delete [] WeightingDensity;
   }

   if ( MPI_Rank == 0 )
   {
      if ( ! isWholeBox  &&  dR2 > TolErrR2 )
      {
         Aux_Message( stderr, "WARNING : In %s, target region is not the entire domain\n", __FUNCTION__ );
         Aux_Message( stderr, "          and dR (%13.7e) > TolErrR (%13.7e),\n", sqrt(dR2), TolErrR );
         Aux_Message( stderr, "          the weighted average center may not have converged yet when the iteration stopped !!\n" );
      }
   }

// return the information about the iteration
   *FinaldR    = sqrt(dR2);
   *FinalNIter = NIter;

} // FUNCTION : Aux_FindWeightedAverageCenter
