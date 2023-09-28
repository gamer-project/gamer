#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_FindWeightedCenter
// Description :  Find the location of the weighted-by-field center within a target region
//
// Note        :  The weighted center is defined as
//                  ( \int w(x,y,z)*r(x,y,z) dxdydz ) / ( \int w(x,y,z) dxdydz ),
//                where w is the weighting density field and r is a vector of the position coordinate.
//                If weighting density field is the mass density, the weighted center is the center of mass.
//
// Parameter   :  WeightedCenter         : the coordinate of the weighted center to be returned
//                Center_ref[]           : the coordinate of center of reference
//                MaxR                   : Maximum radius to specify the region to compute the weighted center
//                MinRho                 : Minimum density to specify the region to compute the weighted center
//                Mode                   : How to select the target region
//                                         (0: all region, 1: by MaxR, 2: by MinRho)
//                WeightingDensityField  : the weighting density field used for compuation as the w(x,y,z) in the above Note
//
// Return      :  WeightedCenter[]
//-------------------------------------------------------------------------------------------------------
void Aux_FindWeightedCenter( double WeightedCenter[], const double Center_ref[], const double MaxR, const double MinRho, const int Mode, const long WeightingDensityField )
{

   const double MaxR2             = SQR( MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

   int   *PID0List = NULL;
   double W_ThisRank, WR_ThisRank[3], W_AllRank, WR_AllRank[3];
   real (*WeightingDensity)[PS1][PS1][PS1];

   W_ThisRank = 0.0;
   for (int d=0; d<3; d++)    WR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ|_PAR_TYPE, PredictParPos_No, NULL_REAL,
                                    SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );
#     endif

//    get the total density on grids
      WeightingDensity = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], WeightingDensity[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, WeightingDensityField, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif


//    calculate the weighted-average center
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         double x, y, z, dx, dy, dz;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;  dz = z - Center_ref[2];
                                       if ( Periodic ) {
                                          if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                          else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                       }
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;  dy = y - Center_ref[1];
                                       if ( Periodic ) {
                                          if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                          else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                       }
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;  dx = x - Center_ref[0];
                                       if ( Periodic ) {
                                          if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                          else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                       }

            bool isIncluded = false;

            switch( Mode )
            {
               case 0:
                   isIncluded = true;
                   break;
               case 1:
                   const double R2 = SQR(dx) + SQR(dy) + SQR(dz);
                   isIncluded = ( R2 < MaxR2 );
                   break;
               //case 2:
               //    isIncluded = ( Dens > MinRho );
               //    break;
               default:
                   Aux_Error( ERROR_INFO, "Invalid mode !!\n");
            
            }

//          only include cells satisfying conditions
            if ( isIncluded )
            {
               const double dw = WeightingDensity[PID][k][j][i]*dv;

               W_ThisRank     += dw;
               WR_ThisRank[0] += dw*x;
               WR_ThisRank[1] += dw*y;
               WR_ThisRank[2] += dw*z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] WeightingDensity;
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks to calculate the weighted center
// --> note that all ranks will get WeightedCenteweighted center]
   MPI_Allreduce( &W_ThisRank, &W_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( WR_ThisRank, WR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)    WeightedCenter[d] = WR_AllRank[d] / W_AllRank;

// map the new center back to the simulation domain
   if ( Periodic )
   for (int d=0; d<3; d++)
   {
      if      ( WeightedCenter[d] >= amr->BoxSize[d] )  WeightedCenter[d] -= amr->BoxSize[d];
      else if ( WeightedCenter[d] < 0.0              )  WeightedCenter[d] += amr->BoxSize[d];

   }

   for (int d=0; d<3; d++)
      if ( WeightedCenter[d] >= amr->BoxSize[d]  ||  WeightedCenter[d] < 0.0 )
         Aux_Error( ERROR_INFO, "WeightedCenter[%d] = %14.7e lies outside the domain !!\n", d, WeightedCenter[d] );

} // FUNCTION : Aux_FindWeightedCenter
