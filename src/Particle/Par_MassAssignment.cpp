#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

static bool WithinRho( const int idx[], const int RhoSize );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_MassAssignment
// Description :  Deposit particle mass to grid
//
// Note        :  1. Three different schemes are supported
//                   (1) NGP : Nearest-Grid-Point
//                   (2) CIC : Could-In-Cell
//                   (3) TSC : Triangular-Shaped-Cloud
//                2. The deposited density field will be stored in Rho
//                   --> This array will be initialized as zero only if "InitZero=true"
//                3. Particles having no contribution to Rho (which has the range EdgeL[d] <= r[d] < EdgeL[d]+RhoSize*dh )
//                   will be ignored
//                4. Particles position will be predicted to the target physical time if PredictPos is on
//                   --> But they will NOT be stored back to the global Pos array
//                   --> Also remember to skip particles waiting for velocity correction since they have time
//                       temporarily set to -dt (moreover, they should already be synchronized with TargetTime)
//
// Parameter   :  ParList      : List of target particle IDs
//                NPar         : Number of particles
//                IntScheme    : Particle interpolation scheme
//                Rho          : Array to store the output density field (assumed to be a cubic array)
//                RhoSize      : Size of Rho in 1D
//                EdgeL        : Left edge of the array Rho
//                dh           : cell size of Rho
//                PredictPos   : true --> predict particle position to TargetTime
//                TargetTime   : Target time for predicting the particle position
//                InitZero     : True --> initialize Rho as zero
//                Periodic     : True --> apply periodic boundary condition
//                PeriodicSize : Number of cells in the periodic box (in the unit of dh)
//                UnitDens     : Assign unit density to each particle regardless of the real particle mass
//                               and cell size
//                               --> Useful for counting the number of particles on each cell (together
//                                   with IntScheme == PAR_INTERP_NGP), which can be used as a refinement
//                                   criterion (OPT__FLAG_NPAR_CELL)
//              
// Return      :  Rho
//-------------------------------------------------------------------------------------------------------
void Par_MassAssignment( const long *ParList, const long NPar, const ParInterp_t IntScheme, real *Rho,
                         const int RhoSize, const double *EdgeL, const double dh, const bool PredictPos,
                         const double TargetTime, const bool InitZero, const bool Periodic, const int PeriodicSize[3],
                         const bool UnitDens )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( NPar > 0  &&  ParList == NULL )      Aux_Error( ERROR_INFO, "ParList == NULL for NPar = %ld !!\n", NPar );
   if ( Rho == NULL )                        Aux_Error( ERROR_INFO, "Rho == NULL !!\n" );
   if ( EdgeL == NULL )                      Aux_Error( ERROR_INFO, "EdgeL == NULL !!\n" );
   if ( PredictPos  &&  TargetTime < 0.0 )   Aux_Error( ERROR_INFO, "TargetTime = %14.7e < 0.0 !!\n", TargetTime );
   if ( Periodic )
      for (int d=0; d<3; d++)
         if ( RhoSize > PeriodicSize[d] )
            Aux_Error( ERROR_INFO, "RhoSize (%d) > PeriodicSize[%d] (%d) !!\n", RhoSize, d, PeriodicSize[d] );
#  endif


// 1. initialization
   if ( InitZero )
      for (int t=0; t<CUBE(RhoSize); t++)    Rho[t] = (real)0.0;

   if ( NPar == 0 )  return;


// 2. copy particle position since they might be modified during the position prediction
   real *Pos[3] = { NULL, NULL, NULL };
   long  ParID;

   for (int d=0; d<3; d++)    Pos[d] = new real [NPar];

   for (long p=0; p<NPar; p++)
   {
      ParID = ParList[p];

      Pos[0][p] = amr->Par->PosX[ParID];
      Pos[1][p] = amr->Par->PosY[ParID];
      Pos[2][p] = amr->Par->PosZ[ParID];
   }


// 3. predict particle position
   real dt;

   if ( PredictPos )
   {
      for (long p=0; p<NPar; p++)
      {
         ParID = ParList[p];

//       skip particles waiting for velocity correction (they should already be synchronized with TargetTime)
         if ( amr->Par->Time[ParID] < (real)0.0 )  continue;

         dt = (real)TargetTime - amr->Par->Time[ParID];

         Pos[0][p] += amr->Par->VelX[ParID]*dt;
         Pos[1][p] += amr->Par->VelY[ParID]*dt;
         Pos[2][p] += amr->Par->VelZ[ParID]*dt;
      }
   } // if ( PredictPos )


// 4. deposit particle mass
   const double _dh  = 1.0 / dh;
   const double _dh3 = CUBE(_dh);

   real (*Rho3D)[RhoSize][RhoSize] = ( real (*)[RhoSize][RhoSize] )Rho;

   int  idx[3];      // array index for Rho
   real ParDens;     // mass density of the cloud

   switch ( IntScheme )
   {
//    4.1 NGP
      case ( PAR_INTERP_NGP ):
      {
         for (long p=0; p<NPar; p++)
         {
            ParID = ParList[p];

//          4.1.1 calculate the nearest grid index
            for (int d=0; d<3; d++)    
            {
               idx[d] = (int)FLOOR( ( Pos[d][p] - EdgeL[d] )*_dh );

//             periodicity
               if ( Periodic )
               {
                  idx[d] = ( idx[d] + PeriodicSize[d] ) % PeriodicSize[d];

#                 ifdef DEBUG_PARTICLE
                  if ( idx[d] < 0  ||  idx[d] >= PeriodicSize[d] )
                     Aux_Error( ERROR_INFO, "incorrect idx[%d] = %d (PeriodicSize = %d) !!\n",
                                d, idx[d], PeriodicSize[d] );
#                 endif
               }
            }

//          4.1.2 assign mass if within the Rho array
//          check for inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Mass[ParID] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, amr->Par->Mass[ParID] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = amr->Par->Mass[ParID]*_dh3;

            if (  WithinRho( idx, RhoSize )  )
               Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens;
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_NGP
      break;


//    4.2 CIC
      case ( PAR_INTERP_CIC ):
      {
         int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
         double dr      [3];     // distance to the center of the left cell
         double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

         for (long p=0; p<NPar; p++)
         {
            ParID = ParList[p];

            for (int d=0; d<3; d++)
            {
//             4.2.1 calculate the array index of the left and right cells
               dr      [d]  = ( Pos[d][p] - EdgeL[d] )*_dh - 0.5;
               idxLR[0][d]  = (int)FLOOR( dr[d] );
               idxLR[1][d]  = idxLR[0][d] + 1;
               dr      [d] -= (double)idxLR[0][d];

//             periodicity
               if ( Periodic )
               {
                  for (int t=0; t<2; t++)
                  {
                     idxLR[t][d] = ( idxLR[t][d] + PeriodicSize[d] ) % PeriodicSize[d];

#                    ifdef DEBUG_PARTICLE
                     if ( idxLR[t][d] < 0  ||  idxLR[t][d] >= PeriodicSize[d] )
                        Aux_Error( ERROR_INFO, "incorrect idxLR[%d][%d] = %d (PeriodicSize = %d) !!\n",
                                   t, d, idxLR[t][d], PeriodicSize[d] );
#                    endif
                  }
               }

//             4.2.2 get the weighting of the nearby 8 cells
               Frac[0][d] = 1.0 - dr[d];
               Frac[1][d] =       dr[d];
            } // for (int d=0; d<3; d++)

//          4.2.3 assign mass if within the Rho array
//          check for inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Mass[ParID] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, amr->Par->Mass[ParID] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = amr->Par->Mass[ParID]*_dh3;

            for (int k=0; k<2; k++) {  idx[2] = idxLR[k][2];
            for (int j=0; j<2; j++) {  idx[1] = idxLR[j][1];
            for (int i=0; i<2; i++) {  idx[0] = idxLR[i][0];

               if (  WithinRho( idx, RhoSize )  )
                  Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens*Frac[i][0]*Frac[j][1]*Frac[k][2];

            }}}
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_CIC
      break;


//    4.3 TSC
      case ( PAR_INTERP_TSC ):
      {
         int    idxLCR[3][3];    // array index of the left (idxLCR[0][d]), central (idxLCR[1][d]) and right (idxLCR[2][d]) cells
         double dr       [3];    // distance to the left edge of the central cell
         double Frac  [3][3];    // weighting of the left (Frac[0][d]), central (Frac[1][d]) and right (Frac[2][d]) cells

         for (long p=0; p<NPar; p++)
         {
            ParID = ParList[p];

            for (int d=0; d<3; d++)
            {
//             4.3.1 calculate the array index of the left, central, and right cells
               dr       [d]  = ( Pos[d][p] - EdgeL[d] )*_dh;
               idxLCR[1][d]  = (int)FLOOR( dr[d] );
               idxLCR[0][d]  = idxLCR[1][d] - 1;
               idxLCR[2][d]  = idxLCR[1][d] + 1;
               dr       [d] -= (double)idxLCR[1][d];

//             periodicity
               if ( Periodic )
               {
                  for (int t=0; t<3; t++)
                  {
                     idxLCR[t][d]  = ( idxLCR[t][d] + PeriodicSize[d] ) % PeriodicSize[d];

#                    ifdef DEBUG_PARTICLE
                     if ( idxLCR[t][d] < 0  ||  idxLCR[t][d] >= PeriodicSize[d] )
                        Aux_Error( ERROR_INFO, "incorrect idxLCR[%d][%d] = %d (PeriodicSize = %d) !!\n",
                                   t, d, idxLCR[t][d], PeriodicSize[d] );
#                    endif
                  }
               }

//             4.3.2 get the weighting of the nearby 27 cells
               Frac[0][d] = 0.5*SQR( 1.0 - dr[d] );
               Frac[1][d] = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
               Frac[2][d] = 0.5*SQR( dr[d] );
            } // for (int d=0; d<3; d++)

//          4.3.3 assign mass if within the Rho array
//          check for inactive particles (which have negative mass)
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Mass[ParID] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, amr->Par->Mass[ParID] );
#           endif

            if ( UnitDens )   ParDens = (real)1.0;
            else              ParDens = amr->Par->Mass[ParID]*_dh3;

            for (int k=0; k<3; k++) {  idx[2] = idxLCR[k][2];
            for (int j=0; j<3; j++) {  idx[1] = idxLCR[j][1];
            for (int i=0; i<3; i++) {  idx[0] = idxLCR[i][0];

               if (  WithinRho( idx, RhoSize )  )
                  Rho3D[ idx[2] ][ idx[1] ][ idx[0] ] += ParDens*Frac[i][0]*Frac[j][1]*Frac[k][2];
            }}}
         } // for (long p=0; p<NPar; p++)
      } // PAR_INTERP_TSC
      break;

      default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
   } // switch ( IntScheme )


// 5. free memory
   for (int d=0; d<3; d++)    delete [] Pos[d];

} // FUNCTION : Par_MassAssignment



//-------------------------------------------------------------------------------------------------------
// Function    :  WithinRho
// Description :  Check whether the input cell lies within the Rho array
//
// Parameter   :  idx      : Cell indices 
//                RhoSize  : Size of the Rho array
//
// Return      :  true --> within the target region
//-------------------------------------------------------------------------------------------------------
bool WithinRho( const int idx[], const int RhoSize )
{

   for (int d=0; d<3; d++)
      if ( idx[d] < 0  ||  idx[d] >= RhoSize )  return false;

   return true;

} // FUNCTION : WithinRho



#endif // #ifdef PARTICLE
