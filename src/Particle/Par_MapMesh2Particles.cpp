#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_MapMesh2Particles
// Description :  Map quantities from mesh onto the particles at their positions
//
// Note        :  1. Input grid of mesh data "Attr" is size of a patch + 2*ParGhost
//                   --> ParGhost may be different for tracer vs. active particles
//                2. EdgeL and EdgeR correspond to edges of Attr, not the patch itself
//                3. CorrectVelocity may be used only for mapping velocity, used to
//                   correct tracer particle trajectories in discontinuous flows
//                   --> See Section 2.2 and Equation 1 of Wittor et al. (2016) MNRAS, 464, 4
//                       (https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.4448W)
//                4. Currently only used for tracer particles
//
// Parameter   :  EdgeL           : Left edge of input grid in 3 dimensions
//                EdgeR           : Right edge of input grid in 3 dimensions
//                _dh             : Inverse of cell size
//                AttrSize3D      : Number of cells on side of input grid
//                Attr            : The input grid of values for the variable to be mapped
//                NPar            : The number of particles belonging to this patch
//                InterpParPos    : The positions of the NPar particles on this patch
//                ParType         : The types of the NPar particles on this patch
//                ParList         : The list of particle IDs on this patch
//                UseTracers      : Whether to map to only tracer particles or only active particles
//                ParAttr         : The array to store the mapped particle attribute
//                CorrectVelocity : If true, particle velocities will be corrected in regions of
//                                  discontinuous flow
//
// Return      :  ParAttr[]
//-------------------------------------------------------------------------------------------------------
void Par_MapMesh2Particles( const double EdgeL[3], const double EdgeR[3],
                            const double _dh, const int AttrSize3D, const real *Attr,
                            const int NPar, real_par *InterpParPos[3],
                            const long_par ParType[], const long ParList[],
                            const bool UseTracers, real_par ParAttr[], const bool CorrectVelocity )
{

   typedef real (*vla)[AttrSize3D][AttrSize3D];
   vla Attr3D = ( vla )Attr;

   const ParInterp_t IntScheme = amr->Par->InterpTracer;

   for (int p=0; p<NPar; p++)
   {
      long ParID = ParList[p];

      if ( UseTracers )
      {
//       skip massive particles
         if ( ParType[ParID] != PTYPE_TRACER )
            continue;
      }
      else
      {
//       skip tracer particles
         if ( ParType[ParID] == PTYPE_TRACER )
            continue;
      }

      switch ( IntScheme ) {

//    1 NGP
      case ( PAR_INTERP_NGP ):
      {
         int idx[3];

//       calculate the nearest grid index
         for (int d=0; d<3; d++)
         {
            idx[d] = int( ( InterpParPos[d][p] - EdgeL[d] )*_dh );

//          prevent from round-off errors (especially for NGP and TSC)
            if ( idx[d] < 0 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeL[d], idx[d] );
#              endif

               idx[d] = 0;
            }
            else if ( idx[d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeR[d], idx[d] );
#              endif

               idx[d] = AttrSize3D - 1;
            }
         } // for (int d=0; d<3; d++)

//       calculate new particle attribute
         ParAttr[p] = (real_par)Attr3D[ idx[2] ][ idx[1] ][ idx[0] ];

      } // PAR_INTERP_NGP
      break;

//    2 CIC
      case ( PAR_INTERP_CIC ):
      {
         int    idxLR[2][3];     // array index of the left (idxLR[0][d]) and right (idxLR[1][d]) cells
         double dr      [3];     // distance to the center of the left cell
         double Frac [2][3];     // weighting of the left (Frac[0][d]) and right (Frac[1][d]) cells

         for (int d=0; d<3; d++)
         {
//          calculate the array index of the left and right cells
            dr      [d] = (double)( InterpParPos[d][p] - (real_par)EdgeL[d] )*_dh - 0.5;
            idxLR[0][d] = int( dr[d] );
            idxLR[1][d] = idxLR[0][d] + 1;

//          prevent from round-off errors
//          (CIC should be clear of this issue unless round-off errors are comparable to dh)
            if ( idxLR[0][d] < 0 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeL[d], idxLR[0][d], idxLR[1][d] );
#              endif

               idxLR[0][d] = 0;
               idxLR[1][d] = 1;
            }
            else if ( idxLR[1][d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeR[d], idxLR[0][d], idxLR[1][d] );
#              endif

               idxLR[0][d] = AttrSize3D - 2;
               idxLR[1][d] = AttrSize3D - 1;
            }

//          get the weighting of the nearby 8 cells
            dr     [d] -= (double)idxLR[0][d];
            Frac[0][d]  = 1.0 - dr[d];
            Frac[1][d]  =       dr[d];
         } // for (int d=0; d<3; d++)

//       calculate attribute
         ParAttr[p] = (real_par)0.0;

         for (int k=0; k<2; k++)
         for (int j=0; j<2; j++)
         for (int i=0; i<2; i++) {
            ParAttr[p] += (real_par)(Attr3D[ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
               *Frac[i][0]*Frac[j][1]*Frac[k][2]);
         }

      } // PAR_INTERP_CIC
      break;

//    3 TSC
      case ( PAR_INTERP_TSC ):
      {
         int    idxLCR[3][3];    // array index of the left/central/right cells (idxLCR[0/1/2][d])
         double dr       [3];    // distance to the left edge of the central cell
         double Frac  [3][3];    // weighting of the left/central/right cells (Frac[0/1/2][d])

         for (int d=0; d<3; d++)
         {
//          calculate the array index of the left, central, and right cells
            dr       [d] = (double)( InterpParPos[d][p] - (real_par)EdgeL[d] )*_dh;
            idxLCR[1][d] = int( dr[d] );
            idxLCR[0][d] = idxLCR[1][d] - 1;
            idxLCR[2][d] = idxLCR[1][d] + 1;

//          prevent from round-off errors (especially for NGP and TSC)
            if ( idxLCR[0][d] < 0 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#              endif

               idxLCR[0][d] = 0;
               idxLCR[1][d] = 1;
               idxLCR[2][d] = 2;
            }
            else if ( idxLCR[2][d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
#              endif

               idxLCR[0][d] = AttrSize3D - 3;
               idxLCR[1][d] = AttrSize3D - 2;
               idxLCR[2][d] = AttrSize3D - 1;
            }

//          get the weighting of the nearby 27 cells
            dr     [d] -= (double)idxLCR[1][d];
            Frac[0][d]  = 0.5*SQR( 1.0 - dr[d] );
            Frac[1][d]  = 0.5*( 1.0 + 2.0*dr[d] - 2.0*SQR(dr[d]) );
            Frac[2][d]  = 0.5*SQR( dr[d] );
         } // for (int d=0; d<3; d++)

//       calculate attribute
         ParAttr[p] = (real_par)0.0;

         for (int k=0; k<3; k++)
         for (int j=0; j<3; j++)
         for (int i=0; i<3; i++) {
            ParAttr[p] += (real_par)(Attr3D[ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
               *Frac[i][0]*Frac[j][1]*Frac[k][2]);
         }

      } // PAR_INTERP_TSC
      break;

      default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      } // switch ( IntScheme )

//    If this is a velocity field and we have asked to correct for the mean velocity,
//    we do it here
      if ( CorrectVelocity ) {

         int idx[3];

//       calculate the nearest grid index
         for (int d=0; d<3; d++)
         {
            idx[d] = int( (double)( InterpParPos[d][p] - (real_par)EdgeL[d] )*_dh );

//          prevent from round-off errors (especially for NGP and TSC)
            if ( idx[d] < 1 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idx-1 %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeL[d], idx[d]-1 );
#              endif

               idx[d] = 1;
            }
            else if ( idx[d] >= AttrSize3D - 1 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real_par)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idx+1 %d) !!\n",
                             d, InterpParPos[d][p], (real_par)EdgeR[d], idx[d]+1 );
#              endif

               idx[d] = AttrSize3D - 2;
            }
         } // for (int d=0; d<3; d++)

//       Now that we have the cell containing the particle, we compute the difference
//       of the velocity of that cell from the mean of the 27 surrounding cells
//       (including itself) and add this to the particle velocity
         double deltav = 0.0;

         for (int k=-1; k<=1; k++)
         for (int j=-1; j<=1; j++)
         for (int i=-1; i<=1; i++) {
            deltav += (double)Attr3D[ idx[2]+k ][ idx[1]+j ][ idx[0]+i ];
         }

         deltav = (double)Attr3D[ idx[2] ][ idx[1] ][ idx[0] ] - deltav/27.0;

         ParAttr[p] += (real_par)deltav;
      } // if ( CorrectVelocity )

   } // for (int p=0; p<NPar; p++)

} // FUNCTION : Par_MapMesh2Particles



#endif // #ifdef PARTICLE
