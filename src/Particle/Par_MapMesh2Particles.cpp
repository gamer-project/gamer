#include "GAMER.h"

#ifdef PARTICLE

void Par_MapMesh2Particles ( const int lv, const int P, const double EdgeL[3],
                             const double EdgeR[3], const int AttrSize3D, const real *Attr,
                             const int NPar, real *InterpParPos[3],
                             const real ParType[], const long ParList[],
                             bool useTracers, real ParAttr[], const int ParGhost )
{

   typedef real (*vla)[AttrSize3D][AttrSize3D][AttrSize3D];
   vla Attr3D = ( vla )Attr;

   const ParInterp_t IntScheme    = amr->Par->Interp;

   const double _dh               = 1.0/amr->dh[lv];

   for (int p=0; p<NPar; p++)
   {
      long ParID = ParList[p];

      if ( useTracers )
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
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idx %d) !!\n",
                             d, InterpParPos[d][p], EdgeL[d], idx[d] );
#              endif

               idx[d] = 0;
            }
            else if ( idx[d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idx %d) !!\n",
                             d, InterpParPos[d][p], EdgeR[d], idx[d] );
#              endif

               idx[d] = AttrSize3D - 1;
            }
         } // for (int d=0; d<3; d++)

//       calculate new particle attribute
         ParAttr[p] = Attr3D[ P ][ idx[2] ][ idx[1] ][ idx[0] ];

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
            dr      [d] = ( InterpParPos[d][p] - EdgeL[d] )*_dh + ParGhost - 0.5;
            idxLR[0][d] = int( dr[d] );
            idxLR[1][d] = idxLR[0][d] + 1;

//          prevent from round-off errors
//          (CIC should be clear of this issue unless round-off errors are comparable to dh)
            if ( idxLR[0][d] < 0 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], EdgeL[d], idxLR[0][d], idxLR[1][d] );
#              endif

               idxLR[0][d] = 0;
               idxLR[1][d] = 1;
            }
            else if ( idxLR[1][d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], EdgeR[d], idxLR[0][d], idxLR[1][d] );
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
         ParAttr[p] = (real)0.0;

         for (int k=0; k<2; k++)
         for (int j=0; j<2; j++)
         for (int i=0; i<2; i++) {
            ParAttr[p] += Attr3D[ P ][ idxLR[k][2] ][ idxLR[j][1] ][ idxLR[i][0] ]
               *Frac[i][0]*Frac[j][1]*Frac[k][2];
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
            dr       [d] = ( InterpParPos[d][p] - EdgeL[d] )*_dh + ParGhost;
            idxLCR[1][d] = int( dr[d] );
            idxLCR[0][d] = idxLCR[1][d] - 1;
            idxLCR[2][d] = idxLCR[1][d] + 1;

//          prevent from round-off errors (especially for NGP and TSC)
            if ( idxLCR[0][d] < 0 )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeL[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeL %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], EdgeL[d], idxLCR[0][d], idxLCR[2][d] );
#              endif

               idxLCR[0][d] = 0;
               idxLCR[1][d] = 1;
               idxLCR[2][d] = 2;
            }
            else if ( idxLCR[2][d] >= AttrSize3D )
            {
#              ifdef DEBUG_PARTICLE
               if (  ! Mis_CompareRealValue( InterpParPos[d][p], (real)EdgeR[d], NULL, false )  )
                  Aux_Error( ERROR_INFO, "index outside the attr array (pos[%d] %14.7e, EdgeR %14.7e, idxL %d, idxR %d) !!\n",
                             d, InterpParPos[d][p], EdgeR[d], idxLCR[0][d], idxLCR[2][d] );
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
         ParAttr[p] = (real)0.0;

         for (int k=0; k<3; k++)
         for (int j=0; j<3; j++)
         for (int i=0; i<3; i++) {

            ParAttr[p] += Attr3D[ P ][ idxLCR[k][2] ][ idxLCR[j][1] ][ idxLCR[i][0] ]
               *Frac[i][0]*Frac[j][1]*Frac[k][2];

         }

      } // PAR_INTERP_TSC
      break;

      default: Aux_Error( ERROR_INFO, "unsupported particle interpolation scheme !!\n" );
      } // switch ( IntScheme )

   } //

} // FUNCTION : Par_MapMesh2Particles

#endif // #ifdef PARTICLE