#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Int_MinMod1D
// Description :  Perform spatial interpolation based on the MinMod limiter
//
// Note        :  1. MinMod limiter:
//                   (1) Slope at each cell is determined by the minimum slope between the left and right slopes
//                   (2) Set slope to zero if the left and right **slopes** have opposite signs
//                   (3) When enabling "OppSignOthOrder", set slope to zero if the left and right **values**
//                       have opposite signs
//                       --> Useful to avoid introducing unphysically large velocity in low-density valleys
//                       --> Because when momentum changes sign (e.g., -1e7, 0, +1e7), even though still
//                           monotonic, the interpolated momentum may become much larger (e.g., O(1e7))
//                           than the central value, which corresponds to extremely large velocity due to
//                           the low central density
//                2. BOTH conservative and monotonic
//                3. 3D interpolation is realized by computing the slopes in all three directions at once
//
// Parameter   :  CData           : Input coarse-grid array
//                CSize           : Size of the CData array
//                CStart          : (x,y,z) starting indices to perform interpolation on the CData array
//                CRange          : Number of grids in each direction to perform interpolation
//                FData           : Output fine-grid array
//                FStart          : (x,y,z) starting indcies to store the interpolation results
//                NComp           : Number of components in the CData and FData array
//                OppSign0thOrder : See the note above
//-------------------------------------------------------------------------------------------------------
void Int_MinMod1D( const real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                   real FData[], const int FSize[3], const int FStart[3], const int NComp,
                   const bool OppSign0thOrder )
{

   const int Cdx = 1;
   const int Cdy = CSize[0];
   const int Cdz = CSize[0]*CSize[1];

   const int Fdx = 1;
   const int Fdy = FSize[0];
   const int Fdz = FSize[0]*FSize[1];

   real Slope_x, Slope_y, Slope_z, LSlope, RSlope;
   int Cx, Cy, Cz, Fx, Fy, Fz, CID, FID, CID0, FID0;


   for (int v=0; v<NComp; v++)
   {
      CID0 = v*CSize[0]*CSize[1]*CSize[2];
      FID0 = v*FSize[0]*FSize[1]*FSize[2];

      for ( Cz=CStart[2], Fz=FStart[2]; Cz<CStart[2]+CRange[2]; Cz++, Fz+=2 )
      for ( Cy=CStart[1], Fy=FStart[1]; Cy<CStart[1]+CRange[1]; Cy++, Fy+=2 )
      for ( Cx=CStart[0], Fx=FStart[0]; Cx<CStart[0]+CRange[0]; Cx++, Fx+=2 )
      {
         CID = CID0 + Cz*Cdz + Cy*Cdy + Cx*Cdx;
         FID = FID0 + Fz*Fdz + Fy*Fdy + Fx*Fdx;

         LSlope = CData[CID    ] - CData[CID-Cdx];
         RSlope = CData[CID+Cdx] - CData[CID    ];
         if ( RSlope*LSlope <= (real)0.0 )   Slope_x = (real)0.0;
         else                                Slope_x = (real)0.25*( FABS(RSlope) < FABS(LSlope) ? RSlope:LSlope );

         if ( OppSign0thOrder  &&  CData[CID-Cdx]*CData[CID+Cdx] < (real)0.0 )   Slope_x = (real)0.0;

         LSlope = CData[CID    ] - CData[CID-Cdy];
         RSlope = CData[CID+Cdy] - CData[CID    ];
         if ( RSlope*LSlope <= (real)0.0 )   Slope_y = (real)0.0;
         else                                Slope_y = (real)0.25*( FABS(RSlope) < FABS(LSlope) ? RSlope:LSlope );

         if ( OppSign0thOrder  &&  CData[CID-Cdy]*CData[CID+Cdy] < (real)0.0 )   Slope_y = (real)0.0;

         LSlope = CData[CID    ] - CData[CID-Cdz];
         RSlope = CData[CID+Cdz] - CData[CID    ];
         if ( RSlope*LSlope <= (real)0.0 )   Slope_z = (real)0.0;
         else                                Slope_z = (real)0.25*( FABS(RSlope) < FABS(LSlope) ? RSlope:LSlope );

         if ( OppSign0thOrder  &&  CData[CID-Cdz]*CData[CID+Cdz] < (real)0.0 )   Slope_z = (real)0.0;


         FData[FID            ] = CData[CID] - Slope_z - Slope_y - Slope_x;
         FData[FID        +Fdx] = CData[CID] - Slope_z - Slope_y + Slope_x;
         FData[FID    +Fdy    ] = CData[CID] - Slope_z + Slope_y - Slope_x;
         FData[FID    +Fdy+Fdx] = CData[CID] - Slope_z + Slope_y + Slope_x;
         FData[FID+Fdz        ] = CData[CID] + Slope_z - Slope_y - Slope_x;
         FData[FID+Fdz    +Fdx] = CData[CID] + Slope_z - Slope_y + Slope_x;
         FData[FID+Fdz+Fdy    ] = CData[CID] + Slope_z + Slope_y - Slope_x;
         FData[FID+Fdz+Fdy+Fdx] = CData[CID] + Slope_z + Slope_y + Slope_x;
      }
   } // for (int v=0; v<NComp; v++)

} // FUNCTION : Int_MinMod1D
