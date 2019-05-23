#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Int_WENO
// Description :  Perform spatial interpolation based on WENO
//
// Note        :  a. The slope at each grid is determined by the minimum slope between the right slope
//                   (difference between the right grid and itself) and the left slope (difference between
//                   itself and the left slope)
//                b. The slope is chosen to be zero if the right and left slopes have different signs
//                c. The interpolation result is BOTH conservative and monotonic
//		  d. 3D interpolation is achieved by performing interpolation along x, y, and z directions
//		     in order --> different from MINMOD1D
//
// Parameter   :  CData       : Input coarse-grid array
//                CSize       : Size of the CData array
//                CStart      : (x,y,z) starting indices to perform interpolation on the CData array
//                CRange      : Number of grids in each direction to perform interpolation
//                FData       : Output fine-grid array
//                FStart      : (x,y,z) starting indcies to store the interpolation results
//                NComp       : Number of components in the CData and FData array
//                UnwrapPhase : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//-------------------------------------------------------------------------------------------------------
void Int_WENO_O3( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData[], const int FSize[3], const int FStart[3], const int NComp,
                  const bool UnwrapPhase )
{
//k=3
// i+0.25
real cR[3][3] = 
{
 { (real)+6.145833e-01, (real)+5.208333e-01, (real)-1.354167e-01 },
 { (real)-1.354167e-01, (real)+1.020833e+00, (real)+1.145833e-01 },
 { (real)+1.145833e-01, (real)-4.791667e-01, (real)+1.364583e+00 },
};

real dR[3] = { 0.17572115384615383, 0.6001311188811189, 0.22414772727272728};

// i-0.25
real cL[3][3] =
{
 { (real)+1.364583e+00, (real)-4.791667e-01, (real)+1.145833e-01 },
 { (real)+1.145833e-01, (real)+1.020833e+00, (real)-1.354167e-01 },
 { (real)-1.354167e-01, (real)+5.208333e-01, (real)+6.145833e-01 },
};

real dL[3] = { 0.22414772727272728, 0.6001311188811189, 0.17572115384615383 };


// interpolation-scheme-dependent parameters
// ===============================================================================
// number of coarse-grid ghost zone
   const int CGhost = 3;
// ===============================================================================


// index stride of the coarse-grid input array
   const int Cdx    = 1;
   const int Cdy    = Cdx*CSize[0];
   const int Cdz    = Cdy*CSize[1];

// index stride of the temporary arrays storing the data after x and y interpolations
   const int Tdx    = 1;
   const int Tdy    = Tdx* CRange[0]*2;
   const int TdzX   = Tdy*(CRange[1]+2*CGhost);    // array after x interpolation
   const int TdzY   = Tdy* CRange[1]*2;            // array after y interpolation

// index stride of the fine-grid output array
   const int Fdx    = 1;
   const int Fdy    = Fdx*FSize[0];
   const int Fdz    = Fdy*FSize[1];

// index stride of different components
   const int CDisp  = CSize[0]*CSize[1]*CSize[2];
   const int FDisp  = FSize[0]*FSize[1]*FSize[2];

   real *CPtr   = CData;
   real *FPtr   = FData;
   real *TDataX = new real [ (CRange[2]+2*CGhost)*TdzX ];   // temporary array after x interpolation
   real *TDataY = new real [ (CRange[2]+2*CGhost)*TdzY ];   // temporary array after y interpolation

   int  Idx_InL3, Idx_InL2, Idx_InL1, Idx_InC, Idx_InR1, Idx_InR2, Idx_InR3, Idx_Out;
   real RStencil[CGhost], LStencil[CGhost];
   real Smoothness[3], AlphaL[3], OmegaL[3], AlphaR[3], OmegaR[3], SumAlphaL, SumAlphaR;
   real Epsilon = (real)1e-5;
   real A = (real)1.0833333333333333;
   real B = (real)0.25;

   for (int v=0; v<NComp; v++)
   {
//    interpolation along x direction
      for (int In_z=CStart[2]-CGhost, Out_z=0;  In_z<CStart[2]+CRange[2]+CGhost;  In_z++, Out_z++)
      for (int In_y=CStart[1]-CGhost, Out_y=0;  In_y<CStart[1]+CRange[1]+CGhost;  In_y++, Out_y++)
      for (int In_x=CStart[0],        Out_x=0;  In_x<CStart[0]+CRange[0];         In_x++, Out_x+=2)
      {
         Idx_Out  = Out_z*TdzX + Out_y*Tdy + Out_x*Tdx;
         Idx_InC  =  In_z*Cdz  +  In_y*Cdy +  In_x*Cdx;
         Idx_InL1 = Idx_InC  - Cdx;
         Idx_InL2 = Idx_InL1 - Cdx;
         Idx_InR1 = Idx_InC  + Cdx;
         Idx_InR2 = Idx_InR1 + Cdx;

         RStencil[0] = cR[0][0]*CPtr[Idx_InC ] + cR[0][1]*CPtr[Idx_InR1] + cR[0][2]*CPtr[Idx_InR2];
         RStencil[1] = cR[1][0]*CPtr[Idx_InL1] + cR[1][1]*CPtr[Idx_InC ] + cR[1][2]*CPtr[Idx_InR1];
         RStencil[2] = cR[2][0]*CPtr[Idx_InL2] + cR[2][1]*CPtr[Idx_InL1] + cR[2][2]*CPtr[Idx_InC ];

         LStencil[0] = cL[0][0]*CPtr[Idx_InC ] + cL[0][1]*CPtr[Idx_InL1] + cL[0][2]*CPtr[Idx_InL2];
         LStencil[1] = cL[1][0]*CPtr[Idx_InL1] + cL[1][1]*CPtr[Idx_InC ] + cL[1][2]*CPtr[Idx_InL1];
         LStencil[2] = cL[2][0]*CPtr[Idx_InL2] + cL[2][1]*CPtr[Idx_InL1] + cL[2][2]*CPtr[Idx_InC ];

         Smoothness[0] = A*SQR(           CPtr[Idx_InC ] - (real)2.0*CPtr[Idx_InR1] +           CPtr[Idx_InR2] )
                        +B*SQR( (real)3.0*CPtr[Idx_InC ] - (real)4.0*CPtr[Idx_InR1] +           CPtr[Idx_InR2] );
      
         Smoothness[1] = A*SQR(           CPtr[Idx_InL1] - (real)2.0*CPtr[Idx_InC ] +           CPtr[Idx_InR1] )
                        +B*SQR(           CPtr[Idx_InL1]                            -           CPtr[Idx_InR1] );

         Smoothness[2] = A*SQR(           CPtr[Idx_InL2] - (real)2.0*CPtr[Idx_InL1] +           CPtr[Idx_InC ] )
                        +B*SQR(           CPtr[Idx_InL2] - (real)4.0*CPtr[Idx_InL1] + (real)3.0*CPtr[Idx_InC ] );

         SumAlphaL = SumAlphaR = (real)0.0;

         for (int i=0;i<3;i++)
         {
           AlphaL[i] = dL[i]/SQR( Epsilon + Smoothness[i] );
           AlphaR[i] = dR[i]/SQR( Epsilon + Smoothness[i] );

           SumAlphaL += AlphaL[i];
           SumAlphaR += AlphaR[i];
         }
       
         TDataX[ Idx_Out ] = TDataX[ Idx_Out + Tdx ] = (real)0.0;
 
         for (int i=0;i<3;i++)
         {
           OmegaL[i] = AlphaL[i] / SumAlphaL;
           OmegaR[i] = AlphaR[i] / SumAlphaR;

           TDataX[ Idx_Out       ] += OmegaL[i]*LStencil[i];
           TDataX[ Idx_Out + Tdx ] += OmegaR[i]*RStencil[i];
         }

      } // for k,j,i



//    interpolation along y direction
      for (int InOut_z=0;             InOut_z<CRange[2]+2*CGhost;  InOut_z++)
      for (int In_y=CGhost, Out_y=0;  In_y   <CGhost+CRange[1];    In_y++, Out_y+=2)
      for (int InOut_x=0;             InOut_x<2*CRange[0];         InOut_x++)
      {
         Idx_Out  = InOut_z*TdzY + Out_y*Tdy + InOut_x*Tdx;
         Idx_InC  = InOut_z*TdzX +  In_y*Tdy + InOut_x*Tdx;
         Idx_InL1 = Idx_InC  - Tdy;
         Idx_InL2 = Idx_InL1 - Tdy;
         Idx_InR1 = Idx_InC  + Tdy;
         Idx_InR2 = Idx_InR1 + Tdy;

         RStencil[0] = cR[0][0]*TDataX[Idx_InC ] + cR[0][1]*TDataX[Idx_InR1] + cR[0][2]*TDataX[Idx_InR2];
         RStencil[1] = cR[1][0]*TDataX[Idx_InL1] + cR[1][1]*TDataX[Idx_InC ] + cR[1][2]*TDataX[Idx_InR1];
         RStencil[2] = cR[2][0]*TDataX[Idx_InL2] + cR[2][1]*TDataX[Idx_InL1] + cR[2][2]*TDataX[Idx_InC ];

         LStencil[0] = cL[0][0]*TDataX[Idx_InC ] + cL[0][1]*TDataX[Idx_InL1] + cL[0][2]*TDataX[Idx_InL2];
         LStencil[1] = cL[1][0]*TDataX[Idx_InL1] + cL[1][1]*TDataX[Idx_InC ] + cL[1][2]*TDataX[Idx_InL1];
         LStencil[2] = cL[2][0]*TDataX[Idx_InL2] + cL[2][1]*TDataX[Idx_InL1] + cL[2][2]*TDataX[Idx_InC ];

         Smoothness[0] = A*SQR(           TDataX[Idx_InC ] - (real)2.0*TDataX[Idx_InR1] +           TDataX[Idx_InR2] )
                        +B*SQR( (real)3.0*TDataX[Idx_InC ] - (real)4.0*TDataX[Idx_InR1] +           TDataX[Idx_InR2] );
      
         Smoothness[1] = A*SQR(           TDataX[Idx_InL1] - (real)2.0*TDataX[Idx_InC ] +           TDataX[Idx_InR1] )
                        +B*SQR(           TDataX[Idx_InL1]                              -           TDataX[Idx_InR1] );

         Smoothness[2] = A*SQR(           TDataX[Idx_InL2] - (real)2.0*TDataX[Idx_InL1] +           TDataX[Idx_InC ] )
                        +B*SQR(           TDataX[Idx_InL2] - (real)4.0*TDataX[Idx_InL1] + (real)3.0*TDataX[Idx_InC ] );

         SumAlphaL = SumAlphaR = (real)0.0;

         for (int i=0;i<3;i++)
         {
           AlphaL[i] = dL[i]/SQR( Epsilon + Smoothness[i] );
           AlphaR[i] = dR[i]/SQR( Epsilon + Smoothness[i] );

           SumAlphaL += AlphaL[i];
           SumAlphaR += AlphaR[i];
         }

         TDataY[ Idx_Out ] = TDataY[ Idx_Out + Tdy ] = (real)0.0;
        
         for (int i=0;i<3;i++)
         {
           OmegaL[i] = AlphaL[i] / SumAlphaL;
           OmegaR[i] = AlphaR[i] / SumAlphaR;

           TDataY[ Idx_Out       ] += OmegaL[i]*LStencil[i];
           TDataY[ Idx_Out + Tdy ] += OmegaR[i]*RStencil[i];
         }

      } // for k,j,i



//    interpolation along z direction
      for (int In_z=CGhost, Out_z=FStart[2];  In_z<CGhost+CRange[2];  In_z++, Out_z+=2)
      for (int In_y=0,      Out_y=FStart[1];  In_y<2*CRange[1];       In_y++, Out_y++)
      for (int In_x=0,      Out_x=FStart[0];  In_x<2*CRange[0];       In_x++, Out_x++)
      {
         Idx_Out  = Out_z*Fdz  + Out_y*Fdy + Out_x*Fdx;
         Idx_InC  =  In_z*TdzY +  In_y*Tdy +  In_x*Tdx;
         Idx_InL1 = Idx_InC  - TdzY;
         Idx_InL2 = Idx_InL1 - TdzY;
         Idx_InR1 = Idx_InC  + TdzY;
         Idx_InR2 = Idx_InR1 + TdzY;

         RStencil[0] = cR[0][0]*TDataY[Idx_InC ] + cR[0][1]*TDataY[Idx_InR1] + cR[0][2]*TDataY[Idx_InR2];
         RStencil[1] = cR[1][0]*TDataY[Idx_InL1] + cR[1][1]*TDataY[Idx_InC ] + cR[1][2]*TDataY[Idx_InR1];
         RStencil[2] = cR[2][0]*TDataY[Idx_InL2] + cR[2][1]*TDataY[Idx_InL1] + cR[2][2]*TDataY[Idx_InC ];

         LStencil[0] = cL[0][0]*TDataY[Idx_InC ] + cL[0][1]*TDataY[Idx_InL1] + cL[0][2]*TDataY[Idx_InL2];
         LStencil[1] = cL[1][0]*TDataY[Idx_InL1] + cL[1][1]*TDataY[Idx_InC ] + cL[1][2]*TDataY[Idx_InL1];
         LStencil[2] = cL[2][0]*TDataY[Idx_InL2] + cL[2][1]*TDataY[Idx_InL1] + cL[2][2]*TDataY[Idx_InC ];

         Smoothness[0] = A*SQR(           TDataY[Idx_InC ] - (real)2.0*TDataY[Idx_InR1] +           TDataY[Idx_InR2] )
                        +B*SQR( (real)3.0*TDataY[Idx_InC ] - (real)4.0*TDataY[Idx_InR1] +           TDataY[Idx_InR2] );
      
         Smoothness[1] = A*SQR(           TDataY[Idx_InL1] - (real)2.0*TDataY[Idx_InC ] +           TDataY[Idx_InR1] )
                        +B*SQR(           TDataY[Idx_InL1]                              -           TDataY[Idx_InR1] );

         Smoothness[2] = A*SQR(           TDataY[Idx_InL2] - (real)2.0*TDataY[Idx_InL1] +           TDataY[Idx_InC ] )
                        +B*SQR(           TDataY[Idx_InL2] - (real)4.0*TDataY[Idx_InL1] + (real)3.0*TDataY[Idx_InC ] );

         SumAlphaL = SumAlphaR = (real)0.0;

         for (int i=0;i<3;i++)
         {
           AlphaL[i] = dL[i]/SQR( Epsilon + Smoothness[i] );
           AlphaR[i] = dR[i]/SQR( Epsilon + Smoothness[i] );

           SumAlphaL += AlphaL[i];
           SumAlphaR += AlphaR[i];
         }

         FPtr[ Idx_Out ] = FPtr[ Idx_Out + Fdz ] = (real)0.0;
        
         for (int i=0;i<3;i++)
         {
           OmegaL[i] = AlphaL[i] / SumAlphaL;
           OmegaR[i] = AlphaR[i] / SumAlphaR;

           FPtr[ Idx_Out       ] += OmegaL[i]*LStencil[i];
           FPtr[ Idx_Out + Fdz ] += OmegaR[i]*RStencil[i];
         }

      } // for k,j,i

      CPtr += CDisp;
      FPtr += FDisp;

   } // for (int v=0; v<NComp; v++)

   delete [] TDataX;
   delete [] TDataY;

} // FUNCTION : Int_WENO
