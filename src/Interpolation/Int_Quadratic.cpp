#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Int_Quadratic
// Description :  Perform spatial interpolation based on the quadratic interpolation
//
// Note	       :  1. The spatial disribution is approximated by a parabola (quadratic polynomial) in each
//                   direction
//                   --> A quadratic polynomial, which passes through the CENTRAL values in three cells,
//                       is used to approximate the spatial distribution
//		  2. The interpolation result is neither conservative nor monotonic
//		  3. 3D interpolation is achieved by performing interpolation along x, y, and z directions
//		     in order
//		  4. The "Monotonic" option is used to ensure that the interpolation results are monotonic
//		     --> A slope limiter is adopted to ensure the monotonicity
//
// Parameter   :  CData	          : Input coarse-grid array
//		          CSize	          : Size of the CData array
//		          CStart          : (x,y,z) starting indices to perform interpolation on the CData array
//		          CRange          : Number of grids in each direction to perform interpolation
//		          FData	          : Output fine-grid array
//		          FStart          : (x,y,z) starting indcies to store the interpolation results
//		          NComp	          : Number of components in the CData and FData array
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                Monotonic       : Ensure that all interpolation results are monotonic
//                MonoCoeff       : Slope limiter coefficient for the option "Monotonic"
//                OppSign0thOrder : See Int_MinMod1D()
//-------------------------------------------------------------------------------------------------------
void Int_Quadratic( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
		            real FData[], const int FSize[3], const int FStart[3], const int NComp,
                    const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff,
                    const bool OppSign0thOrder )
{

// interpolation-scheme-dependent parameters
// ===============================================================================
// number of coarse-grid ghost zone
   const int CGhost = 1;

// interpolation coefficients
   const real R[ 1 + 2*CGhost ] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
   const real L[ 1 + 2*CGhost ] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };
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

// MonoCoeff/4
   const real MonoCoeff_4 = (real)0.25*MonoCoeff;

   real *CPtr   = CData;
   real *FPtr   = FData;
   real *TDataX = new real [ (CRange[2]+2*CGhost)*TdzX ];   // temporary array after x interpolation
   real *TDataY = new real [ (CRange[2]+2*CGhost)*TdzY ];   // temporary array after y interpolation

   int Idx_InL, Idx_InC, Idx_InR, Idx_Out;
   real LSlopeDh_4, RSlopeDh_4, SlopeDh_4, Sign, CDataMax, CDataMin;


   for (int v=0; v<NComp; v++)
   {
//    unwrap phase along x direction
#     if ( MODEL == ELBDM )
      if ( UnwrapPhase )
      {
         for (int k=CStart[2]-CGhost;    k<CStart[2]+CRange[2]+CGhost;  k++)
         for (int j=CStart[1]-CGhost;    j<CStart[1]+CRange[1]+CGhost;  j++)
         for (int i=CStart[0]-CGhost+1;  i<CStart[0]+CRange[0]+CGhost;  i++)
         {
            Idx_InC       = k*Cdz + j*Cdy + i*Cdx;
            Idx_InL       = Idx_InC - Cdx;
            CPtr[Idx_InC] = ELBDM_UnwrapPhase( CPtr[Idx_InL], CPtr[Idx_InC] );
         }
      }
#     endif


//    interpolation along x direction
      for (int In_z=CStart[2]-CGhost, Out_z=0;  In_z<CStart[2]+CRange[2]+CGhost;  In_z++, Out_z++)
      for (int In_y=CStart[1]-CGhost, Out_y=0;  In_y<CStart[1]+CRange[1]+CGhost;  In_y++, Out_y++)
      for (int In_x=CStart[0],        Out_x=0;  In_x<CStart[0]+CRange[0];         In_x++, Out_x+=2)
      {
         Idx_Out = Out_z*TdzX + Out_y*Tdy + Out_x*Tdx;
         Idx_InC =  In_z*Cdz  +  In_y*Cdy +  In_x*Cdx;
         Idx_InL = Idx_InC - Cdx;
         Idx_InR = Idx_InC + Cdx;

         TDataX[ Idx_Out       ] = L[0]*CPtr[Idx_InL] + L[1]*CPtr[Idx_InC] + L[2]*CPtr[Idx_InR];
         TDataX[ Idx_Out + Tdx ] = R[0]*CPtr[Idx_InL] + R[1]*CPtr[Idx_InC] + R[2]*CPtr[Idx_InR];

//       ensure monotonicity
         if ( Monotonic[v] )
         {
            LSlopeDh_4 = CPtr[Idx_InC] - CPtr[Idx_InL];
            RSlopeDh_4 = CPtr[Idx_InR] - CPtr[Idx_InC];

            if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )
            {
               if ( LSlopeDh_4 > (real)0.0 )
               {
                  CDataMax = CPtr[Idx_InR];
                  CDataMin = CPtr[Idx_InL];
               }
               else
               {
                  CDataMax = CPtr[Idx_InL];
                  CDataMin = CPtr[Idx_InR];
               }

               if (  FMAX( TDataX[Idx_Out], TDataX[Idx_Out+Tdx] ) > CDataMax  ||
                     FMIN( TDataX[Idx_Out], TDataX[Idx_Out+Tdx] ) < CDataMin    )
               {
                  SlopeDh_4   = (real)0.125*( CPtr[Idx_InR] - CPtr[Idx_InL] );
                  LSlopeDh_4 *= MonoCoeff_4;
                  RSlopeDh_4 *= MonoCoeff_4;
                  Sign        = SIGN( LSlopeDh_4 );

                  SlopeDh_4 *= Sign;
                  SlopeDh_4  = FMIN( Sign*LSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4  = FMIN( Sign*RSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4 *= Sign;

                  TDataX[ Idx_Out       ] = CPtr[Idx_InC] - SlopeDh_4;
                  TDataX[ Idx_Out + Tdx ] = CPtr[Idx_InC] + SlopeDh_4;
               }
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )

            else
            {
               TDataX[ Idx_Out       ] = CPtr[Idx_InC];
               TDataX[ Idx_Out + Tdx ] = CPtr[Idx_InC];
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 ) ... else
         } // if ( Monotonic[v] )

         if ( OppSign0thOrder  &&  CPtr[Idx_InL]*CPtr[Idx_InR] < (real)0.0 )
         {
            TDataX[ Idx_Out       ] = CPtr[Idx_InC];
            TDataX[ Idx_Out + Tdx ] = CPtr[Idx_InC];
         }
      } // k,j,i


//    unwrap phase along y direction
#     if ( MODEL == ELBDM )
      if ( UnwrapPhase )
      {
         for (int k=0;  k<CRange[2]+2*CGhost;  k++)
         for (int j=1;  j<CRange[1]+2*CGhost;  j++)
         for (int i=0;  i<2*CRange[0];         i++)
         {
            Idx_InC         = k*TdzX + j*Tdy + i*Tdx;
            Idx_InL         = Idx_InC - Tdy;
            TDataX[Idx_InC] = ELBDM_UnwrapPhase( TDataX[Idx_InL], TDataX[Idx_InC] );
         }
      }
#     endif


//    interpolation along y direction
      for (int InOut_z=0;             InOut_z<CRange[2]+2*CGhost;  InOut_z++)
      for (int In_y=CGhost, Out_y=0;  In_y   <CGhost+CRange[1];    In_y++, Out_y+=2)
      for (int InOut_x=0;             InOut_x<2*CRange[0];         InOut_x++)
      {
         Idx_Out = InOut_z*TdzY + Out_y*Tdy + InOut_x*Tdx;
         Idx_InC = InOut_z*TdzX +  In_y*Tdy + InOut_x*Tdx;
         Idx_InL = Idx_InC - Tdy;
         Idx_InR = Idx_InC + Tdy;

         TDataY[ Idx_Out       ] = L[0]*TDataX[Idx_InL] + L[1]*TDataX[Idx_InC] + L[2]*TDataX[Idx_InR];
         TDataY[ Idx_Out + Tdy ] = R[0]*TDataX[Idx_InL] + R[1]*TDataX[Idx_InC] + R[2]*TDataX[Idx_InR];

//       ensure monotonicity
         if ( Monotonic[v] )
         {
            LSlopeDh_4 = TDataX[Idx_InC] - TDataX[Idx_InL];
            RSlopeDh_4 = TDataX[Idx_InR] - TDataX[Idx_InC];

            if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )
            {
               if ( LSlopeDh_4 > (real)0.0 )
               {
                  CDataMax = TDataX[Idx_InR];
                  CDataMin = TDataX[Idx_InL];
               }
               else
               {
                  CDataMax = TDataX[Idx_InL];
                  CDataMin = TDataX[Idx_InR];
               }

               if (  FMAX( TDataY[Idx_Out], TDataY[Idx_Out+Tdy] ) > CDataMax  ||
                     FMIN( TDataY[Idx_Out], TDataY[Idx_Out+Tdy] ) < CDataMin    )
               {
                  SlopeDh_4   = (real)0.125*( TDataX[Idx_InR] - TDataX[Idx_InL] );
                  LSlopeDh_4 *= MonoCoeff_4;
                  RSlopeDh_4 *= MonoCoeff_4;
                  Sign        = SIGN( LSlopeDh_4 );

                  SlopeDh_4 *= Sign;
                  SlopeDh_4  = FMIN( Sign*LSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4  = FMIN( Sign*RSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4 *= Sign;

                  TDataY[ Idx_Out       ] = TDataX[Idx_InC] - SlopeDh_4;
                  TDataY[ Idx_Out + Tdy ] = TDataX[Idx_InC] + SlopeDh_4;
               }
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )

            else
            {
               TDataY[ Idx_Out       ] = TDataX[Idx_InC];
               TDataY[ Idx_Out + Tdy ] = TDataX[Idx_InC];
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 ) ... else
         } // if ( Monotonic[v] )

         if ( OppSign0thOrder  &&  TDataX[Idx_InL]*TDataX[Idx_InR] < (real)0.0 )
         {
            TDataY[ Idx_Out       ] = TDataX[Idx_InC];
            TDataY[ Idx_Out + Tdy ] = TDataX[Idx_InC];
         }
      } // k,j,i


//    unwrap phase along z direction
#     if ( MODEL == ELBDM )
      if ( UnwrapPhase )
      {
         for (int k=1;  k<CRange[2]+2*CGhost;  k++)
         for (int j=0;  j<2*CRange[1];         j++)
         for (int i=0;  i<2*CRange[0];         i++)
         {
            Idx_InC         = k*TdzY + j*Tdy + i*Tdx;
            Idx_InL         = Idx_InC - TdzY;
            TDataY[Idx_InC] = ELBDM_UnwrapPhase( TDataY[Idx_InL], TDataY[Idx_InC] );
         }
      }
#     endif


//    interpolation along z direction
      for (int In_z=CGhost, Out_z=FStart[2];  In_z<CGhost+CRange[2];  In_z++, Out_z+=2)
      for (int In_y=0,      Out_y=FStart[1];  In_y<2*CRange[1];       In_y++, Out_y++)
      for (int In_x=0,      Out_x=FStart[0];  In_x<2*CRange[0];       In_x++, Out_x++)
      {
         Idx_Out = Out_z*Fdz  + Out_y*Fdy + Out_x*Fdx;
         Idx_InC =  In_z*TdzY +  In_y*Tdy +  In_x*Tdx;
         Idx_InL = Idx_InC - TdzY;
         Idx_InR = Idx_InC + TdzY;

         FPtr[ Idx_Out       ] = L[0]*TDataY[Idx_InL] + L[1]*TDataY[Idx_InC] + L[2]*TDataY[Idx_InR];
         FPtr[ Idx_Out + Fdz ] = R[0]*TDataY[Idx_InL] + R[1]*TDataY[Idx_InC] + R[2]*TDataY[Idx_InR];

//       ensure monotonicity
         if ( Monotonic[v] )
         {
            LSlopeDh_4 = TDataY[Idx_InC] - TDataY[Idx_InL];
            RSlopeDh_4 = TDataY[Idx_InR] - TDataY[Idx_InC];

            if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )
            {
               if ( LSlopeDh_4 > (real)0.0 )
               {
                  CDataMax = TDataY[Idx_InR];
                  CDataMin = TDataY[Idx_InL];
               }
               else
               {
                  CDataMax = TDataY[Idx_InL];
                  CDataMin = TDataY[Idx_InR];
               }

               if (  FMAX( FPtr[Idx_Out], FPtr[Idx_Out+Fdz] ) > CDataMax  ||
                     FMIN( FPtr[Idx_Out], FPtr[Idx_Out+Fdz] ) < CDataMin    )
               {
                  SlopeDh_4   = (real)0.125*( TDataY[Idx_InR] - TDataY[Idx_InL] );
                  LSlopeDh_4 *= MonoCoeff_4;
                  RSlopeDh_4 *= MonoCoeff_4;
                  Sign        = SIGN( LSlopeDh_4 );

                  SlopeDh_4 *= Sign;
                  SlopeDh_4  = FMIN( Sign*LSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4  = FMIN( Sign*RSlopeDh_4, SlopeDh_4 );
                  SlopeDh_4 *= Sign;

                  FPtr[ Idx_Out       ] = TDataY[Idx_InC] - SlopeDh_4;
                  FPtr[ Idx_Out + Fdz ] = TDataY[Idx_InC] + SlopeDh_4;
               }
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 )

            else
            {
               FPtr[ Idx_Out       ] = TDataY[Idx_InC];
               FPtr[ Idx_Out + Fdz ] = TDataY[Idx_InC];
            } // if ( LSlopeDh_4*RSlopeDh_4 > (real)0.0 ) ... else
         } // if ( Monotonic[v] )

         if ( OppSign0thOrder  &&  TDataY[Idx_InL]*TDataY[Idx_InR] < (real)0.0 )
         {
            FPtr[ Idx_Out       ] = TDataY[Idx_InC];
            FPtr[ Idx_Out + Fdz ] = TDataY[Idx_InC];
         }
      } // k,j,i

      CPtr += CDisp;
      FPtr += FDisp;

   } // for (int v=0; v<NComp; v++)

   delete [] TDataX;
   delete [] TDataY;

} // FUNCTION : Int_Quadratic
