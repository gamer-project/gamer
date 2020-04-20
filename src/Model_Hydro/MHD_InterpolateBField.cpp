#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )

static void GetFaceCoeff_F( real FaceCoeff[], const real *FInterface, const int Offset, const int Stride, real b[] );
static void GetFaceCoeff_C( real FaceCoeff[], const real *CData, const int Offset, const int Stride[2],
                            const IntScheme_t IntScheme, const bool Monotonic, const real MonoCoeff );
static void GetCellCoeff( real CellCoeff[][11], const real FaceCoeff[][4] );

static inline real GetCellBx( const real CellCoeff[], const real y, const real z );
static inline real GetCellBy( const real CellCoeff[], const real x, const real z );
static inline real GetCellBz( const real CellCoeff[], const real x, const real y );
static inline real GetFaceB_C( const real FaceCoeff[], const real dr0, const real dr1 );

static real (*GetCellB[3])( const real [], const real, const real ) = { GetCellBx, GetCellBy, GetCellBz };


// symbolic constants used for the polynomial coefficient array CellCoeff[]
// --> note that they are in a different order than Balsara 2001, Eqs. (4.6)-(4.8)
#define A0     0
#define Ax     1
#define Ay     2
#define Az     3
#define Axx    4
#define Axy    5
#define Ayz    6
#define Axz    7
#define Axyz   8
#define Axxy   9
#define Axxz  10

#define B0     0
#define Bx     1
#define By     2
#define Bz     3
#define Byy    4
#define Bxy    5
#define Byz    6
#define Bxz    7
#define Bxyz   8
#define Byyz   9
#define Bxyy  10

#define C0     0
#define Cx     1
#define Cy     2
#define Cz     3
#define Czz    4
#define Cxy    5
#define Cyz    6
#define Cxz    7
#define Cxyz   8
#define Cxzz   9
#define Cyzz  10




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_InterpolateBField
// Description :  Perform divergence-preserving interpolation on the magnetic field
//
// Note        :  1. Ref: Balsara 2001, Journal of Computational Physics 174, 614
//                2. Magnetic field is defined at the face center
//                3. Use the input parameter "IntScheme" to determine the adopted interpolation scheme
//                4. Must apply to all three magnetic components at once
//                   --> To preserve the divergence-free condition
//                5. Strictly speaking, the adopted scheme is divergence-preserving instead of divergence-free
//                   --> In other words, the output B field will be divergence-free only if the input B field
//                       is divergence-free
//                6. Assuming (0th/1th/2th) B field array index = (x/y/z) B component
//                   --> In other words, it assumes MAGX/Y/Z = 0/1/2
//
// Parameter   :  CData      : Input coarse-grid array
//                CSize      : Size of CData[]
//                CStart     : (x,y,z) starting indices to perform interpolation on CData[]
//                CRange     : Number of coarse cells along each direction to perform interpolation
//                FData      : Output fine-grid array
//                FSize      : Size of FData[]
//                FStart     : (x,y,z) starting indices to store the interpolation results
//                FInterface : B field to be fixed on the coarse-fine interfaces
//                             --> Set "FInterface[sib] == NULL" for coarse-coarse interfaces
//                IntScheme  : Interpolation scheme
//                Monotonic  : Ensure that all interpolation results are monotonic
//
// Return      :  FData[]
//-------------------------------------------------------------------------------------------------------
void MHD_InterpolateBField( const real **CData, const int CSize[3][3], const int CStart[3][3], const int CRange[3],
                                  real **FData, const int FSize[3][3], const int FStart[3][3],
                            const real *FInterface[6], const IntScheme_t IntScheme, const bool Monotonic )
{

// index stride of CData[] and FData[]: [3][3]=[B component][normal/transverse 1/transverse 2]
   const int dC[3][3]    = { {                       1, CSize[0][0], CSize[0][0]*CSize[0][1] },
                             {             CSize[1][0],           1, CSize[1][0]*CSize[1][1] },
                             { CSize[2][0]*CSize[2][1],           1,             CSize[2][0] } };
   const int dF[3][3]    = { {                       1, FSize[0][0], FSize[0][0]*FSize[0][1] },
                             {             FSize[1][0],           1, FSize[1][0]*FSize[1][1] },
                             { FSize[2][0]*FSize[2][1],           1,             FSize[2][0] } };

// size of FInterface[]
   const int FIntSize[3] = { 2*CRange[0], 2*CRange[1], 2*CRange[2] };

// for determining whether a target coarse cell is adjacent to the index boundaries
   const int IdxLimit[6] = { 0, CRange[0]-1, 0, CRange[1]-1, 0, CRange[2]-1 };

   const real p1_4       = (real)1.0/4.0;
   const real m1_4       = -p1_4;


   real FaceCoeff[6][4];   // polynomial coefficients in Balsara 2001, Eqs. (4.1)-(4.3)
   real CellCoeff[3][11];  // polynomial coefficients in Balsara 2001, Eqs. (4.6)-(4.8)
   real FaceB_F  [6][4];   // B field on the C-F interfaces along different directions

// C/F/I = indices for CData[]/FData[]/FInterface[]
   int i_C[3], j_C[3], k_C[3], i_F[3], j_F[3], k_F[3], i_I, j_I, k_I, IJK[3];


// iterate over all coarse cells to be interpolated
   for (IJK[2]=0; IJK[2]<CRange[2]; IJK[2]++)
   {
      k_I = 2*IJK[2];
      for (int v=0; v<3; v++) {
         k_C[v] = IJK[2] + CStart[v][2];
         k_F[v] = k_I    + FStart[v][2];
      }

   for (IJK[1]=0; IJK[1]<CRange[1]; IJK[1]++)
   {
      j_I = 2*IJK[1];
      for (int v=0; v<3; v++) {
         j_C[v] = IJK[1] + CStart[v][1];
         j_F[v] = j_I    + FStart[v][1];
      }

   for (IJK[0]=0; IJK[0]<CRange[0]; IJK[0]++)
   {
      i_I = 2*IJK[0];
      for (int v=0; v<3; v++) {
         i_C[v] = IJK[0] + CStart[v][0];
         i_F[v] = i_I    + FStart[v][0];
      }

//    1. set the polynomial coefficients on six cell faces
      for (int f=0; f<6; f++)
      {
         const int d = f/2;   // spatial direction

//       C-F interface
         if ( IJK[d] == IdxLimit[f]  &&  FInterface[f] != NULL )
         {
            int Offset, Stride;

            switch ( f )
            {
               case 0: case 1:   Stride = FIntSize[1];   Offset = k_I*Stride + j_I;    break;
               case 2: case 3:   Stride = FIntSize[0];   Offset = k_I*Stride + i_I;    break;
               case 4: case 5:   Stride = FIntSize[0];   Offset = j_I*Stride + i_I;    break;

               default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "f", f );
            }

            GetFaceCoeff_F( FaceCoeff[f], FInterface[f], Offset, Stride, FaceB_F[f] );
         }

//       C-C interface
         else
         {
            const int Stride[2] = { dC[d][1], dC[d][2] };
            const int Offset    = IDX321( i_C[d], j_C[d], k_C[d], CSize[d][0], CSize[d][1] ) + (f&1)*dC[d][0];

            GetFaceCoeff_C( FaceCoeff[f], CData[d], Offset, Stride, IntScheme, Monotonic, INT_MONO_COEFF );
         } // if ( IJK[d] == IdxLimit[f]  &&  FInterface[f] != NULL ) ... else ...
      } // for (int f=0; f<6; f++)


//    2. set the polynomial coefficients of the entire cell
      GetCellCoeff( CellCoeff, FaceCoeff );


//    3. fill in FData[] using the interpolation polynomial
      for (int v=0; v<NCOMP_MAG; v++)
      {
         const int Offset_F0   = IDX321( i_F[v], j_F[v], k_F[v], FSize[v][0], FSize[v][1] );
         const int Offset_F[4] = { Offset_F0,
                                   Offset_F0 + dF[v][1],
                                   Offset_F0 + dF[v][2],
                                   Offset_F0 + dF[v][1] + dF[v][2] };
         int f, Offset_Normal;

//       3-1. left face
         f = 2*v;

//       C-F interface
         if ( IJK[v] == IdxLimit[f]  &&  FInterface[f] != NULL )
         {
            FData[v][ Offset_F[0] ] = FaceB_F[f][0];
            FData[v][ Offset_F[1] ] = FaceB_F[f][1];
            FData[v][ Offset_F[2] ] = FaceB_F[f][2];
            FData[v][ Offset_F[3] ] = FaceB_F[f][3];
         }

//       C-C interface
         else
         {
            FData[v][ Offset_F[0] ] = GetFaceB_C( FaceCoeff[f], m1_4, m1_4 );
            FData[v][ Offset_F[1] ] = GetFaceB_C( FaceCoeff[f], p1_4, m1_4 );
            FData[v][ Offset_F[2] ] = GetFaceB_C( FaceCoeff[f], m1_4, p1_4 );
            FData[v][ Offset_F[3] ] = GetFaceB_C( FaceCoeff[f], p1_4, p1_4 );
         }


//       3-2. central face
         Offset_Normal = dF[v][0];

         FData[v][ Offset_F[0] + Offset_Normal ] = GetCellB[v]( CellCoeff[v], m1_4, m1_4 );
         FData[v][ Offset_F[1] + Offset_Normal ] = GetCellB[v]( CellCoeff[v], p1_4, m1_4 );
         FData[v][ Offset_F[2] + Offset_Normal ] = GetCellB[v]( CellCoeff[v], m1_4, p1_4 );
         FData[v][ Offset_F[3] + Offset_Normal ] = GetCellB[v]( CellCoeff[v], p1_4, p1_4 );


//       3-3. right face
//       --> only necessary on the rightmost face of the target coarse region to avoid redundant assignment
         Offset_Normal = 2*dF[v][0];
         f ++;

         if ( IJK[v] == IdxLimit[f] )
         {
//          C-F interface
            if ( FInterface[f] != NULL )
            {
               FData[v][ Offset_F[0] + Offset_Normal ] = FaceB_F[f][0];
               FData[v][ Offset_F[1] + Offset_Normal ] = FaceB_F[f][1];
               FData[v][ Offset_F[2] + Offset_Normal ] = FaceB_F[f][2];
               FData[v][ Offset_F[3] + Offset_Normal ] = FaceB_F[f][3];
            }

//          C-C interface
            else
            {
               FData[v][ Offset_F[0] + Offset_Normal ] = GetFaceB_C( FaceCoeff[f], m1_4, m1_4 );
               FData[v][ Offset_F[1] + Offset_Normal ] = GetFaceB_C( FaceCoeff[f], p1_4, m1_4 );
               FData[v][ Offset_F[2] + Offset_Normal ] = GetFaceB_C( FaceCoeff[f], m1_4, p1_4 );
               FData[v][ Offset_F[3] + Offset_Normal ] = GetFaceB_C( FaceCoeff[f], p1_4, p1_4 );
            }
         }
      } // for (int v=0; v<NCOMP_MAG; v++)

   }}} // IJK loop

} // FUNCTION : MHD_InterpolateBField



//-------------------------------------------------------------------------------------------------------
// Function    :  GetFaceCoeff_F
// Description :  Evaluate the polynomial coefficients defined in Balsara 2001, Eqs. (4.1)-(4.3)
//                using Eq. (4.5)
//
// Note        :  1. Apply to the coarse-fine interface
//                2. Also return the B field on the target cell face in b[]
//                   --> For copying data to FData[] directly in MHD_InterpolateBField()
//
// Parameter   :  FaceCoeff  : Polynomial coefficients to be returned
//                FInterface : B field on the target coarse-fine interface
//                Offset     : Array offset of FInterface[]
//                Stride     : Array stride of FInterface[]
//                b          : Subset of FInterface[] on the target cell face
//
// Return      :  FaceCoeff[], b
//-------------------------------------------------------------------------------------------------------
void GetFaceCoeff_F( real FaceCoeff[], const real *FInterface, const int Offset, const int Stride, real b[] )
{

   b[0] = FInterface[ Offset              ];
   b[1] = FInterface[ Offset          + 1 ];
   b[2] = FInterface[ Offset + Stride     ];
   b[3] = FInterface[ Offset + Stride + 1 ];

   const real b10  = b[1] - b[0];
   const real b20  = b[2] - b[0];
   const real b31  = b[3] - b[1];
   const real b32  = b[3] - b[2];

   FaceCoeff[0] = (real)0.25*( b[0] + b[1] + b[2] + b[3] );
   FaceCoeff[1] = b10 + b32;
   FaceCoeff[2] = b20 + b31;
   FaceCoeff[3] = (real)2.0*(  ( b32 - b10 ) + ( b31 - b20 )  );

} // FUNCTION : GetFaceCoeff_F



//-------------------------------------------------------------------------------------------------------
// Function    :  GetFaceCoeff_C
// Description :  Evaluate the polynomial coefficients defined in Balsara 2001, Eqs. (4.1)-(4.3)
//                using Eqs. (3.18)-(3.23)
//
// Note        :  1. Apply to the coarse-coarse interface
//
// Parameter   :  FaceCoeff : Polynomial coefficients to be returned
//                Offset    : Array offset of CData[]
//                Stride    : Array stride of CData[]
//                IntScheme : Interpolation scheme
//                            --> INT_MINMOD1D : MinMod
//                                INT_VANLEER  : vanLeer
//                                INT_CQUAD    : conservative quadratic
//                                INT_CQUAR    : conservative quartic
//                Monotonic : Ensure that all interpolation results are monotonic
//                MonoCoeff : Slope limiter coefficient for the option "Monotonic"
//
// Return      :  FaceCoeff[]
//-------------------------------------------------------------------------------------------------------
void GetFaceCoeff_C( real FaceCoeff[], const real *CData, const int Offset, const int Stride[2],
                     const IntScheme_t IntScheme, const bool Monotonic, const real MonoCoeff )
{

   real Cen;
   real L1[2], L2[2], R1[2], R2[2];
   real RSlope, LSlope, Sign;


   Cen = CData[ Offset ];

   for (int t=0; t<2; t++)
   {
      L1[t] = CData[ Offset - 1*Stride[t] ];
      R1[t] = CData[ Offset + 1*Stride[t] ];

      if ( IntScheme == INT_CQUAR ) {
      L2[t] = CData[ Offset - 2*Stride[t] ];
      R2[t] = CData[ Offset + 2*Stride[t] ];
      }
   }

   FaceCoeff[0] = Cen;
   FaceCoeff[3] = (real)0.0;  // cross derivative is always set to zero on a C-C interface


   switch ( IntScheme )
   {
      case INT_MINMOD1D :
      {
         for (int t=0; t<2; t++)
         {
            const int idx = t + 1;

            LSlope = Cen - L1[t];
            RSlope = R1[t] - Cen;

            if ( RSlope*LSlope <= (real)0.0 )   FaceCoeff[idx] = (real)0.0;
            else                                FaceCoeff[idx] = ( FABS(RSlope) < FABS(LSlope) ) ? RSlope : LSlope;
         }
      } // INT_MINMOD1D
      break;


      case INT_VANLEER :
      {
         for (int t=0; t<2; t++)
         {
            const int idx = t + 1;

            LSlope = Cen - L1[t];
            RSlope = R1[t] - Cen;

            if ( RSlope*LSlope <= (real)0.0 )   FaceCoeff[idx] = (real)0.0;
            else                                FaceCoeff[idx] = MonoCoeff*RSlope*LSlope/(LSlope+RSlope);
         }
      } // INT_VANLEER
      break;


      case INT_CQUAD :
      {
         for (int t=0; t<2; t++)
         {
            const int idx = t + 1;

            FaceCoeff[idx] = (real)0.5*( R1[t] - L1[t] );

            if ( Monotonic )
            {
               LSlope = Cen - L1[t];
               RSlope = R1[t] - Cen;

               if ( LSlope*RSlope > (real)0.0 )
               {
                  LSlope *= MonoCoeff;
                  RSlope *= MonoCoeff;
                  Sign    = SIGN( LSlope );

                  FaceCoeff[idx] *= Sign;
                  FaceCoeff[idx]  = FMIN( Sign*LSlope, FaceCoeff[idx] );
                  FaceCoeff[idx]  = FMIN( Sign*RSlope, FaceCoeff[idx] );
                  FaceCoeff[idx] *= Sign;
               }
               else
                  FaceCoeff[idx] = (real)0.0;
            } // if ( Monotonic )
         } // for (int t=0; t<2; t++)
      } // INT_CQUAD
      break;


      case INT_CQUAR :
      {
         const real IntCoeff[2] = { 3.0/32.0, 22.0/32.0 };

         for (int t=0; t<2; t++)
         {
            const int idx = t + 1;

            FaceCoeff[idx] = -IntCoeff[0]*R2[t] + IntCoeff[1]*R1[t] - IntCoeff[1]*L1[t] + IntCoeff[0]*L2[t];

            if ( Monotonic )
            {
               LSlope = Cen - L1[t];
               RSlope = R1[t] - Cen;

               if ( LSlope*RSlope > (real)0.0 )
               {
                  if ( FaceCoeff[idx]*LSlope < (real)0.0 )  FaceCoeff[idx] = (real)0.5*( R1[t] - L1[t] );

                  LSlope *= MonoCoeff;
                  RSlope *= MonoCoeff;
                  Sign    = SIGN( LSlope );

                  FaceCoeff[idx] *= Sign;
                  FaceCoeff[idx]  = FMIN( Sign*LSlope, FaceCoeff[idx] );
                  FaceCoeff[idx]  = FMIN( Sign*RSlope, FaceCoeff[idx] );
                  FaceCoeff[idx] *= Sign;
               }
               else
                  FaceCoeff[idx] = (real)0.0;
            } // if ( Monotonic )
         } // for (int t=0; t<2; t++)
      } // INT_CQUAR
      break;


      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "IntScheme", IntScheme );
         exit(1);
   } // switch ( IntScheme )

} // FUNCTION : GetFaceCoeff_C



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCellCoeff
// Description :  Evaluate the polynomial coefficients defined in Balsara 2001, Eqs. (4.6)-(4.7)
//                using Eqs. (3.32)-(3.38) and (4.15)-(4.18)
//
// Note        :  1. Must apply to all three magnetic components at once
//
// Parameter   :  CellCoeff : Polynomial coefficients to be returned
//                FaceCoeff : Polynomial coefficients on the six cell faces prepared by either
//                            GetFaceCoeff_C() or GetFaceCoeff_F()
//
// Return      :  CellCoeff[]
//-------------------------------------------------------------------------------------------------------
void GetCellCoeff( real CellCoeff[][11], const real FaceCoeff[][4] )
{

   const int A  = MAGX;
   const int B  = MAGY;
   const int C  = MAGZ;

   const int ZE = 0;    // FaceCoeff[][*]: 0th/derivative i/derivative j/cross derivative
   const int DI = 1;
   const int DJ = 2;
   const int CR = 3;

   const int XM = 0;    // FaceCoeff[*][]: for example, XM = B_{x}^{-}
   const int XP = 1;
   const int YM = 2;
   const int YP = 3;
   const int ZM = 4;
   const int ZP = 5;

   real Axyz_16, Bxyz_16, Cxyz_16;


   CellCoeff[A][Ax  ] =           ( FaceCoeff[XP][ZE] - FaceCoeff[XM][ZE] );
   CellCoeff[B][By  ] =           ( FaceCoeff[YP][ZE] - FaceCoeff[YM][ZE] );
   CellCoeff[C][Cz  ] =           ( FaceCoeff[ZP][ZE] - FaceCoeff[ZM][ZE] );

   CellCoeff[A][Ayz ] = (real)0.5*( FaceCoeff[XP][CR] + FaceCoeff[XM][CR] );
   CellCoeff[B][Bxz ] = (real)0.5*( FaceCoeff[YP][CR] + FaceCoeff[YM][CR] );
   CellCoeff[C][Cxy ] = (real)0.5*( FaceCoeff[ZP][CR] + FaceCoeff[ZM][CR] );

   CellCoeff[A][Axyz] =           ( FaceCoeff[XP][CR] - FaceCoeff[XM][CR] );
   CellCoeff[B][Bxyz] =           ( FaceCoeff[YP][CR] - FaceCoeff[YM][CR] );
   CellCoeff[C][Cxyz] =           ( FaceCoeff[ZP][CR] - FaceCoeff[ZM][CR] );

   CellCoeff[A][Axxy] = -(real)0.25*CellCoeff[C][Cxyz];
   CellCoeff[A][Axxz] = -(real)0.25*CellCoeff[B][Bxyz];
   CellCoeff[B][Byyz] = -(real)0.25*CellCoeff[A][Axyz];
   CellCoeff[B][Bxyy] =             CellCoeff[A][Axxy];
   CellCoeff[C][Cyzz] =             CellCoeff[B][Byyz];
   CellCoeff[C][Cxzz] =             CellCoeff[A][Axxz];

   CellCoeff[A][Axy ] =           ( FaceCoeff[XP][DI] - FaceCoeff[XM][DI] );
   CellCoeff[A][Axz ] =           ( FaceCoeff[XP][DJ] - FaceCoeff[XM][DJ] );
   CellCoeff[B][Bxy ] =           ( FaceCoeff[YP][DI] - FaceCoeff[YM][DI] );
   CellCoeff[B][Byz ] =           ( FaceCoeff[YP][DJ] - FaceCoeff[YM][DJ] );
   CellCoeff[C][Cxz ] =           ( FaceCoeff[ZP][DI] - FaceCoeff[ZM][DI] );
   CellCoeff[C][Cyz ] =           ( FaceCoeff[ZP][DJ] - FaceCoeff[ZM][DJ] );

   Axyz_16 = (real)0.0625*CellCoeff[A][Axyz];
   Bxyz_16 = (real)0.0625*CellCoeff[B][Bxyz];
   Cxyz_16 = (real)0.0625*CellCoeff[C][Cxyz];

   CellCoeff[A][Ay  ] = (real)0.5*( FaceCoeff[XP][DI] + FaceCoeff[XM][DI] ) + Cxyz_16;
   CellCoeff[A][Az  ] = (real)0.5*( FaceCoeff[XP][DJ] + FaceCoeff[XM][DJ] ) + Bxyz_16;
   CellCoeff[B][Bx  ] = (real)0.5*( FaceCoeff[YP][DI] + FaceCoeff[YM][DI] ) + Cxyz_16;
   CellCoeff[B][Bz  ] = (real)0.5*( FaceCoeff[YP][DJ] + FaceCoeff[YM][DJ] ) + Axyz_16;
   CellCoeff[C][Cx  ] = (real)0.5*( FaceCoeff[ZP][DI] + FaceCoeff[ZM][DI] ) + Bxyz_16;
   CellCoeff[C][Cy  ] = (real)0.5*( FaceCoeff[ZP][DJ] + FaceCoeff[ZM][DJ] ) + Axyz_16;

   CellCoeff[A][Axx ] = -(real)0.5*( CellCoeff[B][Bxy] + CellCoeff[C][Cxz] );
   CellCoeff[B][Byy ] = -(real)0.5*( CellCoeff[A][Axy] + CellCoeff[C][Cyz] );
   CellCoeff[C][Czz ] = -(real)0.5*( CellCoeff[A][Axz] + CellCoeff[B][Byz] );

   CellCoeff[A][A0  ] = (real)0.5*( FaceCoeff[XP][ZE] + FaceCoeff[XM][ZE] ) - (real)0.25*CellCoeff[A][Axx];
   CellCoeff[B][B0  ] = (real)0.5*( FaceCoeff[YP][ZE] + FaceCoeff[YM][ZE] ) - (real)0.25*CellCoeff[B][Byy];
   CellCoeff[C][C0  ] = (real)0.5*( FaceCoeff[ZP][ZE] + FaceCoeff[ZM][ZE] ) - (real)0.25*CellCoeff[C][Czz];

} // FUNCTION : GetCellCoeff



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCellBx/y/z
// Description :  Get the interpolated Bx/y/z on the x=0/y=0/z=0 plane of a target coarse cell using
//                Eq. (4.6)/(4.7)/(4.8)
//
// Note        :  1. Only applicable to the x=0/y=0/z=0 plane since all x-/y-/z-dependent terms in
//                   Eq. (4.6)/(4.7)/(4.8) have been ignored
//                2. To improve performance, we do not use this function on the coarse cell faces
//                   --> C-C interface: call GetFaceB_C()
//                       C-F interface: directly copy data from FInterface[]
//
// Parameter   :  CellCoeff : Polynomial coefficients of Bx/y/z in a coarse cell
//                x/y/z     : Relative coordinate to the coarse cell center
//                            --> Normalized to the coarse cell width
//
// Return      :  CellBx/y/z
//-------------------------------------------------------------------------------------------------------
inline real GetCellBx( const real CellCoeff[], const real y, const real z )
{

   const real CellBx = CellCoeff[A0] + CellCoeff[Ay]*y + CellCoeff[Az]*z + CellCoeff[Ayz]*y*z;

   return CellBx;

} // FUNCTION : GetCellBx



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCellBy
// Description :  See GetCellBx()
//-------------------------------------------------------------------------------------------------------
inline real GetCellBy( const real CellCoeff[], const real x, const real z )
{

   const real CellBy = CellCoeff[B0] + CellCoeff[Bx]*x + CellCoeff[Bz]*z + CellCoeff[Bxz]*x*z;

   return CellBy;

} // FUNCTION : GetCellBy



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCellBz
// Description :  See GetCellBx()
//-------------------------------------------------------------------------------------------------------
inline real GetCellBz( const real CellCoeff[], const real x, const real y )
{

   const real CellBz = CellCoeff[C0] + CellCoeff[Cx]*x + CellCoeff[Cy]*y + CellCoeff[Cxy]*x*y;

   return CellBz;

} // FUNCTION : GetCellBx



//-------------------------------------------------------------------------------------------------------
// Function    :  GetFaceB_C
// Description :  Get the interpolated B field on the coarse cell face using Eqs. (3.24)-(3.26)
//
// Note        :  1. It's equivalent to calling GetCellB() but should be faster
//                2. Apply to C-C interfaces only
//                   --> Not applicable to F-C interfaces since the cross derivative term is ignored here
//
// Parameter   :  FaceCoeff : Polynomial coefficients on the target cell face
//                dr0/1     : Relative coordinate along the FaceCoeff[1/2] direction
//                            --> Normalized to the coarse cell width
//
// Return      :  FaceB
//-------------------------------------------------------------------------------------------------------
inline real GetFaceB_C( const real FaceCoeff[], const real dr0, const real dr1 )
{

   const real FaceB = FaceCoeff[0] + FaceCoeff[1]*dr0 + FaceCoeff[2]*dr1;

   return FaceB;

} // FUNCTION : GetFaceB_C



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
