#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == SOR  &&  defined USE_PSOLVER_10TO14 )



#define POT_NXT_F    ( PATCH_SIZE+2*POT_GHOST_SIZE           )
#define POT_PAD      ( WARP_SIZE/2 - (POT_NXT_F*2%WARP_SIZE) )
#define POT_NTHREAD  ( RHO_NXT*RHO_NXT*POT_BLOCK_SIZE_Z/2    )
#define POT_USELESS  ( POT_GHOST_SIZE%2                      )


/************************************************************
  Many optimization options for SOR are defined in CUPOT.h
************************************************************/


// variables reside in constant memory
#include "CUPOT_SetConstMem_PoissonSolver.cu"


// parallel reduction routine
#define RED_NTHREAD  POT_NTHREAD
#define RED_SUM

#ifdef SOR_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif


// checks
#ifdef SOR_USE_PADDING
#  if ( WARP_SIZE != 32 )
#     error : ERROR : WARP_SIZE != 32 !!
#  endif
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_PoissonSolver_SOR_10to14cube
// Description :  GPU Poisson solver using the SOR scheme
//
// Note        :  1. Work for POT_GHOST_SIZE = 1, 2, 3 <--> POT_NXT_F = 10, 12, 14
//                   --> For compute capabilities >= 2.0 (which as 48 KB shared memory), it also works
//                       with POT_GHOST_SIZE = 4, 5
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                3. Each patch requires about 3.1*10^6 FLOPS (include the gravity solver)
//                   --> 133 GFLOPS is achieved in one C2050 GPU
//                4. Reference : Numerical Recipes, Chapter 20.5
//                5. Chester Cheng has implemented the SOR_USE_SHUFFLE and SOR_USE_PADDING optimizations, which
//                   greatly improve performance for the case POT_GHOST_SIZE == 5
//                6. Typically, the number of iterations required to reach round-off errors is 20 ~ 25 (single precision)
//
// Padding     :  Below shows how bank conflict is eliminated by padding.
//
//                Example constants :
//                      POT_NXT_F = 18                       // The number of floating point elements per row
//                      POT_PAD   = 16 - (18 * 2 % 32) = 12  // number of floating point elements that needs to be added
//                                                           // within thread groups
//
//                We now show how shared memory (s_FPot array) is accessed by a warp in residual evaluation.
//
//                Before Padding:
//                Thread number   |  Accessed shared memory bank
//                      00 ~ 07   |    | 01 |    | 03 |    | 05 |    | 07 |    | 09 |    | 11 |    | 13 |    | 15 |    |    |
//                      08 ~ 15   |    |    | 02 |    | 04 |    | 06 |    | 08 |    | 10 |    | 12 |    | 14 |    | 16 |    |
//                      16 ~ 23   |    | 05 |    | 07 |    | 09 |    | 11 |    | 13 |    | 15 |    | 17 |    | 19 |    |    |
//                      24 ~ 31   |    |    | 06 |    | 08 |    | 10 |    | 12 |    | 14 |    | 16 |    | 18 |    | 20 |    |
//
//                After Padding:
//                Thread number   |  Accessed shared memory bank
//                      00 ~ 07   |    | 01 |    | 03 |    | 05 |    | 07 |    | 09 |    | 11 |    | 13 |    | 15 |    |    |
//                      08 ~ 15   |    |    | 02 |    | 04 |    | 06 |    | 08 |    | 10 |    | 12 |    | 14 |    | 16 |    |
//                                ----------------- PAD 12 FLOATING POINTS HERE !!!!! ---------------------------------------
//                      16 ~ 23   |    | 17 |    | 19 |    | 21 |    | 23 |    | 25 |    | 27 |    | 29 |    | 31 |    |    |
//                      24 ~ 31   |    |    | 18 |    | 20 |    | 22 |    | 24 |    | 26 |    | 28 |    | 30 |    | 00 |    |
//
//
//                Additional Notes for Padding:
//                      1. When threads 08 ~ 15 access the elements below them (+y direction), we have to skip the padded
//                         elements. Same for when threads 16~23 access the elements above them (-y direction).
//                      2. For every warp we need to pad #PAD_POT floating point elements. Each xy plane has 4 warps working
//                         on it, so for each xy plane we need to pad #4*PAD_POT floating point elements.
//
//
// Parameter   :  g_Rho_Array       : Global memory array to store the input density
//                g_Pot_Array_In    : Global memory array storing the input "coarse-grid" potential for
//                                    interpolation
//                g_Pot_Array_Out   : Global memory array to store the output potential
//                Min_Iter          : Minimum # of iterations for SOR
//                Max_Iter          : Maximum # of iterations for SOR
//                Omega_6           : Omega / 6
//                Const             : (Coefficient in front of the RHS in the Poisson eq.) / dh^2
//                IntScheme         : Interpolation scheme for potential
//                                    --> currently supported schemes include
//                                        INT_CQUAD : conservative quadratic interpolation
//                                        INT_QUAD  : quadratic interpolation
//---------------------------------------------------------------------------------------------------
__global__ void CUPOT_PoissonSolver_SOR_10to14cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6,
                                                    const real Const, const IntScheme_t IntScheme )
{

   const uint bid       = blockIdx.x;
   const uint tid_x     = threadIdx.x;
   const uint tid_y     = threadIdx.y;
   const uint tid_z     = threadIdx.z;
   const uint bdim_x    = blockDim.x;
   const uint bdim_y    = blockDim.y;
   const uint bdim_z    = blockDim.z;
   const uint ID        = __umul24( tid_z, __umul24(bdim_x,bdim_y) ) + __umul24( tid_y, bdim_x ) + tid_x;
   const uint dx        = 1;
   const uint dy        = POT_NXT_F;
   const uint dz        = POT_NXT_F*POT_NXT_F;
   const uint DispEven  = ( tid_y + tid_z ) & 1;
   const uint DispOdd   = DispEven^1;
   const uint DispFlip  = bdim_z & 1;
   const uint RhoID0    = __umul24( tid_z, RHO_NXT*RHO_NXT ) + __umul24( tid_y, RHO_NXT )+ ( tid_x << 1 );
   const uint dRhoID    = __umul24( bdim_z, RHO_NXT*RHO_NXT );
#  ifdef SOR_USE_PADDING
   const uint dPotID    = __umul24( bdim_z, POT_NXT_F*POT_NXT_F + POT_PAD*4 );
   const uint warpID    = ID % WARP_SIZE;
   const uint pad_dy_0  = ( warpID >=  8 && warpID <= 15 ) ? dy + POT_PAD : dy;    //
   const uint pad_dy_1  = ( warpID >= 16 && warpID <= 23 ) ? dy + POT_PAD : dy;    // please refer to the Padding notes above!
   const uint pad_dz    = dz + POT_PAD*4;                                          //
   const uint pad_pot   = ( tid_y < 2 ) ? 0 : POT_PAD*((tid_y-2)/4 + 1);
#  else
   const uint dPotID    = __umul24( bdim_z, POT_NXT_F*POT_NXT_F );
   const uint pad_dy_0  = dy;
   const uint pad_dy_1  = dy;
   const uint pad_dz    = dz;
   const uint pad_pot   = 0;
#  endif
   const uint PotID0    = pad_pot + __umul24( 1+tid_z, pad_dz ) + __umul24( 1+tid_y, dy ) + ( tid_x << 1 ) + 1;

   uint ip, im, jp, jm, kp, km, t, s_index;
   uint PotID, RhoID, DispPotID, DispRhoID, Disp;
   real Residual, Residual_Total_Old, Residual_ThreadSum;

   __shared__ real s_Residual_Total;

#  ifdef SOR_USE_PADDING
   __shared__ real s_FPot[ POT_NXT_F*POT_NXT_F*POT_NXT_F + POT_PAD*4*POT_NXT_F ];
#  else
   __shared__ real s_FPot[ POT_NXT_F*POT_NXT_F*POT_NXT_F ];
#  endif

#  ifdef SOR_CPOT_SHARED
   __shared__ real s_CPot[ POT_NXT  *POT_NXT  *POT_NXT   ];
#  endif

#  ifdef SOR_RHO_SHARED
   __shared__ real s_Rho_Array[ RHO_NXT*RHO_NXT*RHO_NXT ];
#  endif


// a1. load the fine-grid density into the shared memory
// -----------------------------------------------------------------------------------------------------------
#  ifdef SOR_RHO_SHARED
   t = ID;
   do {  s_Rho_Array[t] = g_Rho_Array[bid][t];  t += (POT_NTHREAD);}     while ( t < RHO_NXT*RHO_NXT*RHO_NXT);
   __syncthreads();
#  else
   const real *s_Rho_Array = g_Rho_Array[bid];
#  endif


// a2. load the coarse-grid potential into the shared memory
// -----------------------------------------------------------------------------------------------------------
#  ifdef SOR_CPOT_SHARED
   t = ID;
   do {  s_CPot[t] = g_Pot_Array_In[bid][t];    t += POT_NTHREAD; }     while ( t < POT_NXT*POT_NXT*POT_NXT );
   __syncthreads();
#  else
   const real *s_CPot = g_Pot_Array_In[bid];
#  endif



// b. evaluate the "fine-grid" potential by interpolation (as the initial guess and the B.C.)
// -----------------------------------------------------------------------------------------------------------
   const int N_CSlice = POT_NTHREAD / ( (POT_NXT-2)*(POT_NXT-2) );

   if ( ID < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   {
      const real Const_8   = 1.0/8.0;
      const real Const_64  = 1.0/64.0;
      const real Const_512 = 1.0/512.0;

      const int Cdx  = 1;
      const int Cdy  = POT_NXT;
      const int Cdz  = POT_NXT*POT_NXT;
      const int CIDx = 1 + ID % ( POT_NXT-2 );
      const int CIDy = 1 + (  ID % ( (POT_NXT-2)*(POT_NXT-2) )  ) / ( POT_NXT-2 );
      const int CIDz = 1 + ID / ( (POT_NXT-2)*(POT_NXT-2) );
      int       CID  = __mul24( CIDz, Cdz ) + __mul24( CIDy, Cdy ) + __mul24( CIDx, Cdx );
      const int Fdx  = 1;
      const int Fdy  = POT_NXT_F;
      const int FIDx = ( (CIDx-1)<<1 ) - POT_USELESS;
      const int FIDy = ( (CIDy-1)<<1 ) - POT_USELESS;
      int       FIDz = ( (CIDz-1)<<1 ) - POT_USELESS;
#     ifdef SOR_USE_PADDING
      const int Fpad = ( FIDy < 3 ) ? 0 : POT_PAD*((FIDy-3)/4 + 1);    // padding logic
      const int Fdz  = POT_NXT_F*POT_NXT_F + POT_PAD*4;                // added padding
#     else
      const int Fpad = 0;
      const int Fdz  = POT_NXT_F*POT_NXT_F;
#     endif
      int       FID  = Fpad + __mul24( FIDz, Fdz ) + __mul24( FIDy, Fdy ) + __mul24( FIDx, Fdx );

      real TempFPot1, TempFPot2, TempFPot3, TempFPot4, TempFPot5, TempFPot6, TempFPot7, TempFPot8;
      real Slope_00, Slope_01, Slope_02, Slope_03, Slope_04, Slope_05, Slope_06, Slope_07;
      real Slope_08, Slope_09, Slope_10, Slope_11, Slope_12;
      int  Idx, Idy, Idz, ii, jj, kk;


      for (int z=CIDz; z<POT_NXT-1; z+=N_CSlice)
      {
         switch ( IntScheme )
         {
            /*
            case INT_CENTRAL :
            {
               Slope_00 = (real)0.125 * ( s_CPot[CID+Cdx] - s_CPot[CID-Cdx] );
               Slope_01 = (real)0.125 * ( s_CPot[CID+Cdy] - s_CPot[CID-Cdy] );
               Slope_02 = (real)0.125 * ( s_CPot[CID+Cdz] - s_CPot[CID-Cdz] );

               TempFPot1 = s_CPot[CID] - Slope_00 - Slope_01 - Slope_02;
               TempFPot2 = s_CPot[CID] + Slope_00 - Slope_01 - Slope_02;
               TempFPot3 = s_CPot[CID] - Slope_00 + Slope_01 - Slope_02;
               TempFPot4 = s_CPot[CID] + Slope_00 + Slope_01 - Slope_02;
               TempFPot5 = s_CPot[CID] - Slope_00 - Slope_01 + Slope_02;
               TempFPot6 = s_CPot[CID] + Slope_00 - Slope_01 + Slope_02;
               TempFPot7 = s_CPot[CID] - Slope_00 + Slope_01 + Slope_02;
               TempFPot8 = s_CPot[CID] + Slope_00 + Slope_01 + Slope_02;
            }
            break; // INT_CENTRAL
            */


            case INT_CQUAD :
            {
               Slope_00 = Const_8   * ( s_CPot[CID+Cdx        ] - s_CPot[CID-Cdx        ] );
               Slope_01 = Const_8   * ( s_CPot[CID    +Cdy    ] - s_CPot[CID    -Cdy    ] );
               Slope_02 = Const_8   * ( s_CPot[CID        +Cdz] - s_CPot[CID        -Cdz] );

               Slope_03 = Const_64  * ( s_CPot[CID+Cdx    -Cdz] - s_CPot[CID-Cdx    -Cdz] );
               Slope_04 = Const_64  * ( s_CPot[CID    +Cdy-Cdz] - s_CPot[CID    -Cdy-Cdz] );
               Slope_05 = Const_64  * ( s_CPot[CID+Cdx-Cdy    ] - s_CPot[CID-Cdx-Cdy    ] );
               Slope_06 = Const_64  * ( s_CPot[CID+Cdx+Cdy    ] - s_CPot[CID-Cdx+Cdy    ] );
               Slope_07 = Const_64  * ( s_CPot[CID+Cdx    +Cdz] - s_CPot[CID-Cdx    +Cdz] );
               Slope_08 = Const_64  * ( s_CPot[CID    +Cdy+Cdz] - s_CPot[CID    -Cdy+Cdz] );

               Slope_09 = Const_512 * ( s_CPot[CID+Cdx-Cdy-Cdz] - s_CPot[CID-Cdx-Cdy-Cdz] );
               Slope_10 = Const_512 * ( s_CPot[CID+Cdx+Cdy-Cdz] - s_CPot[CID-Cdx+Cdy-Cdz] );
               Slope_11 = Const_512 * ( s_CPot[CID+Cdx-Cdy+Cdz] - s_CPot[CID-Cdx-Cdy+Cdz] );
               Slope_12 = Const_512 * ( s_CPot[CID+Cdx+Cdy+Cdz] - s_CPot[CID-Cdx+Cdy+Cdz] );


               TempFPot1 = - Slope_00 - Slope_01 - Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                           + Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot2 = + Slope_00 - Slope_01 - Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                           - Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot3 = - Slope_00 + Slope_01 - Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                           + Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot4 = + Slope_00 + Slope_01 - Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                           - Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot5 = - Slope_00 - Slope_01 + Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                           - Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];

               TempFPot6 = + Slope_00 - Slope_01 + Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                           + Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot7 = - Slope_00 + Slope_01 + Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                           - Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12 + s_CPot[CID];

               TempFPot8 = + Slope_00 + Slope_01 + Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                           + Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12 + s_CPot[CID];
            }
            break; // INT_CQUAD


            case INT_QUAD :
            {
               TempFPot1 = TempFPot2 = TempFPot3 = TempFPot4 = (real)0.0;
               TempFPot5 = TempFPot6 = TempFPot7 = TempFPot8 = (real)0.0;

               for (int dk=-1; dk<=1; dk++)  {  Idz = dk+1;    kk = __mul24( dk, Cdz );
               for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
               for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

                  TempFPot1 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mm[Idy] * c_Mm[Idx];
                  TempFPot2 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mm[Idy] * c_Mp[Idx];
                  TempFPot3 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mp[Idy] * c_Mm[Idx];
                  TempFPot4 += s_CPot[CID+kk+jj+ii] * c_Mm[Idz] * c_Mp[Idy] * c_Mp[Idx];
                  TempFPot5 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mm[Idy] * c_Mm[Idx];
                  TempFPot6 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mm[Idy] * c_Mp[Idx];
                  TempFPot7 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mp[Idy] * c_Mm[Idx];
                  TempFPot8 += s_CPot[CID+kk+jj+ii] * c_Mp[Idz] * c_Mp[Idy] * c_Mp[Idx];

               }}}
            }
            break; // INT_QUAD

         } // switch ( IntScheme )



//       save data to the shared-memory array.
//       Currently this part is highly diverged. However, since the interpolation takes much less time than the
//       SOR iteration does, we have not yet tried to optimize this part
         if ( FIDz >= 0 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   s_FPot[FID            ] = TempFPot1;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   s_FPot[FID+Fdx        ] = TempFPot2;
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   s_FPot[FID    +Fdy    ] = TempFPot3;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   s_FPot[FID+Fdx+Fdy    ] = TempFPot4;
         }

         if ( FIDz <= POT_NXT_F-2 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   s_FPot[FID        +Fdz] = TempFPot5;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   s_FPot[FID+Fdx    +Fdz] = TempFPot6;
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   s_FPot[FID    +Fdy+Fdz] = TempFPot7;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   s_FPot[FID+Fdx+Fdy+Fdz] = TempFPot8;
         }

         CID  += __mul24(   N_CSlice, Cdz );
         FID  += __mul24( 2*N_CSlice, Fdz );
         FIDz += 2*N_CSlice;

      } // for (int z=CIDz; z<POT_NXT-1; z+=N_CSlice)
   } // if ( ID < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   __syncthreads();



// c. use the SOR scheme to evaluate potential
// -----------------------------------------------------------------------------------------------------------
   Residual_Total_Old = __FLT_MAX__;

   for (uint Iter=0; Iter<Max_Iter; Iter++)
   {
//    (c1). evaluate residual, update potential
//    ==============================================================================
      Residual_ThreadSum = (real)0.0;
      Disp               = DispEven;

      for (uint pass=0; pass<2; pass++)    // pass = (0,1) <--> (even,odd) step
      {
         PotID = PotID0;
         RhoID = RhoID0;

         for (uint z=tid_z; z<RHO_NXT; z+=bdim_z)
         {
            DispPotID = PotID + Disp;
            DispRhoID = RhoID + Disp;

            ip = DispPotID + dx;
            jp = DispPotID + pad_dy_0;
            kp = DispPotID + pad_dz;
            im = DispPotID - dx;
            jm = DispPotID - pad_dy_1;
            km = DispPotID - pad_dz;


//          evaluate the residual
            Residual = (  s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                          - (real)6.0*s_FPot[DispPotID] - Const*s_Rho_Array[DispRhoID]  );

//          update potential
            s_FPot[DispPotID] += Omega_6*Residual;

//          store the sum of the residuals evaluated by the same thread in a per-thread register
            Residual_ThreadSum += FABS( Residual );

            PotID += dPotID;
            RhoID += dRhoID;
            Disp   = Disp^DispFlip;

         } // for (int ZLoop=0; ZLoop<RHO_NXT; ZLoop+=bdim_z)

         Disp = DispOdd;

         __syncthreads();

      } // for (int pass=0; pass<2; pass++)


      if ( Iter+1 >= Min_Iter  &&  Iter % SOR_MOD_REDUCTION == 0 )
      {
//       (c2). perform parallel reduction to get the one-norm of residual
//       ==============================================================================
#        ifdef SOR_USE_SHUFFLE
         Residual_ThreadSum = BlockReduction_Shuffle ( Residual_ThreadSum );
#        else
         Residual_ThreadSum = BlockReduction_WarpSync( Residual_ThreadSum );
#        endif

//       broadcast to all threads
         if ( ID == 0 )    s_Residual_Total = Residual_ThreadSum;
         __syncthreads();


//       (c3). termination criterion
//       ==============================================================================
         if ( s_Residual_Total > Residual_Total_Old )    break;

         Residual_Total_Old = s_Residual_Total;

      } // if ( Iter+1 >= Min_Iter  &&  Iter % SOR_MOD_REDUCTION == 0 )

      __syncthreads();

   } // for (uint Iter=0; Iter<Max_Iter; Iter++)


// d. store potential back to the global memory
// -----------------------------------------------------------------------------------------------------------
   t = ID;

   do
   {
#     ifdef SOR_USE_PADDING
#     define GHOST_DIFF ( POT_GHOST_SIZE - GRA_GHOST_SIZE - 1 )
      uint dy_global  = t % (GRA_NXT*GRA_NXT) / GRA_NXT;
      uint pad_global = ((dy_global + GHOST_DIFF) < 2) ? 0 : POT_PAD*((dy_global + GHOST_DIFF-2)/4 + 1);
#     else
      uint pad_global = 0;
#     endif

      s_index =   __umul24(  t/(GRA_NXT*GRA_NXT)         + POT_GHOST_SIZE - GRA_GHOST_SIZE,  pad_dz  )
                + __umul24(  t%(GRA_NXT*GRA_NXT)/GRA_NXT + POT_GHOST_SIZE - GRA_GHOST_SIZE,  dy  )
                +            t%(GRA_NXT        )         + POT_GHOST_SIZE - GRA_GHOST_SIZE
                             + pad_global;

      g_Pot_Array_Out[bid][t] = s_FPot[s_index];

      t += POT_NTHREAD;
   }
   while ( t < GRA_NXT*GRA_NXT*GRA_NXT );

} // FUNCTION : CUPOT_PoissonSolver_SOR_10to14cube



#endif // #if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == SOR  &&  defined USE_PSOLVER_10TO14 )
