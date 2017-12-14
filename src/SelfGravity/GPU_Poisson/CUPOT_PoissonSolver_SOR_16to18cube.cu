#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == SOR  &&  !defined USE_PSOLVER_10TO14 )



#define POT_NXT_F    ( PATCH_SIZE+2*POT_GHOST_SIZE        )
#define POT_NTHREAD  ( RHO_NXT*RHO_NXT*POT_BLOCK_SIZE_Z/2 )
#define POT_USELESS  ( POT_GHOST_SIZE%2                   )
#define POT_USELESS2 ( POT_USELESS^(GRA_GHOST_SIZE&1)     )
#define POT_NXT_INT  ( (POT_NXT-2)*2                      )
#define ip           ( PotCen + Disp15                    )
#define im           ( PotCen + Disp14                    )
#define jp           ( PotCen + POT_NXT_F/2               )
#define jm           ( PotCen - POT_NXT_F/2               )
#define kp           ( PotCen + dz                        )
#define km           ( PotCen - dz                        )

// additional "__syncthreads" function must be called for POT_GHOST_SIZE == 4 or the emulation mode
#if ( POT_GHOST_SIZE == 4  ||  defined  __DEVICE_EMULATION__ )
   #define SYNCTHREADS()   __syncthreads()
#else
   #define SYNCTHREADS()
#endif

// variables reside in constant memory
#include "CUPOT_PoissonSolver_SetConstMem.cu"




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_PoissonSolver_SOR_16to18cube
// Description :  GPU Poisson solver using the SOR scheme
//
// Note        :  a. Work for POT_GHOST_SIZE = 4, 5 <--> POT_NXT_F = 16, 18
//                b. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                c. We do not use automatic arrays to prevent from using the high-latency local memory
//                   --> unroll loops manually ...
//                d. Reference : Numerical Recipes, Chapter 20.5
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
__global__ void CUPOT_PoissonSolver_SOR_16to18cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ],
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ],
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6,
                                                    const real Const, const IntScheme_t IntScheme )
{

   const uint bx      = blockIdx.x;
   const uint tx      = threadIdx.x;
   const uint ty      = threadIdx.y;
   const uint ID0     = ty*blockDim.x + tx;                    // the 1-D index of thread
   const uint Disp1   = ty&1;                                  // ty = (odd,even) <--> Disp1  = ( 1, 0)
   const uint Disp2   = Disp1^1;                               // ty = (odd,even) <--> Disp2  = ( 0, 1)
   const int  Disp3   = -1 + (int)(Disp2<<1);                  // ty = (odd,even) <--> Disp3  = (-1,+1)
   const int  Disp7   = -Disp1;                                // ty = (odd,even) <--> Disp7  = (-1, 0)
   const int  Disp12  = -Disp3;                                // ty = (odd,even) <--> Disp12 = (+1,-1)
   const int  Disp13  = -Disp2;                                // ty = (odd,even) <--> Disp12 = ( 0,-1)
   const uint dz      = POT_NXT_F*POT_NXT_F/2;
   const uint RhoCen0 = ty*RHO_NXT + (tx<<1);                  // the index of rho
   const uint PotCen0 = dz + __umul24(1+ty, POT_NXT_F/2) + tx; // the index of the left potential
   const uint FloorPow2 = 1<<(31-__clz(POT_NTHREAD) );         // FloorPow2: largest power-of-two value not
   const uint Remain    = POT_NTHREAD - FloorPow2;             //            greater than POT_NTHREAD

   real BPot_xy1, BPot_xy2, BPot_yz1, BPot_yz2, BPot_xz1, BPot_xz2;     // boundary potential stored in registers
   real RPot0, RPot1, RPot2, RPot3, RPot4, RPot5, RPot6, RPot7, RPot8;  // internal potential stored in registers
   real RPot9, RPot10, RPot11, RPot12, RPot13;                          // internal potential stored in registers
#  if ( POT_GHOST_SIZE == 5 )
   real RPot14, RPot15;                                                 // internal potential stored in registers
#  endif
   real Residual, Residual_Total_Old, Temp, Temp1, Temp2;
   uint ID1, ID2, ID3, RhoCen, PotCen, SendID, RecvID, Disp4, Disp5;
   int Disp6, Disp8, Disp9, Disp10, Disp11, Disp14, Disp15, Disp16;
   int ii, jj, Idx, Idy;

   __shared__ real s_FPot  [POT_NXT_F*POT_NXT_F*POT_NXT_F/2];
   __shared__ real s_CPot1 [POT_NXT    *POT_NXT    ];
   __shared__ real s_CPot2 [POT_NXT    *POT_NXT    ];
   __shared__ real s_CPot3 [POT_NXT    *POT_NXT    ];
   __shared__ real s_IntPot[POT_NXT_INT*POT_NXT_INT];
   __shared__ real s_Residual_Total[POT_NTHREAD];

   real *s_CPot_z1, *s_CPot_z2, *s_CPot_z3, *s_Temp;


// a. evaluate the "fine-grid" potential by interpolation (as the initial guess and the B.C.)
// ---------------------------------------------------------------------------------------------------------
   const real Const_8   = 1.0/8.0;
   const real Const_64  = 1.0/64.0;
   const real Const_512 = 1.0/512.0;

   const int Cdx  = 1;
   const int Cdy  = POT_NXT;
   const int CIDx = 1 + ID0 % ( POT_NXT-2 );
   const int CIDy = 1 + ID0 / ( POT_NXT-2 );
   const int CID  = __mul24( CIDy, Cdy ) + __mul24( CIDx, Cdx );

   const int Fdx  = 1;
   const int Fdy  = POT_NXT_INT;
   const int FIDx = (CIDx-1)*2;
   const int FIDy = (CIDy-1)*2;
   const int FID  = __mul24( FIDy, Fdy ) + __mul24( FIDx, Fdx );

   real Slope_00, Slope_01, Slope_02, Slope_03, Slope_04, Slope_05, Slope_06, Slope_07;
   real Slope_08, Slope_09, Slope_10, Slope_11, Slope_12;


// first we load three slices of the coarse-grid potential into the shared memory
   s_CPot_z1 = s_CPot1;
   s_CPot_z2 = s_CPot2;
   s_CPot_z3 = s_CPot3;

   for (uint t=ID0; t<POT_NXT*POT_NXT; t+=POT_NTHREAD)
   {
      ID1          = t + 0*POT_NXT*POT_NXT;
      s_CPot_z1[t] = g_Pot_Array_In[bx][ID1];
   }
   __syncthreads();

   for (uint t=ID0; t<POT_NXT*POT_NXT; t+=POT_NTHREAD)
   {
      ID1          = t + 1*POT_NXT*POT_NXT;
      s_CPot_z2[t] = g_Pot_Array_In[bx][ID1];
   }
   __syncthreads();

   for (uint t=ID0; t<POT_NXT*POT_NXT; t+=POT_NTHREAD)
   {
      ID1          = t + 2*POT_NXT*POT_NXT;
      s_CPot_z3[t] = g_Pot_Array_In[bx][ID1];
   }
   __syncthreads();


// (a1). interpolation : the lowest z plane
// ===========================================================================
   if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
   {
      switch ( IntScheme )
      {
         /*
         case INT_CENTRAL :
         {
            Slope_00 = (real)0.125 * ( s_CPot_z2[CID+Cdx] - s_CPot_z2[CID-Cdx] );
            Slope_01 = (real)0.125 * ( s_CPot_z2[CID+Cdy] - s_CPot_z2[CID-Cdy] );
            Slope_02 = (real)0.125 * ( s_CPot_z3[CID    ] - s_CPot_z1[CID    ] );

#if ( POT_GHOST_SIZE == 4 )   // lower plane

            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 - Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 - Slope_02;

#else                         // upper plane

            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 + Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 + Slope_02;

#endif
         }
         break; // INT_CENTRAL
         */


         case INT_CQUAD :
         {
            Slope_00 = Const_8   * ( s_CPot_z2[CID+Cdx    ] - s_CPot_z2[CID-Cdx    ] );
            Slope_01 = Const_8   * ( s_CPot_z2[CID+Cdy    ] - s_CPot_z2[CID-Cdy    ] );
            Slope_02 = Const_8   * ( s_CPot_z3[CID        ] - s_CPot_z1[CID        ] );

            Slope_03 = Const_64  * ( s_CPot_z1[CID+Cdx    ] - s_CPot_z1[CID-Cdx    ] );
            Slope_04 = Const_64  * ( s_CPot_z1[CID    +Cdy] - s_CPot_z1[CID    -Cdy] );
            Slope_05 = Const_64  * ( s_CPot_z2[CID+Cdx-Cdy] - s_CPot_z2[CID-Cdx-Cdy] );
            Slope_06 = Const_64  * ( s_CPot_z2[CID+Cdx+Cdy] - s_CPot_z2[CID-Cdx+Cdy] );
            Slope_07 = Const_64  * ( s_CPot_z3[CID+Cdx    ] - s_CPot_z3[CID-Cdx    ] );
            Slope_08 = Const_64  * ( s_CPot_z3[CID    +Cdy] - s_CPot_z3[CID    -Cdy] );

            Slope_09 = Const_512 * ( s_CPot_z1[CID+Cdx-Cdy] - s_CPot_z1[CID-Cdx-Cdy] );
            Slope_10 = Const_512 * ( s_CPot_z1[CID+Cdx+Cdy] - s_CPot_z1[CID-Cdx+Cdy] );
            Slope_11 = Const_512 * ( s_CPot_z3[CID+Cdx-Cdy] - s_CPot_z3[CID-Cdx-Cdy] );
            Slope_12 = Const_512 * ( s_CPot_z3[CID+Cdx+Cdy] - s_CPot_z3[CID-Cdx+Cdy] );


#if ( POT_GHOST_SIZE == 4 )   // lower plane

            s_IntPot[FID        ] = - Slope_00 - Slope_01 - Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 - Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 - Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 - Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

#else                         // upper plane

            s_IntPot[FID        ] = - Slope_00 - Slope_01 + Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 + Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 + Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 + Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

#endif
         }
         break; // INT_CQUAD


         case INT_QUAD :
         {
            s_IntPot[FID        ] = (real)0.0;
            s_IntPot[FID+Fdx    ] = (real)0.0;
            s_IntPot[FID    +Fdy] = (real)0.0;
            s_IntPot[FID+Fdx+Fdy] = (real)0.0;

            for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
            for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

#if ( POT_GHOST_SIZE == 4 )   // lower plane

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mp[Idx];

#else                         // upper plane

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mp[Idx];

#endif
            }} // for di,dj
         }
         break; // INT_QUAD

      } // switch ( IntScheme )
   } // if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

   __syncthreads();


// store data into shared memory
   ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+Disp1;
   ID2         = __umul24( 1            +ty, POT_NXT_F/2 ) + tx                   +Disp1;
   s_FPot[ID2] = s_IntPot[ID1];


// store data into registers
   ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+Disp2;
   BPot_xy1    = s_IntPot[ID1];

   __syncthreads();


// for POT_USELESS == 0, no z plane is useless --> one more z plane (upper plane) to store
#if ( POT_GHOST_SIZE == 4 )

   if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
   {
      switch ( IntScheme )
      {
         /*
         case INT_CENTRAL :
         {
            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 + Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 + Slope_02;
         }
         break; // INT_CENTRAL
         */


         case INT_CQUAD :
         {
            s_IntPot[FID        ] = - Slope_00 - Slope_01 + Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 + Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 + Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 + Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];
         }
         break; // INT_CQUAD


         case INT_QUAD :
         {
            s_IntPot[FID        ] = (real)0.0;
            s_IntPot[FID+Fdx    ] = (real)0.0;
            s_IntPot[FID    +Fdy] = (real)0.0;
            s_IntPot[FID+Fdx+Fdy] = (real)0.0;

            for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
            for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mp[Idx];

            }} // for di,dj
         }
         break; // INT_QUAD

      } // switch ( IntScheme )
   } // if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

   __syncthreads();


// store the internal potential into shared memory
   ID1         = __umul24( 1+ty, POT_NXT_INT ) + 1+(tx<<1)+Disp2;
   ID2         = 1*POT_NXT_F*POT_NXT_F/2 + __umul24( 1+ty, POT_NXT_F/2) + tx+Disp2;
   s_FPot[ID2] = s_IntPot[ID1];


// store the internal potential into registers
   ID1         = __umul24( 1+ty, POT_NXT_INT ) + 1+(tx<<1)+Disp1;
   RPot0       = s_IntPot[ID1];


// store the boundary potential into shared memory or registers
   if ( ID0 < RHO_NXT/2 )
   {
      ID3         = ID0;

//    shared memory: -yz plane
      ID1         = __umul24( 1+(ID3<<1), POT_NXT_INT );
      ID2         = POT_NXT_F*POT_NXT_F/2 + __umul24( 1+(ID3<<1), POT_NXT_F/2 );
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: +yz plane
      ID1         = __umul24( 2+(ID3<<1), POT_NXT_INT ) + POT_NXT_INT-1;
      ID2         = POT_NXT_F*POT_NXT_F/2 + __umul24( 2+(ID3<<1), POT_NXT_F/2 ) + POT_NXT_F/2-1;
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: -xz plane
      ID1         = 1+(ID3<<1);
      ID2         = POT_NXT_F*POT_NXT_F/2 + ID3;
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: +xz plane
      ID1         = (POT_NXT_INT-1)*POT_NXT_INT + 2+(ID3<<1);
      ID2         = POT_NXT_F*POT_NXT_F/2 + (POT_NXT_F-1)*POT_NXT_F/2 + ID3+1;
      s_FPot[ID2] = s_IntPot[ID1];


//    registers: -yz plane
      ID1         = __umul24( 2+(ID3<<1), POT_NXT_INT );
      BPot_yz1    = s_IntPot[ID1];

//    registers: +yz plane
      ID1         = __umul24( 1+(ID3<<1), POT_NXT_INT ) + POT_NXT_INT-1;
      BPot_yz2    = s_IntPot[ID1];

//    registers: -xz plane
      ID1         = 2+(ID3<<1);
      BPot_xz1    = s_IntPot[ID1];

//    registers: +xz plane
      ID1         = (POT_NXT_INT-1)*POT_NXT_INT + 1+(ID3<<1);
      BPot_xz2    = s_IntPot[ID1];

   } // if ( ID0 < RHO_NXT/2 )

   __syncthreads();

#endif // ( POT_GHOST_SIZE == 4 )


// (a2). interpolation : central z planes
// ===========================================================================
   for (uint Cz=0; Cz<=POT_NXT-5; Cz++)
   {
      s_Temp    = s_CPot_z1;
      s_CPot_z1 = s_CPot_z2;
      s_CPot_z2 = s_CPot_z3;
      s_CPot_z3 = s_Temp;


//    load one slice of the coarse-grid potential into the shared memory
      for (uint t=ID0; t<POT_NXT*POT_NXT; t+=POT_NTHREAD)
      {
         ID1          = t + __umul24( Cz+3, POT_NXT*POT_NXT );
         s_CPot_z3[t] = g_Pot_Array_In[bx][ID1];
      }
      __syncthreads();


//    interpolation
      if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
      {
         switch ( IntScheme )
         {
            /*
            case INT_CENTRAL :
            {
               Slope_00 = (real)0.125 * ( s_CPot_z2[CID+Cdx] - s_CPot_z2[CID-Cdx] );
               Slope_01 = (real)0.125 * ( s_CPot_z2[CID+Cdy] - s_CPot_z2[CID-Cdy] );
               Slope_02 = (real)0.125 * ( s_CPot_z3[CID    ] - s_CPot_z1[CID    ] );
            }
            break; // INT_CENTRAL
            */


            case INT_CQUAD :
            {
               Slope_00 = Const_8   * ( s_CPot_z2[CID+Cdx    ] - s_CPot_z2[CID-Cdx    ] );
               Slope_01 = Const_8   * ( s_CPot_z2[CID+Cdy    ] - s_CPot_z2[CID-Cdy    ] );
               Slope_02 = Const_8   * ( s_CPot_z3[CID        ] - s_CPot_z1[CID        ] );

               Slope_03 = Const_64  * ( s_CPot_z1[CID+Cdx    ] - s_CPot_z1[CID-Cdx    ] );
               Slope_04 = Const_64  * ( s_CPot_z1[CID    +Cdy] - s_CPot_z1[CID    -Cdy] );
               Slope_05 = Const_64  * ( s_CPot_z2[CID+Cdx-Cdy] - s_CPot_z2[CID-Cdx-Cdy] );
               Slope_06 = Const_64  * ( s_CPot_z2[CID+Cdx+Cdy] - s_CPot_z2[CID-Cdx+Cdy] );
               Slope_07 = Const_64  * ( s_CPot_z3[CID+Cdx    ] - s_CPot_z3[CID-Cdx    ] );
               Slope_08 = Const_64  * ( s_CPot_z3[CID    +Cdy] - s_CPot_z3[CID    -Cdy] );

               Slope_09 = Const_512 * ( s_CPot_z1[CID+Cdx-Cdy] - s_CPot_z1[CID-Cdx-Cdy] );
               Slope_10 = Const_512 * ( s_CPot_z1[CID+Cdx+Cdy] - s_CPot_z1[CID-Cdx+Cdy] );
               Slope_11 = Const_512 * ( s_CPot_z3[CID+Cdx-Cdy] - s_CPot_z3[CID-Cdx-Cdy] );
               Slope_12 = Const_512 * ( s_CPot_z3[CID+Cdx+Cdy] - s_CPot_z3[CID-Cdx+Cdy] );
            }
            break; // INT_CQUAD

         } // switch ( IntScheme )
      } // if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

#if ( POT_GHOST_SIZE == 4 )
      Disp16 = Disp2;
#else
      Disp16 = Disp1;
#endif


//    since the amount of shared memory is exhausted, we can only save one z plane at a time
      for (int UpDown=0; UpDown<2; UpDown++)
      {
         const real Sign = (real)-1.0 + (real)2.0*(real)UpDown;

         if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
         {
            switch ( IntScheme )
            {
               /*
               case INT_CENTRAL :
               {
                  s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 + Sign*Slope_02;
                  s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 + Sign*Slope_02;
                  s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 + Sign*Slope_02;
                  s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 + Sign*Slope_02;
               }
               break; // INT_CENTRAL
               */


               case INT_CQUAD :
               {
                  s_IntPot[FID        ] = - Slope_00 - Slope_01 - Slope_05 + Slope_06
                                          - Sign*( - Slope_02 - Slope_03 - Slope_04 + Slope_07 + Slope_08
                                          - Slope_09 + Slope_10 + Slope_11 - Slope_12 )
                                          + s_CPot_z2[CID];

                  s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 + Slope_05 - Slope_06
                                          - Sign*( - Slope_02 + Slope_03 - Slope_04 - Slope_07 + Slope_08
                                          + Slope_09 - Slope_10 - Slope_11 + Slope_12 )
                                          + s_CPot_z2[CID];

                  s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 + Slope_05 - Slope_06
                                          - Sign*( - Slope_02 - Slope_03 + Slope_04 + Slope_07 - Slope_08
                                          + Slope_09 - Slope_10 - Slope_11 + Slope_12 )
                                          + s_CPot_z2[CID];

                  s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 - Slope_05 + Slope_06
                                          - Sign*( - Slope_02 + Slope_03 + Slope_04 - Slope_07 - Slope_08
                                          - Slope_09 + Slope_10 + Slope_11 - Slope_12 )
                                          + s_CPot_z2[CID];
               }
               break; // INT_CQUAD


               case INT_QUAD :
               {
                  s_IntPot[FID        ] = (real)0.0;
                  s_IntPot[FID+Fdx    ] = (real)0.0;
                  s_IntPot[FID    +Fdy] = (real)0.0;
                  s_IntPot[FID+Fdx+Fdy] = (real)0.0;

                  if ( UpDown == 0 )
                  {
                     for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
                     for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

                        s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mp[Idx];

                     }} // for di,dj
                  } // if ( UpDown == 0 )

                  else
                  {
                     for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
                     for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

                        s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mm[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mp[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mm[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mp[Idx];
                        s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mp[Idx];

                     }} // for di,dj
                  } // if ( UpDown == 0 ) ... else ...
               }
               break; // INT_QUAD

            } // switch ( IntScheme )
         } //if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

         __syncthreads();


//       store the internal potential into shared memory
         ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+(Disp16^1);
         ID2         = __umul24( 2-POT_USELESS+UpDown+Cz*2, POT_NXT_F*POT_NXT_F/2 )
                       + __umul24( 1+ty, POT_NXT_F/2) + tx+(Disp16^1);
         s_FPot[ID2] = s_IntPot[ID1];


//       store the internal potential into registers
         ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+Disp16;

#if ( POT_GHOST_SIZE == 4 )

         if ( UpDown == 0 )
         {
            switch ( Cz )
            {
               case 0:  RPot1  = s_IntPot[ID1];  break;
               case 1:  RPot3  = s_IntPot[ID1];  break;
               case 2:  RPot5  = s_IntPot[ID1];  break;
               case 3:  RPot7  = s_IntPot[ID1];  break;
               case 4:  RPot9  = s_IntPot[ID1];  break;
               case 5:  RPot11 = s_IntPot[ID1];  break;
            }
         }

         else
         {
            switch ( Cz )
            {
               case 0:  RPot2  = s_IntPot[ID1];  break;
               case 1:  RPot4  = s_IntPot[ID1];  break;
               case 2:  RPot6  = s_IntPot[ID1];  break;
               case 3:  RPot8  = s_IntPot[ID1];  break;
               case 4:  RPot10 = s_IntPot[ID1];  break;
               case 5:  RPot12 = s_IntPot[ID1];  break;
            }
         }

#else

         if ( UpDown == 0 )
         {
            switch ( Cz )
            {
               case 0:  RPot0  = s_IntPot[ID1];  break;
               case 1:  RPot2  = s_IntPot[ID1];  break;
               case 2:  RPot4  = s_IntPot[ID1];  break;
               case 3:  RPot6  = s_IntPot[ID1];  break;
               case 4:  RPot8  = s_IntPot[ID1];  break;
               case 5:  RPot10 = s_IntPot[ID1];  break;
               case 6:  RPot12 = s_IntPot[ID1];  break;
               case 7:  RPot14 = s_IntPot[ID1];  break;
            }
         }

         else
         {
            switch ( Cz )
            {
               case 0:  RPot1  = s_IntPot[ID1];  break;
               case 1:  RPot3  = s_IntPot[ID1];  break;
               case 2:  RPot5  = s_IntPot[ID1];  break;
               case 3:  RPot7  = s_IntPot[ID1];  break;
               case 4:  RPot9  = s_IntPot[ID1];  break;
               case 5:  RPot11 = s_IntPot[ID1];  break;
               case 6:  RPot13 = s_IntPot[ID1];  break;
               case 7:  RPot15 = s_IntPot[ID1];  break;
            }
         }

#endif // #if ( POT_GHOST_SIZE == 4 )


//       store the boundary potential into shared memory or registers
         if (  ID0 >= (2*Cz+UpDown+1-POT_USELESS)*RHO_NXT/2  &&  ID0 < (2*Cz+UpDown+2-POT_USELESS)*RHO_NXT/2  )
         {
            ID3 = ID0 % (RHO_NXT/2) ;

//          shared memory: -yz plane
            ID1         = __umul24( 2-UpDown+2*POT_USELESS*UpDown+(ID3<<1), POT_NXT_INT ) + POT_USELESS;
            ID2         = __umul24( 2-POT_USELESS+UpDown+Cz*2, POT_NXT_F*POT_NXT_F/2 )
                          + __umul24( 2-UpDown+2*POT_USELESS*UpDown-POT_USELESS+(ID3<<1), POT_NXT_F/2 );
            s_FPot[ID2] = s_IntPot[ID1];

//          shared memory: +yz plane
            ID1         = __umul24( 1+UpDown-2*POT_USELESS*UpDown+2*POT_USELESS+(ID3<<1), POT_NXT_INT )
                          + POT_NXT_INT-1-POT_USELESS;
            ID2         = __umul24( 2-POT_USELESS+UpDown+Cz*2, POT_NXT_F*POT_NXT_F/2 )
                          + __umul24( 1+UpDown-2*POT_USELESS*UpDown+POT_USELESS+(ID3<<1), POT_NXT_F/2 )
                          + POT_NXT_F/2-1;
            s_FPot[ID2] = s_IntPot[ID1];

//          shared memory: -xz plane
            ID1         = POT_USELESS*POT_NXT_INT + 2-UpDown+2*POT_USELESS*UpDown+(ID3<<1);
            ID2         = __umul24( 2-POT_USELESS+UpDown+Cz*2, POT_NXT_F*POT_NXT_F/2 )
                          + 1-UpDown-POT_USELESS+2*POT_USELESS*UpDown+ID3;
            s_FPot[ID2] = s_IntPot[ID1];

//          shared memory: +xz plane
            ID1         = (POT_NXT_INT-1-POT_USELESS)*POT_NXT_INT
                           + 1+UpDown-2*POT_USELESS*UpDown+2*POT_USELESS+(ID3<<1);
            ID2         = __umul24( 2-POT_USELESS+UpDown+Cz*2, POT_NXT_F*POT_NXT_F/2 )
                          + (POT_NXT_F-1)*POT_NXT_F/2 + UpDown-2*POT_USELESS*UpDown+POT_USELESS+ID3;
            s_FPot[ID2] = s_IntPot[ID1];


//          registers: -yz plane
            ID1         = __umul24( 1+UpDown-2*POT_USELESS*UpDown+2*POT_USELESS+(ID3<<1), POT_NXT_INT )
                          + POT_USELESS;
            BPot_yz1    = s_IntPot[ID1];

//          registers: +yz plane
            ID1         = __umul24( 2-UpDown+2*POT_USELESS*UpDown+(ID3<<1), POT_NXT_INT )
                          + POT_NXT_INT-1-POT_USELESS;
            BPot_yz2    = s_IntPot[ID1];

//          registers: -xz plane
            ID1         = POT_USELESS*POT_NXT_INT + 1+UpDown-2*POT_USELESS*UpDown+2*POT_USELESS+(ID3<<1);
            BPot_xz1    = s_IntPot[ID1];

//          registers: +xz plane
            ID1         = (POT_NXT_INT-1-POT_USELESS)*POT_NXT_INT + 2-UpDown+2*POT_USELESS*UpDown+(ID3<<1);
            BPot_xz2    = s_IntPot[ID1];

         } // if (  ID0 >= (2*Cz+UpDown+1-POT_USELESS)*RHO_NXT/2 && ID0 < (2*Cz+UpDown+2-POT_USELESS)*RHO_NXT/2  )

         Disp16 = Disp16^1;

         __syncthreads();

      } //  for (int UpDown=0; UpDown<2; UpDown++)

   } // for (uint Cz=0; Cz<=POT_NXT-5; Cz++)


// (a3). interpolation : the highest z plane
// ===========================================================================

// load one slice of the coarse-grid potential into shared memory
   s_Temp    = s_CPot_z1;
   s_CPot_z1 = s_CPot_z2;
   s_CPot_z2 = s_CPot_z3;
   s_CPot_z3 = s_Temp;

   for (uint t=ID0; t<POT_NXT*POT_NXT; t+=POT_NTHREAD)
   {
      ID1          = t + (POT_NXT-1)*POT_NXT*POT_NXT;
      s_CPot_z3[t] = g_Pot_Array_In[bx][ID1];
   }
   __syncthreads();


// interpolation
   if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
   {
      switch ( IntScheme )
      {
         /*
         case INT_CENTRAL :
         {
            Slope_00 = (real)0.125 * ( s_CPot_z2[CID+Cdx] - s_CPot_z2[CID-Cdx] );
            Slope_01 = (real)0.125 * ( s_CPot_z2[CID+Cdy] - s_CPot_z2[CID-Cdy] );
            Slope_02 = (real)0.125 * ( s_CPot_z3[CID    ] - s_CPot_z1[CID    ] );

#if ( POT_GHOST_SIZE == 4 )   // upper plane

            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 + Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 + Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 + Slope_02;

#else                         // lower plane

            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 - Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 - Slope_02;

#endif
         }
         break; // INT_CENTRAL
         */


         case INT_CQUAD :
         {
            Slope_00 = Const_8   * ( s_CPot_z2[CID+Cdx    ] - s_CPot_z2[CID-Cdx    ] );
            Slope_01 = Const_8   * ( s_CPot_z2[CID+Cdy    ] - s_CPot_z2[CID-Cdy    ] );
            Slope_02 = Const_8   * ( s_CPot_z3[CID        ] - s_CPot_z1[CID        ] );

            Slope_03 = Const_64  * ( s_CPot_z1[CID+Cdx    ] - s_CPot_z1[CID-Cdx    ] );
            Slope_04 = Const_64  * ( s_CPot_z1[CID    +Cdy] - s_CPot_z1[CID    -Cdy] );
            Slope_05 = Const_64  * ( s_CPot_z2[CID+Cdx-Cdy] - s_CPot_z2[CID-Cdx-Cdy] );
            Slope_06 = Const_64  * ( s_CPot_z2[CID+Cdx+Cdy] - s_CPot_z2[CID-Cdx+Cdy] );
            Slope_07 = Const_64  * ( s_CPot_z3[CID+Cdx    ] - s_CPot_z3[CID-Cdx    ] );
            Slope_08 = Const_64  * ( s_CPot_z3[CID    +Cdy] - s_CPot_z3[CID    -Cdy] );

            Slope_09 = Const_512 * ( s_CPot_z1[CID+Cdx-Cdy] - s_CPot_z1[CID-Cdx-Cdy] );
            Slope_10 = Const_512 * ( s_CPot_z1[CID+Cdx+Cdy] - s_CPot_z1[CID-Cdx+Cdy] );
            Slope_11 = Const_512 * ( s_CPot_z3[CID+Cdx-Cdy] - s_CPot_z3[CID-Cdx-Cdy] );
            Slope_12 = Const_512 * ( s_CPot_z3[CID+Cdx+Cdy] - s_CPot_z3[CID-Cdx+Cdy] );


#if ( POT_GHOST_SIZE == 4 )   // upper plane

            s_IntPot[FID        ] = - Slope_00 - Slope_01 + Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 + Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 + Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 + Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

#else                         // lower plane

            s_IntPot[FID        ] = - Slope_00 - Slope_01 - Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 - Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 - Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 - Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

#endif
         }
         break; // INT_CQUAD


         case INT_QUAD :
         {
            s_IntPot[FID        ] = (real)0.0;
            s_IntPot[FID+Fdx    ] = (real)0.0;
            s_IntPot[FID    +Fdy] = (real)0.0;
            s_IntPot[FID+Fdx+Fdy] = (real)0.0;

            for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
            for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

#if ( POT_GHOST_SIZE == 4 )   // upper plane

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mp[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mp[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mp[2] * Mp[Idy] * Mp[Idx];

#else                         // lower plane

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mp[Idx];

#endif
            }} // for di,dj
         }
         break; // INT_QUAD

      } // switch ( IntScheme )
   } // if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

   __syncthreads();


// store the boundary potential into shared memory
   ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+Disp2;
   ID2         = (POT_NXT_F-1)*POT_NXT_F*POT_NXT_F/2 + __umul24( 1+ty, POT_NXT_F/2 ) + tx+Disp2;
   s_FPot[ID2] = s_IntPot[ID1];


// store the boundary potential into registers
   ID1         = __umul24( 1+POT_USELESS+ty, POT_NXT_INT ) + 1+POT_USELESS+(tx<<1)+Disp1;
   BPot_xy2    = s_IntPot[ID1];

   __syncthreads();


// for POT_USELESS == 0, no z plane is useless --> one more z plane (lower plane) to store
#if ( POT_GHOST_SIZE == 4 )

   if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )
   {
      switch ( IntScheme )
      {
         /*
         case INT_CENTRAL :
         {
            s_IntPot[FID        ] = s_CPot_z2[CID] - Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID+Fdx    ] = s_CPot_z2[CID] + Slope_00 - Slope_01 - Slope_02;
            s_IntPot[FID    +Fdy] = s_CPot_z2[CID] - Slope_00 + Slope_01 - Slope_02;
            s_IntPot[FID+Fdx+Fdy] = s_CPot_z2[CID] + Slope_00 + Slope_01 - Slope_02;
         }
         break; // INT_CENTRAL
         */


         case INT_CQUAD :
         {
            s_IntPot[FID        ] = - Slope_00 - Slope_01 - Slope_02 - Slope_03 - Slope_04 - Slope_05 + Slope_06
                                    + Slope_07 + Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx    ] = + Slope_00 - Slope_01 - Slope_02 + Slope_03 - Slope_04 + Slope_05 - Slope_06
                                    - Slope_07 + Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID    +Fdy] = - Slope_00 + Slope_01 - Slope_02 - Slope_03 + Slope_04 + Slope_05 - Slope_06
                                    + Slope_07 - Slope_08 + Slope_09 - Slope_10 - Slope_11 + Slope_12
                                    + s_CPot_z2[CID];

            s_IntPot[FID+Fdx+Fdy] = + Slope_00 + Slope_01 - Slope_02 + Slope_03 + Slope_04 - Slope_05 + Slope_06
                                    - Slope_07 - Slope_08 - Slope_09 + Slope_10 + Slope_11 - Slope_12
                                    + s_CPot_z2[CID];
         }
         break; // INT_CQUAD


         case INT_QUAD :
         {
            s_IntPot[FID        ] = (real)0.0;
            s_IntPot[FID+Fdx    ] = (real)0.0;
            s_IntPot[FID    +Fdy] = (real)0.0;
            s_IntPot[FID+Fdx+Fdy] = (real)0.0;

            for (int dj=-1; dj<=1; dj++)  {  Idy = dj+1;    jj = __mul24( dj, Cdy );
            for (int di=-1; di<=1; di++)  {  Idx = di+1;    ii = __mul24( di, Cdx );

               s_IntPot[FID        ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID        ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx    ] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mm[Idy] * Mp[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID    +Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mm[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z1[CID+jj+ii] * Mm[0] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z2[CID+jj+ii] * Mm[1] * Mp[Idy] * Mp[Idx];
               s_IntPot[FID+Fdx+Fdy] += s_CPot_z3[CID+jj+ii] * Mm[2] * Mp[Idy] * Mp[Idx];

            }} // for di,dj
         }
         break; // INT_QUAD

      } // switch ( IntScheme )
   } // if ( ID0 < (POT_NXT-2)*(POT_NXT-2) )

   __syncthreads();


// store the internal potential into shared memory
   ID1         = __umul24( 1+ty, POT_NXT_INT ) + 1+(tx<<1)+Disp1;
   ID2         = (POT_NXT_F-2)*POT_NXT_F*POT_NXT_F/2 + __umul24( 1+ty, POT_NXT_F/2) + tx+Disp1;
   s_FPot[ID2] = s_IntPot[ID1];


// store the internal potential into registers
   ID1         = __umul24( 1+ty, POT_NXT_INT ) + 1+(tx<<1)+Disp2;
   RPot13      = s_IntPot[ID1];

// store the boundary potential into shared memory or registers
   if (  ID0 >= 13*RHO_NXT/2  &&  ID0 < 14*RHO_NXT/2  )
   {
      ID3         = ID0 % (RHO_NXT/2);

//    shared memory: -yz plane
      ID1         = __umul24( 2+(ID3<<1), POT_NXT_INT );
      ID2         = (POT_NXT_F-2)*POT_NXT_F*POT_NXT_F/2 + __umul24( 2+(ID3<<1), POT_NXT_F/2 );
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: +yz plane
      ID1         = __umul24( 1+(ID3<<1), POT_NXT_INT ) + POT_NXT_INT-1;
      ID2         = (POT_NXT_F-2)*POT_NXT_F*POT_NXT_F/2 + __umul24( 1+(ID3<<1), POT_NXT_F/2 ) + POT_NXT_F/2-1;
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: -xz plane
      ID1         = 2+(ID3<<1);
      ID2         = (POT_NXT_F-2)*POT_NXT_F*POT_NXT_F/2 + ID3+1;
      s_FPot[ID2] = s_IntPot[ID1];

//    shared memory: +xz plane
      ID1         = (POT_NXT_INT-1)*POT_NXT_INT + 1+(ID3<<1);
      ID2         = (POT_NXT_F-2)*POT_NXT_F*POT_NXT_F/2 + (POT_NXT_F-1)*POT_NXT_F/2 + ID3;
      s_FPot[ID2] = s_IntPot[ID1];


//    registers: -yz plane
      ID1         = __umul24( 1+(ID3<<1), POT_NXT_INT );
      BPot_yz1    = s_IntPot[ID1];

//    registers: +yz plane
      ID1         = __umul24( 2+(ID3<<1), POT_NXT_INT ) + POT_NXT_INT-1;
      BPot_yz2    = s_IntPot[ID1];

//    registers: -xz plane
      ID1         = 1+(ID3<<1);
      BPot_xz1    = s_IntPot[ID1];

//    registers: +xz plane
      ID1         = (POT_NXT_INT-1)*POT_NXT_INT + 2+(ID3<<1);
      BPot_xz2    = s_IntPot[ID1];

   } // if (  ID0 >= 13*RHO_NXT/2  &&  ID0 < 14*RHO_NXT/2  )

   __syncthreads();

#endif // ( POT_GHOST_SIZE == 4 )



// b. use the SOR scheme to evaluate potential
// ---------------------------------------------------------------------------------------------------------
   Residual_Total_Old = __FLT_MAX__;

   for (int Iter=0; Iter<Max_Iter; Iter++)
   {
      s_Residual_Total[ID0] = (real)0.0;

      PotCen   = PotCen0 + Disp1;
      RhoCen   = RhoCen0 + Disp1;
      Disp4    = Disp1;
      Disp5    = Disp2;
      Disp6    = Disp3;
      Disp8    = Disp7;
      Disp9    = Disp2;
      Disp10   = Disp13;
      Disp11   = Disp1;
      Disp14   = Disp8;
      Disp15   = Disp9;

      for (int pass=0; pass<2; pass++)
      {

//       (b1). evaluate residual, update potential
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//       z = 0 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot0 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot0 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 1 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot1 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot1 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 2 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot2 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot2 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 3 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot3 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot3 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 4 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot4 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot4 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 5 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot5 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot5 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 6 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot6 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot6 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 7 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot7 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot7 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 8 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot8 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot8 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 9 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot9 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot9 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 10 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot10 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot10 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 11 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot11 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot11 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 12 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot12 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot12 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 13 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot13 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot13 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

#if ( POT_GHOST_SIZE == 5 )

         PotCen += dz - Disp6;
         RhoCen += RHO_NXT*RHO_NXT - Disp6;
         Disp14  = Disp8;
         Disp15  = Disp9;


//       z = 14 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot14 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot14 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

         PotCen += dz + Disp6;
         RhoCen += RHO_NXT*RHO_NXT + Disp6;
         Disp14  = Disp10;
         Disp15  = Disp11;


//       z = 15 plane
// =======================================
//       evaluate residual
         Residual = ( s_FPot[kp] + s_FPot[km] + s_FPot[jp] + s_FPot[jm] + s_FPot[ip] + s_FPot[im]
                      - (real)6.0*RPot15 - Const*g_Rho_Array[bx][RhoCen] );

//       update potential
         RPot15 += Omega_6*Residual;

//       save the absolute value of residual of each grid into a shared array for evaluating the sum
         s_Residual_Total[ID0] += FABS( Residual );

#endif // #if ( POT_GHOST_SIZE == 5 )


         __syncthreads();



//       (b2). exchange the potential stored in the shared memory and registers
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//       (b2-1). exchange the boundary potential
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//       -xy plane
// =======================================
         SendID         = __umul24(1+ty, POT_NXT_F/2) + tx+Disp5;
         RecvID         = __umul24(1+ty, POT_NXT_F/2) + tx+Disp4;
         Temp           = BPot_xy1;
         BPot_xy1       = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       +xy plane
// =======================================
         SendID         = (POT_NXT_F-1)*POT_NXT_F*POT_NXT_F/2 + __umul24(1+ty, POT_NXT_F/2) + tx+Disp4;
         RecvID         = (POT_NXT_F-1)*POT_NXT_F*POT_NXT_F/2 + __umul24(1+ty, POT_NXT_F/2) + tx+Disp5;
         Temp           = BPot_xy2;
         BPot_xy2       = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       -yz plane (boundary potential in the +-yz planes will be pasted after the exchange of internal potential)
// =======================================
         RecvID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + __umul24(1+(tx<<1)+Disp4, POT_NXT_F/2);
         Temp1          = BPot_yz1;
         BPot_yz1       = s_FPot[RecvID];

//       +yz plane
// =======================================
         RecvID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + __umul24(1+(tx<<1)+Disp5, POT_NXT_F/2)
                          + (POT_NXT_F/2)-1;
         Temp2          = BPot_yz2;
         BPot_yz2       = s_FPot[RecvID];

//       -xz plane
// =======================================
         SendID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + tx+Disp5;
         RecvID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + tx+Disp4;
         Temp           = BPot_xz1;
         BPot_xz1       = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       +xz plane
// =======================================
         SendID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + (POT_NXT_F-1)*POT_NXT_F/2 + tx+Disp4;
         RecvID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + (POT_NXT_F-1)*POT_NXT_F/2 + tx+Disp5;
         Temp           = BPot_xz2;
         BPot_xz2       = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;


         __syncthreads();



//       (b2-2). exchange the internal potential stored in shared memory and registers
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//       z = 0 plane
// =======================================
         SendID         = PotCen0 +  0*dz + Disp4;
         RecvID         = PotCen0 +  0*dz + Disp5;
         Temp           = RPot0;
         RPot0          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 1 plane
// =======================================
         SendID         = PotCen0 +  1*dz + Disp5;
         RecvID         = PotCen0 +  1*dz + Disp4;
         Temp           = RPot1;
         RPot1          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 2 plane
// =======================================
         SendID         = PotCen0 +  2*dz + Disp4;
         RecvID         = PotCen0 +  2*dz + Disp5;
         Temp           = RPot2;
         RPot2          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 3 plane
// =======================================
         SendID         = PotCen0 +  3*dz + Disp5;
         RecvID         = PotCen0 +  3*dz + Disp4;
         Temp           = RPot3;
         RPot3          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 4 plane
// =======================================
         SendID         = PotCen0 +  4*dz + Disp4;
         RecvID         = PotCen0 +  4*dz + Disp5;
         Temp           = RPot4;
         RPot4          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 5 plane
// =======================================
         SendID         = PotCen0 +  5*dz + Disp5;
         RecvID         = PotCen0 +  5*dz + Disp4;
         Temp           = RPot5;
         RPot5          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 6 plane
// =======================================
         SendID         = PotCen0 +  6*dz + Disp4;
         RecvID         = PotCen0 +  6*dz + Disp5;
         Temp           = RPot6;
         RPot6          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 7 plane
// =======================================
         SendID         = PotCen0 +  7*dz + Disp5;
         RecvID         = PotCen0 +  7*dz + Disp4;
         Temp           = RPot7;
         RPot7          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 8 plane
// =======================================
         SendID         = PotCen0 +  8*dz + Disp4;
         RecvID         = PotCen0 +  8*dz + Disp5;
         Temp           = RPot8;
         RPot8          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 9 plane
// =======================================
         SendID         = PotCen0 +  9*dz + Disp5;
         RecvID         = PotCen0 +  9*dz + Disp4;
         Temp           = RPot9;
         RPot9          = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 10 plane
// =======================================
         SendID         = PotCen0 + 10*dz + Disp4;
         RecvID         = PotCen0 + 10*dz + Disp5;
         Temp           = RPot10;
         RPot10         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 11 plane
// =======================================
         SendID         = PotCen0 + 11*dz + Disp5;
         RecvID         = PotCen0 + 11*dz + Disp4;
         Temp           = RPot11;
         RPot11         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 12 plane
// =======================================
         SendID         = PotCen0 + 12*dz + Disp4;
         RecvID         = PotCen0 + 12*dz + Disp5;
         Temp           = RPot12;
         RPot12         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 13 plane
// =======================================
         SendID         = PotCen0 + 13*dz + Disp5;
         RecvID         = PotCen0 + 13*dz + Disp4;
         Temp           = RPot13;
         RPot13         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;


#if ( POT_GHOST_SIZE == 5 )

//       z = 14 plane
// =======================================
         SendID         = PotCen0 + 14*dz + Disp4;
         RecvID         = PotCen0 + 14*dz + Disp5;
         Temp           = RPot14;
         RPot14         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

//       z = 15 plane
// =======================================
         SendID         = PotCen0 + 15*dz + Disp5;
         RecvID         = PotCen0 + 15*dz + Disp4;
         Temp           = RPot15;
         RPot15         = s_FPot[RecvID];
         SYNCTHREADS();
         s_FPot[SendID] = Temp;

#endif

         __syncthreads();



//       (b2-3). copy the +-yz-plane boundary potential stored in the temparory registers
//               back to shared memory
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//       -yz plane
// =======================================
         SendID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + __umul24(1+(tx<<1)+Disp5, POT_NXT_F/2);
         s_FPot[SendID] = Temp1;

//       +yz plane
// =======================================
         SendID         = __umul24(1+ty, POT_NXT_F*POT_NXT_F/2) + __umul24(1+(tx<<1)+Disp4, POT_NXT_F/2)
                          + (POT_NXT_F/2)-1;
         s_FPot[SendID] = Temp2;


         __syncthreads();



//       (b2-4). reset parameters for pass == 1 (the odd step)
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         PotCen   = PotCen0 + Disp2;
         RhoCen   = RhoCen0 + Disp2;
         Disp4    = Disp2;
         Disp5    = Disp1;
         Disp6    = Disp12;
         Disp8    = Disp13;
         Disp9    = Disp1;
         Disp10   = Disp7;
         Disp11   = Disp2;

      } // for (int pass=0; pass<2; pass++)



//    (b3). perform the reduction operation to get the sum of all residuals
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//    sum up the elements larger than FloorPow2 to ensure that the number of remaining elements is power-of-two
      if ( ID0 < Remain )  s_Residual_Total[ID0] += s_Residual_Total[ ID0 + FloorPow2 ];

//    parallel reduction
#     if ( POT_NTHREAD >= 1024 )
#     error : ERROR : POT_NTHREAD must < 1024 !!
#     endif

#     if ( POT_NTHREAD >= 512 )
      if ( ID0 < 256 )  s_Residual_Total[ID0] += s_Residual_Total[ ID0 + 256 ];  __syncthreads();
#     endif

#     if ( POT_NTHREAD >= 256 )
      if ( ID0 < 128 )  s_Residual_Total[ID0] += s_Residual_Total[ ID0 + 128 ];  __syncthreads();
#     endif

#     if ( POT_NTHREAD >= 128 )
      if ( ID0 <  64 )  s_Residual_Total[ID0] += s_Residual_Total[ ID0 +  64 ];  __syncthreads();
#     endif

//    adopting warp-synchronous mechanism
      if ( ID0 < 32 )
      {
//       declare volatile pointer to ensure that the operations are not reordered
         volatile real *s_Sum = s_Residual_Total;

         s_Sum[ID0] += s_Sum[ID0+32];  // here we have assumed that POT_NTHREAD >= 64
         s_Sum[ID0] += s_Sum[ID0+16];
         s_Sum[ID0] += s_Sum[ID0+ 8];
         s_Sum[ID0] += s_Sum[ID0+ 4];
         s_Sum[ID0] += s_Sum[ID0+ 2];
         s_Sum[ID0] += s_Sum[ID0+ 1];
      }
      __syncthreads();


//    (b4). termination criterion
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if ( Iter+1 >= Min_Iter  &&  s_Residual_Total[0] > Residual_Total_Old )    break;

      Residual_Total_Old = s_Residual_Total[0];

      __syncthreads();

   } // for (int Iter=0; Iter<Max_Iter; Iter++)



// c. store potential back to the global memory
// ---------------------------------------------------------------------------------------------------------

// (c1). internal potential stored in shared memory   -->   global memory  (store one z-slice at a time)
   const uint y = ( ID0 % (GRA_NXT*GRA_NXT/2) ) / (GRA_NXT/2);
   const uint x = ( ID0 % (GRA_NXT/2) );

   if ( ID0 < GRA_NXT*GRA_NXT/2 )
   {
      for (uint z=0; z<GRA_NXT; z++)
      {
         Disp4 = (y+z)&1;
         ID1   = __umul24( z, GRA_NXT*GRA_NXT ) + __umul24( y, GRA_NXT )
                 + 2*x+Disp4-2*(POT_USELESS2&Disp4)+POT_USELESS2;
         ID2   =   __umul24( z+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F*POT_NXT_F/2 )
                 + __umul24( y+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F/2 )
                 + x+(POT_GHOST_SIZE-GRA_GHOST_SIZE)/2+POT_USELESS2-(POT_USELESS2&Disp4);

         g_Pot_Array_Out[bx][ID1] = s_FPot[ID2];
      }
   }
   __syncthreads();


// (c2). internal potential stored in the registers   -->   global memory
#if ( POT_GHOST_SIZE == 4 )

   if (  ty >= POT_GHOST_SIZE-GRA_GHOST_SIZE-1  &&  ty <= POT_GHOST_SIZE+GRA_NXT-2  )
   {
#     if ( GRA_GHOST_SIZE == 2 )
      s_FPot[ PotCen0 +  1*dz ] = RPot1;
#     endif
#     if ( GRA_GHOST_SIZE >= 1 )
      s_FPot[ PotCen0 +  2*dz ] = RPot2;
#     endif
      s_FPot[ PotCen0 +  3*dz ] = RPot3;
      s_FPot[ PotCen0 +  4*dz ] = RPot4;
      s_FPot[ PotCen0 +  5*dz ] = RPot5;
      s_FPot[ PotCen0 +  6*dz ] = RPot6;
      s_FPot[ PotCen0 +  7*dz ] = RPot7;
      s_FPot[ PotCen0 +  8*dz ] = RPot8;
      s_FPot[ PotCen0 +  9*dz ] = RPot9;
      s_FPot[ PotCen0 + 10*dz ] = RPot10;
#     if ( GRA_GHOST_SIZE >= 1 )
      s_FPot[ PotCen0 + 11*dz ] = RPot11;
#     endif
#     if ( GRA_GHOST_SIZE == 2 )
      s_FPot[ PotCen0 + 12*dz ] = RPot12;
#     endif
   }
   __syncthreads();

   if ( ID0 < GRA_NXT*GRA_NXT/2 )
   {
      for (uint z=0; z<GRA_NXT; z++)
      {
         Disp4 = (y+z)&1;
         ID1   = __umul24( z, GRA_NXT*GRA_NXT ) + __umul24( y, GRA_NXT )
//               + 2*x+(Disp4^(GRA_GHOST_SIZE/2));
                 + 2*x + (  Disp4 ^ ( 1-(GRA_GHOST_SIZE&1) )  );
         ID2   =   __umul24( z+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F*POT_NXT_F/2 )
                 + __umul24( y+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F/2 )
//               + x+1-(Disp4&(GRA_GHOST_SIZE/2));
                 + x + 2 -(GRA_GHOST_SIZE+1)/2 - ( Disp4 & (1-GRA_GHOST_SIZE&1) );

         g_Pot_Array_Out[bx][ID1] = s_FPot[ID2];
      }
   }
   __syncthreads();

#else // #if ( POT_GHOST_SIZE == 4 )

   if (  ty >= POT_GHOST_SIZE-GRA_GHOST_SIZE-1  &&  ty <= POT_GHOST_SIZE+GRA_NXT-2  )
   {
#     if ( GRA_GHOST_SIZE == 2 )
      s_FPot[ PotCen0 +  2*dz ] = RPot2;
#     endif
#     if ( GRA_GHOST_SIZE >= 1 )
      s_FPot[ PotCen0 +  3*dz ] = RPot3;
#     endif
      s_FPot[ PotCen0 +  4*dz ] = RPot4;
      s_FPot[ PotCen0 +  5*dz ] = RPot5;
      s_FPot[ PotCen0 +  6*dz ] = RPot6;
      s_FPot[ PotCen0 +  7*dz ] = RPot7;
      s_FPot[ PotCen0 +  8*dz ] = RPot8;
      s_FPot[ PotCen0 +  9*dz ] = RPot9;
      s_FPot[ PotCen0 + 10*dz ] = RPot10;
      s_FPot[ PotCen0 + 11*dz ] = RPot11;
#     if ( GRA_GHOST_SIZE >= 1 )
      s_FPot[ PotCen0 + 12*dz ] = RPot12;
#     endif
#     if ( GRA_GHOST_SIZE == 2 )
      s_FPot[ PotCen0 + 13*dz ] = RPot13;
#     endif
   }
   __syncthreads();

   if ( ID0 < GRA_NXT*GRA_NXT/2 )
   {
      for (uint z=0; z<GRA_NXT; z++)
      {
         Disp4 = (y+z)&1;
         ID1   = __umul24( z, GRA_NXT*GRA_NXT ) + __umul24( y, GRA_NXT )
//               + 2*x+( (Disp4^1) ^ (GRA_GHOST_SIZE/2) );
                 + 2*x + (  (Disp4^1) ^ ( 1-(GRA_GHOST_SIZE&1) )  );
         ID2   =   __umul24( z+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F*POT_NXT_F/2 )
                 + __umul24( y+POT_GHOST_SIZE-GRA_GHOST_SIZE, POT_NXT_F/2 )
//               + x+1+( (Disp4^1) & (GRA_GHOST_SIZE&1) );
                 + x + 2 - (GRA_GHOST_SIZE+1)/2 + ( (Disp4^1) & (GRA_GHOST_SIZE&1) );

         g_Pot_Array_Out[bx][ID1] = s_FPot[ID2];
      }
   }

#endif // #if ( POT_GHOST_SIZE == 4 ) ... else ...

} // FUNCTION : CUPOT_PoissonSolver_SOR_16to18cube



#endif // #if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == SOR  &&  !defined USE_PSOLVER_10TO14 )
