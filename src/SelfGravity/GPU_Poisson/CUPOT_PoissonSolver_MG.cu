#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == MG )



#define POT_NXT_F    ( PATCH_SIZE+2*POT_GHOST_SIZE )  
#define POT_NTHREAD  ( POT_BLOCK_SIZE_X            )
#define POT_USELESS  ( POT_GHOST_SIZE%2            )

#if   ( POT_NXT_F == 18 )

//    for POT_NXT_F == 18, we reuse the same shared memory array due to the lack of shared memory
#     define REUSE_SHARED

#     define MAX_NLV             4U
#     define NBOTTOM_SMOOTH      1U
#     define NGRID_LV0          18U
#     define NGRID_LV1           9U
#     define NGRID_LV2           5U
#     define NGRID_LV3           3U
#elif ( POT_NXT_F == 16 )
#     define MAX_NLV             3U
#     define NBOTTOM_SMOOTH      7U
#     define NGRID_LV0          16U
#     define NGRID_LV1           8U
#     define NGRID_LV2           4U
#elif ( POT_NXT_F == 14 )
#     define MAX_NLV             3U
#     define NBOTTOM_SMOOTH      7U
#     define NGRID_LV0          14U
#     define NGRID_LV1           7U
#     define NGRID_LV2           4U
#elif ( POT_NXT_F == 12 )
#     define MAX_NLV             3U
#     define NBOTTOM_SMOOTH      1U
#     define NGRID_LV0          12U
#     define NGRID_LV1           6U
#     define NGRID_LV2           3U
#elif ( POT_NXT_F == 10 )
#     define MAX_NLV             3U
#     define NBOTTOM_SMOOTH      1U
#     define NGRID_LV0          10U
#     define NGRID_LV1           5U
#     define NGRID_LV2           3U
#else
#error ERROR : not supported POT_NXT_F
#endif

#if ( MAX_NLV != 3  &&  MAX_NLV != 4 )
#error ERROR : MAX_NLV != 3 or 4
#endif

// variables reside in constant memory
#include "CUPOT_PoissonSolver_SetConstMem.cu"

// prototype
static __device__ void LoadRho( const real *g_Rho, real *s_Rho, const real Poi_Coeff, const uint g_Idx0 );
static __device__ void Smoothing( real *Sol, const real *RHS, const real dh, const uint NGrid, const uint Idx0 );
static __device__ void ComputeDefect( const real *Sol, const real *RHS, real *Def, const real dh, 
                                      const uint NGrid, const uint Idx0 );
static __device__ void Restrict( const real *FData, real *CData, const uint NGrid_F, const uint NGrid_C,
                                 const uint Idx0 );
static __device__ void Prolongate_and_Correct( const real *CData, real *FData, const uint NGrid_C, 
                                               const uint NGrid_F, const uint FIdx0 );
static __device__ void EstimateError( const real *Sol, const real *RHS, const real dh, real *s_Error, 
                                      real *s_SolSum, const uint tid );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_PoissonSolver_MG
// Description :  GPU Poisson solver using the multigrid scheme 
//
// Note        :  a. Work for POT_GHOST_SIZE = 1, 2, 3, 4, 5 <--> POT_NXT_F = 10, 12, 14, 16, 18
//                b. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                c. Reference : Numerical Recipes, Chapter 20.6
//
// Parameter   :  g_Rho_Array       : Global memory array storing the input density 
//                g_Pot_Array_In    : Global memory array storing the input "coarse-grid" potential for ]
//                                    interpolation
//                g_Pot_Array_Out   : Global memory array to store the output potential
//                dh_Min            : Grid size of the input data
//                Max_Iter          : Maximum number of iterations for multigrid
//                NPre_Smooth       : Number of pre-smoothing steps for multigrid
//                NPost_Smooth      : Number of post-smoothing steps for multigrid
//                Tolerated_Error   : Maximum tolerated error for multigrid
//                Poi_Coeff         : Coefficient in front of the RHS in the Poisson eq.
//                IntScheme         : Interpolation scheme for potential
//                                    --> currently supported schemes include
//                                        INT_CQUAD : conservative quadratic interpolation 
//                                        INT_QUAD  : quadratic interpolation 
//---------------------------------------------------------------------------------------------------
__global__ void CUPOT_PoissonSolver_MG( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ], 
                                        const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ], 
                                              real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                        const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                                        const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                                        const IntScheme_t IntScheme )
{

   const uint bid = blockIdx.x;
   const uint tid = threadIdx.x; 
   const uint dy  = POT_NXT_F;
   const uint dz  = POT_NXT_F*POT_NXT_F;

   int  Iter;
   uint t, s_Idx;
   real dh[MAX_NLV];


// set the grid sizes at all different levels
#  if   ( MAX_NLV == 4 )
   const uint NGrid[MAX_NLV] = { NGRID_LV0, NGRID_LV1, NGRID_LV2, NGRID_LV3 };
#  elif ( MAX_NLV == 3 )
   const uint NGrid[MAX_NLV] = { NGRID_LV0, NGRID_LV1, NGRID_LV2 };
#  endif

   dh[0] = dh_Min;
   for (uint Lv=1U; Lv<MAX_NLV; Lv++)   dh[Lv] = dh_Min * ( NGrid[0] - 1U ) / ( NGrid[Lv] - 1U );


// allocate shared memory 
   __shared__ real s_Sol_Lv0[ NGRID_LV0*NGRID_LV0*NGRID_LV0 ];
#  ifndef FLOAT8
   __shared__ real s_RHS_Lv0[ NGRID_LV0*NGRID_LV0*NGRID_LV0 ];
   __shared__ real s_SolSum[POT_NTHREAD];
   __shared__ real s_Error [POT_NTHREAD];
#  else // shared memory is too small for double precision --> use global memory instead
   __shared__ real *s_RHS_Lv0;
   __shared__ real *s_SolSum;
   __shared__ real *s_Error;
   if ( tid == 0 )
   {
      s_RHS_Lv0 = (real*)malloc( sizeof(real*)*NGRID_LV0*NGRID_LV0*NGRID_LV0 );
      s_SolSum  = (real*)malloc( sizeof(real*)*POT_NTHREAD );
      s_Error   = (real*)malloc( sizeof(real*)*POT_NTHREAD );
   }
   __syncthreads();
#  ifdef GAMER_DEBUG
   if ( tid == 0 ) 
   {
      if ( s_RHS_Lv0 == NULL )
      {
         printf( "ERROR : dynamic global memory allocation for \"%s\" failed at block %d in \"%s\" !!\n",
                 "s_RHS_Lv0", bid, __FUNCTION__ );
         return; 
      }
      if ( s_SolSum == NULL )
      {
         printf( "ERROR : dynamic global memory allocation for \"%s\" failed at block %d in \"%s\" !!\n",
                 "s_SolSum", bid, __FUNCTION__ );
         return; 
      }
      if ( s_Error == NULL )
      {
         printf( "ERROR : dynamic global memory allocation for \"%s\" failed at block %d in \"%s\" !!\n",
                 "s_Error", bid, __FUNCTION__ );
         return; 
      }
   }
#  endif
#  endif // #ifndef FLOAT8 ... else ...
   real *s_Def_Lv0 = s_RHS_Lv0;  // s_Def_Lv0, s_CPot and RHS_Lv0 share the same shared-memory array
   real *s_CPot    = s_RHS_Lv0;

#  ifdef REUSE_SHARED 
// reuse the shared-memory arrays due to the lack of shared memory
   real *s_Sol_Lv1 = s_RHS_Lv0;
   real *s_Sol_Lv2 = s_Sol_Lv1 + NGRID_LV1*NGRID_LV1*NGRID_LV1;
   real *s_Sol_Lv3 = s_Sol_Lv2 + NGRID_LV2*NGRID_LV2*NGRID_LV2;
   real *s_RHS_Lv2 = s_Sol_Lv3 + NGRID_LV3*NGRID_LV3*NGRID_LV3;
   real *s_RHS_Lv3 = s_RHS_Lv2 + NGRID_LV2*NGRID_LV2*NGRID_LV2;
   real *s_Def_Lv1 = s_RHS_Lv3 + NGRID_LV3*NGRID_LV3*NGRID_LV3;
   real *s_Def_Lv2 = s_Def_Lv1 + NGRID_LV1*NGRID_LV1*NGRID_LV1;
   real *s_Def_Lv3 = s_Def_Lv2 + NGRID_LV2*NGRID_LV2*NGRID_LV2;

// use global memory for s_RHS_Lv1 because s_RHS_Lv1 and s_Def_Lv0 cannot share the same memory space
   __shared__ real *s_RHS_Lv1;
   if ( tid == 0 )   s_RHS_Lv1 = (real*) malloc( sizeof(real)*NGRID_LV1*NGRID_LV1*NGRID_LV1 );
   __syncthreads();
#  ifdef GAMER_DEBUG
   if ( tid == 0  &&  s_RHS_Lv1 == NULL )
   {
      printf( "ERROR : dynamic global memory allocation for \"%s\" failed at block %d in \"%s\" !!\n",
              "s_RHS_Lv1", bid, __FUNCTION__ );
      return; 
   }
#  endif
#  else // #ifdef REUSE_SHARED ... else ...
   __shared__ real s_Sol_Lv1[ NGRID_LV1*NGRID_LV1*NGRID_LV1 ];
   __shared__ real s_Sol_Lv2[ NGRID_LV2*NGRID_LV2*NGRID_LV2 ];
   __shared__ real s_RHS_Lv1[ NGRID_LV1*NGRID_LV1*NGRID_LV1 ];
   __shared__ real s_RHS_Lv2[ NGRID_LV2*NGRID_LV2*NGRID_LV2 ];
   __shared__ real s_Def_Lv1[ NGRID_LV1*NGRID_LV1*NGRID_LV1 ];
   __shared__ real s_Def_Lv2[ NGRID_LV2*NGRID_LV2*NGRID_LV2 ];
#  if ( MAX_NLV == 4 )
   __shared__ real s_Sol_Lv3[ NGRID_LV3*NGRID_LV3*NGRID_LV3 ];
   __shared__ real s_RHS_Lv3[ NGRID_LV3*NGRID_LV3*NGRID_LV3 ];
   __shared__ real s_Def_Lv3[ NGRID_LV3*NGRID_LV3*NGRID_LV3 ];
#  endif
#  endif // #ifdef REUSE_SHARED ... else ...

#  if   ( MAX_NLV == 4 )
   real *s_Sol[MAX_NLV] = { s_Sol_Lv0, s_Sol_Lv1, s_Sol_Lv2, s_Sol_Lv3 };
   real *s_RHS[MAX_NLV] = { s_RHS_Lv0, s_RHS_Lv1, s_RHS_Lv2, s_RHS_Lv3 };
   real *s_Def[MAX_NLV] = { s_Def_Lv0, s_Def_Lv1, s_Def_Lv2, s_Def_Lv3 };
#  elif ( MAX_NLV == 3 )
   real *s_Sol[MAX_NLV] = { s_Sol_Lv0, s_Sol_Lv1, s_Sol_Lv2 };
   real *s_RHS[MAX_NLV] = { s_RHS_Lv0, s_RHS_Lv1, s_RHS_Lv2 };
   real *s_Def[MAX_NLV] = { s_Def_Lv0, s_Def_Lv1, s_Def_Lv2 };
#  endif

   s_Error[0] = __FLT_MAX__;



// a. load the coarse-grid potential into shared memory
// -----------------------------------------------------------------------------------------------------------
   t = tid;
   while ( t < POT_NXT*POT_NXT*POT_NXT )
   {
      s_CPot[t] = g_Pot_Array_In[bid][t];    
      t += POT_NTHREAD; 
   }
   __syncthreads();



// b. evaluate the "fine-grid" potential by interpolation (as the initial guess and the B.C.)
// -----------------------------------------------------------------------------------------------------------
   const int N_CSlice = POT_NTHREAD / ( (POT_NXT-2)*(POT_NXT-2) );

   if ( tid < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   {
      const real Const_8   = 1.0/8.0;
      const real Const_64  = 1.0/64.0;
      const real Const_512 = 1.0/512.0;

      const int Cdx  = 1;
      const int Cdy  = POT_NXT;
      const int Cdz  = POT_NXT*POT_NXT;
      const int CIDx = 1 + tid % ( POT_NXT-2 );
      const int CIDy = 1 + (  tid % ( (POT_NXT-2)*(POT_NXT-2) )  ) / ( POT_NXT-2 );
      const int CIDz = 1 + tid / ( (POT_NXT-2)*(POT_NXT-2) );
      int       CID  = __mul24( CIDz, Cdz ) + __mul24( CIDy, Cdy ) + __mul24( CIDx, Cdx );

      const int Fdx  = 1;
      const int Fdy  = POT_NXT_F;
      const int Fdz  = POT_NXT_F*POT_NXT_F;
      const int FIDx = ( (CIDx-1)<<1 ) - POT_USELESS;
      const int FIDy = ( (CIDy-1)<<1 ) - POT_USELESS;
      int       FIDz = ( (CIDz-1)<<1 ) - POT_USELESS;
      int       FID  = __mul24( FIDz, Fdz ) + __mul24( FIDy, Fdy ) + __mul24( FIDx, Fdx );

      real TempFPot1, TempFPot2, TempFPot3, TempFPot4, TempFPot5, TempFPot6, TempFPot7, TempFPot8;
      real Slope_00, Slope_01, Slope_02, Slope_03, Slope_04, Slope_05, Slope_06, Slope_07;
      real Slope_08, Slope_09, Slope_10, Slope_11, Slope_12;
      int Idx, Idy, Idz, ii, jj, kk;


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
         
                  TempFPot1 += s_CPot[CID+kk+jj+ii] * Mm[Idz] * Mm[Idy] * Mm[Idx];
                  TempFPot2 += s_CPot[CID+kk+jj+ii] * Mm[Idz] * Mm[Idy] * Mp[Idx];
                  TempFPot3 += s_CPot[CID+kk+jj+ii] * Mm[Idz] * Mp[Idy] * Mm[Idx];
                  TempFPot4 += s_CPot[CID+kk+jj+ii] * Mm[Idz] * Mp[Idy] * Mp[Idx];
                  TempFPot5 += s_CPot[CID+kk+jj+ii] * Mp[Idz] * Mm[Idy] * Mm[Idx];
                  TempFPot6 += s_CPot[CID+kk+jj+ii] * Mp[Idz] * Mm[Idy] * Mp[Idx];
                  TempFPot7 += s_CPot[CID+kk+jj+ii] * Mp[Idz] * Mp[Idy] * Mm[Idx];
                  TempFPot8 += s_CPot[CID+kk+jj+ii] * Mp[Idz] * Mp[Idy] * Mp[Idx];

               }}}
            }
            break; // INT_QUAD
         
         } // switch ( IntScheme )



//       save data to the shared-memory array. 
//       Currently this part is highly diverge. However, since the interpolation takes much less time than the 
//       Poisson solver does, we have not yet tried to optimize this part
         if ( FIDz >= 0 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   s_Sol_Lv0[FID            ] = TempFPot1;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   s_Sol_Lv0[FID+Fdx        ] = TempFPot2;
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   s_Sol_Lv0[FID    +Fdy    ] = TempFPot3;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   s_Sol_Lv0[FID+Fdx+Fdy    ] = TempFPot4;
         }
         
         if ( FIDz <= POT_NXT_F-2 )
         {
            if ( FIDx >= 0            &&  FIDy >= 0           )   s_Sol_Lv0[FID        +Fdz] = TempFPot5;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy >= 0           )   s_Sol_Lv0[FID+Fdx    +Fdz] = TempFPot6;
            if ( FIDx >= 0            &&  FIDy <= POT_NXT_F-2 )   s_Sol_Lv0[FID    +Fdy+Fdz] = TempFPot7;
            if ( FIDx <= POT_NXT_F-2  &&  FIDy <= POT_NXT_F-2 )   s_Sol_Lv0[FID+Fdx+Fdy+Fdz] = TempFPot8;
         }
         
         CID  += __mul24(   N_CSlice, Cdz );
         FID  += __mul24( 2*N_CSlice, Fdz );
         FIDz += 2*N_CSlice;

      } // for (int z=CIDz; z<POT_NXT-1; z+=N_CSlice)
   } // if ( tid < N_CSlice*(POT_NXT-2)*(POT_NXT-2) )
   __syncthreads();



// c1. initialize s_Def_Lv{0-3} as zero (just to make sure that the boundary cells of s_Def_Lv{0-3} are zero)
//    (note that s_Def_Lv0 and s_CPot share the same array)
// -----------------------------------------------------------------------------------------------------------
#  ifndef REUSE_SHARED 
   t = tid;    while ( t < NGRID_LV0*NGRID_LV0*NGRID_LV0 )  {  s_Def_Lv0[t] = (real)0.0;  t += POT_NTHREAD; }
   t = tid;    while ( t < NGRID_LV1*NGRID_LV1*NGRID_LV1 )  {  s_Def_Lv1[t] = (real)0.0;  t += POT_NTHREAD; }
   t = tid;    while ( t < NGRID_LV2*NGRID_LV2*NGRID_LV2 )  {  s_Def_Lv2[t] = (real)0.0;  t += POT_NTHREAD; }
#  if ( MAX_NLV == 4 )
   t = tid;    while ( t < NGRID_LV3*NGRID_LV3*NGRID_LV3 )  {  s_Def_Lv3[t] = (real)0.0;  t += POT_NTHREAD; }
#  endif
   __syncthreads();
#  endif

// c2 load density into shared memory
   LoadRho( g_Rho_Array[bid], s_RHS_Lv0, Poi_Coeff, tid );



// d. use the MG scheme to evaluate potential
// -----------------------------------------------------------------------------------------------------------
   Iter = 0;

   while ( Iter < Max_Iter  &&  s_Error[0] > Tolerated_Error )
   {
//    V-cycle : finer --> coarser grids
      for (int Lv=0; Lv<MAX_NLV-1; Lv++)
      {
//       pre-smoothing (apply relaxation to compute solution/correction)
         for (int PreStep=0; PreStep<NPre_Smooth; PreStep++)
         Smoothing( s_Sol[Lv], s_RHS[Lv], dh[Lv], NGrid[Lv], tid );

//       compute defect
         ComputeDefect( s_Sol[Lv], s_RHS[Lv], s_Def[Lv], dh[Lv], NGrid[Lv], tid );

//       restrict defect (use as the RHS at the next level)
         Restrict( s_Def[Lv], s_RHS[Lv+1], NGrid[Lv], NGrid[Lv+1], tid );

//       initialize the correction at the next level to zero
         t = tid;  
         while ( t < NGrid[Lv+1]*NGrid[Lv+1]*NGrid[Lv+1] )  {  s_Sol[Lv+1][t] = (real)0.0;  t += POT_NTHREAD; }
         __syncthreads();
      }


//    calculate the correction at the bottom level
      for (int BottomStep=0; BottomStep<NBOTTOM_SMOOTH; BottomStep++)
      Smoothing( s_Sol[MAX_NLV-1], s_RHS[MAX_NLV-1], dh[MAX_NLV-1], NGrid[MAX_NLV-1], tid );


//    V-cycle : coarser --> finer grids
      for (int Lv=MAX_NLV-2; Lv>=0; Lv--)
      {
//       prolongate correction (from Lv+1 to Lv) and correct solution/correction at Lv
         Prolongate_and_Correct( s_Sol[Lv+1], s_Sol[Lv], NGrid[Lv+1], NGrid[Lv], tid );

//       load s_RHS_Lv0 from global memory since s_RHS_Lv0 and s_Def_Lv0 share the same shared-memory array
         if ( Lv == 0 )    LoadRho( g_Rho_Array[bid], s_RHS_Lv0, Poi_Coeff, tid );

//       post-smoothing (apply relaxation to compute solution/correction again)
         for (int PostStep=0; PostStep<NPost_Smooth; PostStep++)
         Smoothing( s_Sol[Lv], s_RHS[Lv], dh[Lv], NGrid[Lv], tid );
      }

//    estimate error
      EstimateError( s_Sol[0], s_RHS[0], dh[0], s_Error, s_SolSum, tid );
      Iter ++;

//    if ( tid == 0  &&  bid == 0 )    printf( "Patch %3d, Iter %3d: Error = %13.7e\n", bid, Iter, s_Error[0] );

   } // while ( Iter < Max_Iter  &&  Error > Tolerated_Error )



// e. store potential back to global memory 
// -----------------------------------------------------------------------------------------------------------
   t = tid;
   while ( t < GRA_NXT*GRA_NXT*GRA_NXT )
   { 
      s_Idx =   __umul24(  t/(GRA_NXT*GRA_NXT)         + POT_GHOST_SIZE - GRA_GHOST_SIZE,  dz  ) 
              + __umul24(  t%(GRA_NXT*GRA_NXT)/GRA_NXT + POT_GHOST_SIZE - GRA_GHOST_SIZE,  dy  )
              +            t%(GRA_NXT        )         + POT_GHOST_SIZE - GRA_GHOST_SIZE;

      g_Pot_Array_Out[bid][t] = s_Sol_Lv0[s_Idx];

      t += POT_NTHREAD; 
   }

// free memory
#  ifdef REUSE_SHARED 
   if ( tid == 0 )   free( s_RHS_Lv1 );
#  endif
#  ifdef FLOAT8
   if ( tid == 0 )
   {
      free( s_RHS_Lv0 );
      free( s_SolSum  );
      free( s_Error   );
   }
#  endif

} // FUNCTION : CUPOT_PoissonSolver_MG



//-------------------------------------------------------------------------------------------------------
// Function    :  Smoothing
// Description :  Use Gauss-Seidel method for smoothing
//
// Note        :  1. B.C. should be stored in the input array "Sol"
//                2. "__syncthreads" is invoked in the end of this function 
//
// Parameter   :  Sol   : 1D array to store the output solution to the Poisson equation
//                RHS   : 1D array storing the RHS of the Poisson equation
//                dh    : Grid size
//                NGrid : Number of cells in each spatial direction for Sol and RHS
//                Idx0  : Starting cell index
//-------------------------------------------------------------------------------------------------------
__device__ void Smoothing( real *Sol, const real *RHS, const real dh, const uint NGrid, const uint Idx0 )
{

   const real dh2      = dh*dh;
   const real One_Six  = (real)1.0/(real)6.0;
   const uint NGrid_m2 = NGrid - 2U;
   const uint NGrid2   = __umul24( NGrid, NGrid );
   const uint NGrid_2  = (NGrid_m2+1)>>1U;      // if NGrid is an odd number, one cell is padded
   const uint NGrid2_2 = __umul24( NGrid_m2, NGrid_2 );
   const uint NInner_2 = __umul24( NGrid_m2, NGrid2_2 );
   const uint di       = 1U;
   const uint dj       = NGrid;
   const uint dk       = NGrid2;

   uint i, j, k, ip, jp, kp, im, jm, km, ijk, pass_flip;
   uint Idx;


// odd-even ordering
   for (uint pass=0; pass<2U; pass++)
   {
      Idx       = Idx0;
      pass_flip = pass & 1U;

      while ( Idx < NInner_2 )
      {
         i = Idx%NGrid_2<<1U;
         j = Idx%NGrid2_2/NGrid_2;
         k = Idx/NGrid2_2;

         i += 1U + ( (j&1U)^(k&1U)^pass_flip );
         j++;
         k++;

         ijk = __umul24( k, NGrid2 ) + __umul24( j, NGrid ) + i;
         ip  = ijk + di;
         jp  = ijk + dj;
         kp  = ijk + dk;
         im  = ijk - di;
         jm  = ijk - dj;
         km  = ijk - dk;

//       update solution
//###OPTIMIZATION: try to optimize out this conditional operation (for the case that NGrid is odd)
#        if (  ( POT_NXT_F & (POT_NXT_F-1) ) != 0  )
         if ( i <= NGrid_m2 )
#        endif
         Sol[ijk] = One_Six*( Sol[kp] + Sol[km] + Sol[jp] + Sol[jm] + Sol[ip] + Sol[im] - dh2*RHS[ijk] );

         Idx += POT_NTHREAD;
      } // while ( Idx < NInner_2 )

      __syncthreads();

   } // for (int pass=0; pass<2; pass++)

} // FUNCTION : Smoothing



//-------------------------------------------------------------------------------------------------------
// Function    :  ComputeDefect
// Description :  Compute negative defect defined as "-(Laplacian(Sol)-RHS)"
//
// Note        :  1. B.C. should be stored in the input array "Sol"
//                2. It is assumed that the boundary values of Def have already been initialized as zero
//                   (unless REUSE_SHARED is defined)
//                3. "__syncthreads" is invoked in the end of this function 
//
// Parameter   :  Sol            : 1D array storing the input solution to the Poisson equation
//                RHS            : 1D array storing the RHS of the Poisson equation
//                Def            : 1D array to store the output defect
//                dh             : Grid size
//                NGrid          : Number of cells in each spatial direction for Sol and RHS
//                Idx0           : Starting cell index
//-------------------------------------------------------------------------------------------------------
__device__ void ComputeDefect( const real *Sol, const real *RHS, real *Def, const real dh, const uint NGrid, 
                               const uint Idx0 )
{

   const real _dh2       = (real)-1.0/(dh*dh);
   const uint NGrid2     = __umul24( NGrid, NGrid );
   const uint NGrid_m2   = NGrid - 2U;
   const uint NGrid_m2_2 = __umul24( NGrid_m2, NGrid_m2 );
   const uint NGrid_m2_3 = __umul24( NGrid_m2, NGrid_m2_2 );
   const uint di         = 1U;
   const uint dj         = NGrid;
   const uint dk         = NGrid2;

   uint i, j, k, ip, jp, kp, im, jm, km, ijk;
   uint Idx = Idx0;


   while ( Idx < NGrid_m2_3 )
   {
      i   = 1U + Idx%NGrid_m2;
      j   = 1U + Idx%NGrid_m2_2/NGrid_m2;
      k   = 1U + Idx/NGrid_m2_2;
      ijk = __umul24( k, NGrid2 ) + __umul24( j, NGrid ) + i;
      ip  = ijk + di;
      jp  = ijk + dj;
      kp  = ijk + dk;
      im  = ijk - di;
      jm  = ijk - dj;
      km  = ijk - dk;

      Def[ijk] = _dh2*( Sol[kp] + Sol[km] + Sol[jp] + Sol[jm] + Sol[ip] + Sol[im] - (real)6.0*Sol[ijk] )
                 + RHS[ijk];

      Idx += POT_NTHREAD;
   } // while ( Idx < NGrid_m2_3 )


// set the boundary values as zeros
#  ifdef REUSE_SHARED 
   Idx = Idx0;
   while ( Idx < NGrid2 )
   {  
      i = Idx%NGrid;
      j = Idx/NGrid;

      Def[ __umul24(       0U, NGrid2 ) + __umul24(        j, NGrid ) +        i ] = (real)0.0;
      Def[ __umul24( NGrid-1U, NGrid2 ) + __umul24(        j, NGrid ) +        i ] = (real)0.0;
      Def[ __umul24(        j, NGrid2 ) + __umul24(        i, NGrid ) +       0U ] = (real)0.0;
      Def[ __umul24(        j, NGrid2 ) + __umul24(        i, NGrid ) + NGrid-1U ] = (real)0.0;
      Def[ __umul24(        j, NGrid2 ) + __umul24(       0U, NGrid ) +        i ] = (real)0.0;
      Def[ __umul24(        j, NGrid2 ) + __umul24( NGrid-1U, NGrid ) +        i ] = (real)0.0;

      Idx += POT_NTHREAD;
   }
#  endif

   __syncthreads();

} // FUNCTION : ComputeDefect



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadRho
// Description :  Load density field from the global memory to the shared memory
//
// Note        :  1. "__syncthreads" is invoked in the end of this function 
//                2. Loaded data will be multiplied by "Poi_Coeff"
//                3. The size of the shared-memory array "s_Rho" is "RHO_NXT+2" (padded with one zero on each
//                   side in each direction)
//
// Parameter   :  g_Rho       : Global memory array storing the input density
//                s_Rho       : Shared memory array to store the density
//                Poi_Coeff   : Coefficient in front of density in the Poisson equation (4*Pi*Newton_G*a)
//                g_Idx0      : Starting read index from the global memory
//-------------------------------------------------------------------------------------------------------
__device__ void LoadRho( const real *g_Rho, real *s_Rho, const real Poi_Coeff, const uint g_Idx0 )
{

   uint  s_Idx, g_Idx = g_Idx0;
   uint3 g_Idx3D;

   while ( g_Idx < RHO_NXT*RHO_NXT*RHO_NXT )
   {
      g_Idx3D.x = g_Idx%RHO_NXT;
      g_Idx3D.y = g_Idx%(RHO_NXT*RHO_NXT)/RHO_NXT;
      g_Idx3D.z = g_Idx/(RHO_NXT*RHO_NXT);
      s_Idx     = __umul24(  __umul24( g_Idx3D.z+1U, NGRID_LV0 ) + g_Idx3D.y+1U, NGRID_LV0  ) + g_Idx3D.x+1U;

      s_Rho[s_Idx] = Poi_Coeff*g_Rho[g_Idx];

      g_Idx += POT_NTHREAD;
   }

   __syncthreads();

} // LoadRho



//-------------------------------------------------------------------------------------------------------
// Function    :  Restrict
// Description :  Restrict the input fine-grid data to get the coarse-grid data 
//
// Note        :  1. We assume that the input arrays follow the "finite-difference" fashion, in which the data
//                   are defined in the cell intersections instead of cell averages
//                   --> N^3 cells define a 3D grid with the size equal to (N-1)^3
//                2. Fine-grid and coarse-grid data at boundaries are assumed to be zero (because defect at 
//                   boundaries are always zero)
//                3. "__syncthreads" is invoked in the end of this function 
//
// Parameter   :  FData    : 1D array storing the input fine-grid data
//                CData    : 1D array to store the output coarse-grid data
//                NGrid_F  : Number of fine-grid cells in each spatial direction
//                NGrid_C  : Number of coarse-grid cells in each spatial direction
//                CIdx0    : Starting coarse-grid index
//-------------------------------------------------------------------------------------------------------
__device__ void Restrict( const real *FData, real *CData, const uint NGrid_F, const uint NGrid_C,
                          const uint CIdx0 )
{

   const real Ratio       = real(NGrid_F-1) / real(NGrid_C-1);
   const uint NGrid_F2    = __umul24( NGrid_F, NGrid_F );
   const uint NGrid_C2    = __umul24( NGrid_C, NGrid_C );
   const uint NGrid_Cm2   = NGrid_C - 2U;
   const uint NGrid_Cm2_2 = __umul24( NGrid_Cm2, NGrid_Cm2 );
   const uint NGrid_Cm2_3 = __umul24( NGrid_Cm2, NGrid_Cm2_2 );
   const uint Fdi         = 1U;
   const uint Fdj         = NGrid_F;
   const uint Fdk         = NGrid_F2;

   uint Ci, Cj, Ck, Fi, Fj, Fk, Cijk, Fijk;
   real x, y, z, Coeff_xm, Coeff_xc, Coeff_xp, Coeff_ym, Coeff_yc, Coeff_yp, Coeff_zm, Coeff_zc, Coeff_zp;
   uint CIdx = CIdx0;


   while ( CIdx < NGrid_Cm2_3 )
   {
      Ci       = 1U + CIdx%NGrid_Cm2;
      Cj       = 1U + CIdx%NGrid_Cm2_2/NGrid_Cm2;
      Ck       = 1U + CIdx/NGrid_Cm2_2;

      x        = Ci*Ratio;
      y        = Cj*Ratio;
      z        = Ck*Ratio;
      Fi       = uint( x + (real)0.5 );
      Fj       = uint( y + (real)0.5 );
      Fk       = uint( z + (real)0.5 );

      Cijk     = __umul24( Ck, NGrid_C2 ) + __umul24( Cj, NGrid_C ) + Ci;
      Fijk     = __umul24( Fk, NGrid_F2 ) + __umul24( Fj, NGrid_F ) + Fi;

      Coeff_xm = (real)0.5 * ( Fi + (real)0.5 - x );
      Coeff_ym = (real)0.5 * ( Fj + (real)0.5 - y );
      Coeff_zm = (real)0.5 * ( Fk + (real)0.5 - z );

      Coeff_xp = (real)0.5 - Coeff_xm;
      Coeff_yp = (real)0.5 - Coeff_ym;
      Coeff_zp = (real)0.5 - Coeff_zm;

      Coeff_xc = (real)0.5;
      Coeff_yc = (real)0.5;
      Coeff_zc = (real)0.5;
//    Coeff_xc = (real)1.0 - Coeff_xm - Coeff_xp;
//    Coeff_yc = (real)1.0 - Coeff_ym - Coeff_yp;
//    Coeff_zc = (real)1.0 - Coeff_zm - Coeff_zp;


//###OPTIMIZATION : follow the same strategy adopted in "Int_Quadratic"
      CData[Cijk] =       Coeff_zm * Coeff_ym * Coeff_xm * FData[ Fijk - Fdk - Fdj - Fdi ]
                       +  Coeff_zm * Coeff_ym * Coeff_xc * FData[ Fijk - Fdk - Fdj       ]
                       +  Coeff_zm * Coeff_ym * Coeff_xp * FData[ Fijk - Fdk - Fdj + Fdi ]
                       +  Coeff_zm * Coeff_yc * Coeff_xm * FData[ Fijk - Fdk       - Fdi ]
                       +  Coeff_zm * Coeff_yc * Coeff_xc * FData[ Fijk - Fdk             ]
                       +  Coeff_zm * Coeff_yc * Coeff_xp * FData[ Fijk - Fdk       + Fdi ]
                       +  Coeff_zm * Coeff_yp * Coeff_xm * FData[ Fijk - Fdk + Fdj - Fdi ]
                       +  Coeff_zm * Coeff_yp * Coeff_xc * FData[ Fijk - Fdk + Fdj       ]
                       +  Coeff_zm * Coeff_yp * Coeff_xp * FData[ Fijk - Fdk + Fdj + Fdi ]

                       +  Coeff_zc * Coeff_ym * Coeff_xm * FData[ Fijk       - Fdj - Fdi ]
                       +  Coeff_zc * Coeff_ym * Coeff_xc * FData[ Fijk       - Fdj       ]
                       +  Coeff_zc * Coeff_ym * Coeff_xp * FData[ Fijk       - Fdj + Fdi ]
                       +  Coeff_zc * Coeff_yc * Coeff_xm * FData[ Fijk             - Fdi ]
                       +  Coeff_zc * Coeff_yc * Coeff_xc * FData[ Fijk                   ]
                       +  Coeff_zc * Coeff_yc * Coeff_xp * FData[ Fijk             + Fdi ]
                       +  Coeff_zc * Coeff_yp * Coeff_xm * FData[ Fijk       + Fdj - Fdi ]
                       +  Coeff_zc * Coeff_yp * Coeff_xc * FData[ Fijk       + Fdj       ]
                       +  Coeff_zc * Coeff_yp * Coeff_xp * FData[ Fijk       + Fdj + Fdi ]

                       +  Coeff_zp * Coeff_ym * Coeff_xm * FData[ Fijk + Fdk - Fdj - Fdi ]
                       +  Coeff_zp * Coeff_ym * Coeff_xc * FData[ Fijk + Fdk - Fdj       ]
                       +  Coeff_zp * Coeff_ym * Coeff_xp * FData[ Fijk + Fdk - Fdj + Fdi ]
                       +  Coeff_zp * Coeff_yc * Coeff_xm * FData[ Fijk + Fdk       - Fdi ]
                       +  Coeff_zp * Coeff_yc * Coeff_xc * FData[ Fijk + Fdk             ]
                       +  Coeff_zp * Coeff_yc * Coeff_xp * FData[ Fijk + Fdk       + Fdi ]
                       +  Coeff_zp * Coeff_yp * Coeff_xm * FData[ Fijk + Fdk + Fdj - Fdi ]
                       +  Coeff_zp * Coeff_yp * Coeff_xc * FData[ Fijk + Fdk + Fdj       ]
                       +  Coeff_zp * Coeff_yp * Coeff_xp * FData[ Fijk + Fdk + Fdj + Fdi ];

//    coefficient adopted in Enzo which seems to give faster convergence rate
//    CData[Cijk] *= (real)0.52*Ratio;

      CIdx += POT_NTHREAD;

   } // while ( CIdx < NGrid_Cm2_3 )

   __syncthreads();

} // FUNCTION : Restrict



//-------------------------------------------------------------------------------------------------------
// Function    :  Prolongate_and_Correct
// Description :  Prolongate the input coarse-grid correction to correct the fine-grid solution/correction
//
// Note        :  1. We assume that the input arrays follow the "finite-difference" fashion, in which the data
//                   are defined in the cell intersections instead of cell averages
//                   --> N^3 cells define a 3D grid with the size equal to (N-1)^3
//                2. Boundary data of FData_1D are not corrected (since solution/correction at boundaries
//                   should be fixed
//
// Parameter   :  CData    : 1D array storing the input coarse-grid data
//                FData    : 1D array to store the output fine-grid data
//                NGrid_C  : Number of coarse-grid cells in each spatial direction
//                NGrid_F  : Number of fine-grid cells in each spatial direction
//                FIdx0    : Starting fine-grid index
//-------------------------------------------------------------------------------------------------------
__device__ void Prolongate_and_Correct( const real *CData, real *FData, const uint NGrid_C, const uint NGrid_F, 
                                        const uint FIdx0 )
{

   const real Ratio       = real(NGrid_C-1) / real(NGrid_F-1);
   const uint NGrid_C2    = __umul24( NGrid_C, NGrid_C );
   const uint NGrid_F2    = __umul24( NGrid_F, NGrid_F );
   const uint NGrid_Fm2   = NGrid_F - 2U;
   const uint NGrid_Fm2_2 = __umul24( NGrid_Fm2, NGrid_Fm2 );
   const uint NGrid_Fm2_3 = __umul24( NGrid_Fm2, NGrid_Fm2_2 );
   const uint Cdi         = 1U;
   const uint Cdj         = NGrid_C;
   const uint Cdk         = NGrid_C2;

   uint Ci, Cj, Ck, Fi, Fj, Fk, Cijk, Fijk;
   real x, y, z, Coeff_xm, Coeff_xp, Coeff_ym, Coeff_yp, Coeff_zm, Coeff_zp;
   uint FIdx = FIdx0;


   while ( FIdx < NGrid_Fm2_3 )
   {
      Fi       = 1U + FIdx%NGrid_Fm2;
      Fj       = 1U + FIdx%NGrid_Fm2_2/NGrid_Fm2;
      Fk       = 1U + FIdx/NGrid_Fm2_2;

      x        = Fi*Ratio;
      y        = Fj*Ratio;
      z        = Fk*Ratio;
      Ci       = uint( x );
      Cj       = uint( y );
      Ck       = uint( z );

      Cijk     = __umul24( Ck, NGrid_C2 ) + __umul24( Cj, NGrid_C ) + Ci;
      Fijk     = __umul24( Fk, NGrid_F2 ) + __umul24( Fj, NGrid_F ) + Fi;

      Coeff_xm = Ci + 1U - x;
      Coeff_ym = Cj + 1U - y;
      Coeff_zm = Ck + 1U - z;
      Coeff_xp = (real)1.0 - Coeff_xm;
      Coeff_yp = (real)1.0 - Coeff_ym;
      Coeff_zp = (real)1.0 - Coeff_zm;

      FData[Fijk] +=    Coeff_zm * Coeff_ym * Coeff_xm * CData[ Cijk                   ]
                     +  Coeff_zm * Coeff_ym * Coeff_xp * CData[ Cijk             + Cdi ]
                     +  Coeff_zm * Coeff_yp * Coeff_xm * CData[ Cijk       + Cdj       ]
                     +  Coeff_zp * Coeff_ym * Coeff_xm * CData[ Cijk + Cdk             ]
                     +  Coeff_zm * Coeff_yp * Coeff_xp * CData[ Cijk       + Cdj + Cdi ]
                     +  Coeff_zp * Coeff_yp * Coeff_xm * CData[ Cijk + Cdk + Cdj       ]
                     +  Coeff_zp * Coeff_ym * Coeff_xp * CData[ Cijk + Cdk       + Cdi ]
                     +  Coeff_zp * Coeff_yp * Coeff_xp * CData[ Cijk + Cdk + Cdj + Cdi ];

      FIdx += POT_NTHREAD;

   } // while ( FIdx < NGrid_Fm2_3 )

   __syncthreads();

} // FUNCTION : Prolongate_and_Correct



//-------------------------------------------------------------------------------------------------------
// Function    :  EstimateError
// Description :  Estimate the L1 error 
//
// Note        :  1. "__syncthreads" is invoked in the end of this function 
//                2. Shared-memory arrays "s_Error" and "s_SolSum" are used for GPU reduction
//
// Parameter   :  Sol      : 1D array storing the input solution to the Poisson equation
//                RHS      : 1D array storing the RHS of the Poisson equation
//                dh       : Grid size
//                s_Error  : Shared-memory array to store the L1 error
//                s_SolSum : Shared-memroy array to store the sum of solution
//                tid      : Thread index
//-------------------------------------------------------------------------------------------------------
__device__ void EstimateError( const real *Sol, const real *RHS, const real dh, real *s_Error, real *s_SolSum,
                               const uint tid )
{

#  define NGRID_M2 ( NGRID_LV0 - 2U )

   const real dh2       = dh*dh;
   const real _dh2      = (real)-1.0/dh2;
   const uint di        = 1U;
   const uint dj        = NGRID_LV0;
   const uint dk        = NGRID_LV0*NGRID_LV0;
   const uint FloorPow2 = 1<<(31-__clz(POT_NTHREAD) ); // largest power-of-two value not greater than POT_NTHREAD
   const uint Remain    = POT_NTHREAD - FloorPow2;

   uint i, j, k, ip, jp, kp, im, jm, km, ijk; 
   uint Idx = tid;

   s_Error [tid] = (real)0.0;
   s_SolSum[tid] = (real)0.0;


// 1. calculate defect
   while ( Idx < NGRID_M2*NGRID_M2*NGRID_M2 )
   {
      i   = 1U + Idx%NGRID_M2;
      j   = 1U + Idx%(NGRID_M2*NGRID_M2)/NGRID_M2;
      k   = 1U + Idx/(NGRID_M2*NGRID_M2);
      ijk = __umul24( k, NGRID_LV0*NGRID_LV0 ) + __umul24( j, NGRID_LV0 ) + i;
      ip  = ijk + di;
      jp  = ijk + dj;
      kp  = ijk + dk;
      im  = ijk - di;
      jm  = ijk - dj;
      km  = ijk - dk;

      s_Error [tid] += FABS(  _dh2*( Sol[kp]+Sol[km]+Sol[jp]+Sol[jm]+Sol[ip]+Sol[im]-(real)6.0*Sol[ijk] ) 
                              + RHS[ijk]  );
      s_SolSum[tid] += FABS( Sol[ijk] );

      Idx += POT_NTHREAD;
   } // while ( Idx < NGRID_M2*NGRID_M2*NGRID_M2 )

   __syncthreads();


// 2. perform the reduction operation to get the L1 error

// first sum up the elements larger than FloorPow2 to ensure that the number of remaining elements is power-of-two
   if ( tid < Remain )  
   {
      s_Error [tid] += s_Error [ tid + FloorPow2 ];
      s_SolSum[tid] += s_SolSum[ tid + FloorPow2 ];
   }

// parallel reduction
#  if ( POT_NTHREAD >= 1024 )
#  error : ERROR : POT_NTHREAD must < 1024 !!
#  endif

#  if ( POT_NTHREAD >= 512 )
   if ( tid < 256 )  
   {
      s_Error [tid] += s_Error [ tid + 256 ];
      s_SolSum[tid] += s_SolSum[ tid + 256 ];
   }
   __syncthreads();
#  endif

#  if ( POT_NTHREAD >= 256 ) 
   if ( tid < 128 )
   {
      s_Error [tid] += s_Error [ tid + 128 ];
      s_SolSum[tid] += s_SolSum[ tid + 128 ];
   }
   __syncthreads();
#  endif

#  if ( POT_NTHREAD >= 128 )
   if ( tid <  64 )
   {
      s_Error [tid] += s_Error [ tid + 64 ];
      s_SolSum[tid] += s_SolSum[ tid + 64 ];
   }
   __syncthreads();
#  endif

// adopting warp-synchronous mechanism
   if ( tid < 32 ) 
   {  
//    declare volatile pointer to ensure that the operations are not reordered
      volatile real *s_vErr = s_Error;   
      volatile real *s_vSol = s_SolSum;   

      s_vErr[tid] += s_vErr[tid+32];   // here we have assumed that POT_NTHREAD >= 64
      s_vErr[tid] += s_vErr[tid+16];
      s_vErr[tid] += s_vErr[tid+ 8];
      s_vErr[tid] += s_vErr[tid+ 4];
      s_vErr[tid] += s_vErr[tid+ 2];
      s_vErr[tid] += s_vErr[tid+ 1];

      s_vSol[tid] += s_vSol[tid+32];
      s_vSol[tid] += s_vSol[tid+16];
      s_vSol[tid] += s_vSol[tid+ 8];
      s_vSol[tid] += s_vSol[tid+ 4];
      s_vSol[tid] += s_vSol[tid+ 2];
      s_vSol[tid] += s_vSol[tid+ 1];

      s_vErr[tid] = dh2*s_vErr[tid]/s_vSol[tid];
   }
   __syncthreads();

#  undef NGRID_M2

} // FUNCTION : EstimateError



#endif // #if ( defined GRAVITY  &&  defined GPU  &&  POT_SCHEME == MG )
