#include "CUPOT.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"


// parallel reduction routine
#define RED_NTHREAD  ( DT_GRA_BLOCK_SIZE )
#define RED_MAX

#ifdef DT_GRA_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif


// variables reside in constant memory
__constant__ double c_ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_SetConstMem_dtSolver_HydroGravity
// Description :  Set the constant memory used by CUPOT_dtSolver_HydroGravity()
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUPOT_SetConstMem_dtSolver_HydroGravity( double h_ExtAcc_AuxArray[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( c_ExtAcc_AuxArray, h_ExtAcc_AuxArray, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_SetConstMem_dtSolver_HydroGravity

#endif // #ifdef __CUDACC__




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_dtSolver_HydroGravity
// Description :  Estimate the evolution time-step (dt) required for the hydro gravity solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the free-fall time of the maximum gravitational acceleration
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_dt_Array        : Array to store the minimum dt in each target patch
//                g_Pot_Array       : Array storing the prepared potential data of each target patch
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                NPatchGroup       : Number of target patch groups (for CPU only)
//                dh                : Cell size
//                Safety            : dt safety factor
//                P5_Gradient       : Use 5-point stencil to evaluate the potential gradient
//                GravityType       : Types of gravity --> self-gravity, external gravity, both
//                c_ExtAcc_AuxArray : Auxiliary array for adding external acceleration (for CPU only)
//                                    --> When using GPU, this array is stored in the constant memory and does
//                                        not need to be passed as a function argument
//                                        --> Declared on top of this file with the prefix "c_" to
//                                            highlight that this is a constant variable on GPU
//                ExtAcc_Time       : Physical time for adding the external acceleration
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_dtSolver_HydroGravity( real g_dt_Array[], const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                  const double g_Corner_Array[][3],
                                  const real dh, const real Safety, const bool P5_Gradient,
                                  const OptGravityType_t GravityType,
                                  const double ExtAcc_Time )
#else
void CPU_dtSolver_HydroGravity  ( real g_dt_Array[], const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                  const double g_Corner_Array[][3], const int NPatchGroup,
                                  const real dh, const real Safety, const bool P5_Gradient,
                                  const OptGravityType_t GravityType, const double c_ExtAcc_AuxArray[],
                                  const double ExtAcc_Time )
#endif
{

// check
#  ifdef GAMER_DEBUG
   if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  &&  ExtAcc_Time < 0.0 )
      printf( "ERROR : incorrect ExtAcc_Time (%14.7e) !!\n", ExtAcc_Time );
#  endif


   const real dh2         = (real)2.0*dh;
   const real Gra_Const   = ( P5_Gradient ) ? (real)-1.0/((real)12.0*dh) : (real)-1.0/((real)2.0*dh);
   const int  PS1_sqr     = SQR(PS1);
   const int  didx_pot[3] = { 1, GRA_NXT, SQR(GRA_NXT) };


// load potential from global to shared memory to improve the GPU performance
#  ifdef __CUDACC__
   __shared__ real s_Pot[ CUBE(GRA_NXT) ];

   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      for (int t=threadIdx.x; t<CUBE(GRA_NXT); t+=DT_GRA_BLOCK_SIZE)
         s_Pot[t] = g_Pot_Array[blockIdx.x][t];
   }

   __syncthreads();
#  endif // #ifdef __CUDACC__


// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//     to work on different patches
#  ifdef __CUDACC__
   const int P = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatchGroup*8; P++)
#  endif
   {
//    point to the potential array of the target patch
#     ifdef __CUDACC__
      const real *const Pot = s_Pot;
#     else
      const real *const Pot = g_Pot_Array[P];
#     endif

      real AccMax = (real)0.0;


//    loop over all cells of the target patch
      CGPU_LOOP( t, CUBE(PS1) )
      {
         const int i_ext   = t % PS1;
         const int j_ext   = t % PS1_sqr / PS1;
         const int k_ext   = t / PS1_sqr;

         const int i_pot   = i_ext + GRA_GHOST_SIZE;
         const int j_pot   = j_ext + GRA_GHOST_SIZE;
         const int k_pot   = k_ext + GRA_GHOST_SIZE;
         const int idx_pot = IDX321( i_pot, j_pot, k_pot, GRA_NXT, GRA_NXT );

         real Acc[3] = { 0.0, 0.0, 0.0 };

//       external gravity
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         {
            double x, y, z;

            x = g_Corner_Array[P][0] + double(i_ext*dh);
            y = g_Corner_Array[P][1] + double(j_ext*dh);
            z = g_Corner_Array[P][2] + double(k_ext*dh);

            ExternalAcc( Acc, x, y, z, ExtAcc_Time, c_ExtAcc_AuxArray );
         }

//       self-gravity
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
         {
            const int idx_xp1 = idx_pot + didx_pot[0];
            const int idx_yp1 = idx_pot + didx_pot[1];
            const int idx_zp1 = idx_pot + didx_pot[2];
            const int idx_xm1 = idx_pot - didx_pot[0];
            const int idx_ym1 = idx_pot - didx_pot[1];
            const int idx_zm1 = idx_pot - didx_pot[2];

            if ( P5_Gradient )   // 5-point gradient
            {
               const real Const_8 = (real)8.0;

               const int idx_xp2 = idx_xp1 + didx_pot[0];
               const int idx_yp2 = idx_yp1 + didx_pot[1];
               const int idx_zp2 = idx_zp1 + didx_pot[2];
               const int idx_xm2 = idx_xm1 - didx_pot[0];
               const int idx_ym2 = idx_ym1 - didx_pot[1];
               const int idx_zm2 = idx_zm1 - didx_pot[2];

               Acc[0] += Gra_Const * ( - Pot[idx_xp2] + Pot[idx_xm2] + Const_8*Pot[idx_xp1] - Const_8*Pot[idx_xm1] );
               Acc[1] += Gra_Const * ( - Pot[idx_yp2] + Pot[idx_ym2] + Const_8*Pot[idx_yp1] - Const_8*Pot[idx_ym1] );
               Acc[2] += Gra_Const * ( - Pot[idx_zp2] + Pot[idx_zm2] + Const_8*Pot[idx_zp1] - Const_8*Pot[idx_zm1] );
            }

            else                 // 3-point gradient
            {
               Acc[0] += Gra_Const * ( Pot[idx_xp1] - Pot[idx_xm1] );
               Acc[1] += Gra_Const * ( Pot[idx_yp1] - Pot[idx_ym1] );
               Acc[2] += Gra_Const * ( Pot[idx_zp1] - Pot[idx_zm1] );
            }
         } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )

//       get the maximum acceleration
         for (int d=0; d<3; d++)    AccMax = FMAX( AccMax, FABS(Acc[d]) );
      } // CGPU_LOOP( t, CUBE(PS1) )


//    get the minimum dt
//    perform parallel reduction to get the maximum acceleration in each thread block
//    --> store in the thread 0
#     ifdef __CUDACC__
#     ifdef DT_GRA_USE_SHUFFLE
      AccMax = BlockReduction_Shuffle ( AccMax );
#     else
      AccMax = BlockReduction_WarpSync( AccMax );
#     endif
      if ( threadIdx.x == 0 )
#     endif // #ifdef __CUDACC__
      g_dt_Array[P] = Safety*SQRT( dh2/AccMax );

   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_dtSolver_HydroGravity



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
