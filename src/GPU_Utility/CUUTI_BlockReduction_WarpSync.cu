#include "Macro.h"

#ifdef GPU


// check
// one must define RED_NTHREAD for the reduction kernel in advance since we use the static shared memory
#ifndef RED_NTHREAD
#  error : ERROR : RED_NTHREAD is not defined in BlockReduction_WarpSync !!
#endif


// define the reduction operation here
#if   defined RED_SUM
#  define RED( a, b )   ( (a) + (b) )
#elif defined RED_MAX
#  define RED( a, b )   MAX( (a), (b) )
#elif defined RED_MIN
#  define RED( a, b )   MIN( (a), (b) )
#else
#  error : undefined reduction operation !!
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  BlockReduction_WarpSync
// Description :  GPU reduction within each thread block using the implicit warp synchronization
//
// Note        :  1. Mainly used for the Fermi GPUs
//                   --> For Kepler or later GPUs, it's recommended to use BlockReduction_Shuffle() instead
//                2. Reference: http://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/reduction/doc/reduction.pdf
//                   --> By Mark Harris
//                3. Assuming warp size == 32
//                4. Must define RED_NTHREAD in advance since we use the static shared memory
//                   --> RED_NTHREAD must < 2048 && >= 64
//                5. Must define either RED_SUM, RED_MAX, or RED_MIN in advance to determine the reduction operation
//                6. Only thread 0 will hold the correct result after calling this function
//
// Parameter   :  val : Per-thread value for the reduction
//
// Return value:  Reduction of "val"
//---------------------------------------------------------------------------------------------------
__inline__ __device__
real BlockReduction_WarpSync( real val )
{

   const uint tid_x     = threadIdx.x;
   const uint tid_y     = threadIdx.y;
   const uint tid_z     = threadIdx.z;
   const uint bdim_x    = blockDim.x;
   const uint bdim_y    = blockDim.y;
   const uint ID        = __umul24( tid_z, __umul24(bdim_x,bdim_y) ) + __umul24( tid_y, bdim_x ) + tid_x;
   const uint FloorPow2 = 1 << ( 31-__clz(RED_NTHREAD) );   // largest power-of-two value not greater than RED_NTHREAD
   const uint Remain    = RED_NTHREAD - FloorPow2;

   __shared__ real s_Reduction[RED_NTHREAD];


// store values for the reduction to the shared memory
   s_Reduction[ID] = val;
   __syncthreads();


// perform reduction for the elements larger than FloorPow2 to ensure that the number of remaining elements is power-of-two
   if ( ID < Remain )   s_Reduction[ID] = RED( s_Reduction[ID], s_Reduction[ ID + FloorPow2 ] );
   __syncthreads();


// parallel reduction with the shared memory
#  if ( RED_NTHREAD >= 2048 )
#  error : ERROR : RED_NTHREAD >= 2048 !!
#  endif

#  if ( RED_NTHREAD >= 1024 )
   if ( ID < 512 )   s_Reduction[ID] = RED( s_Reduction[ID], s_Reduction[ ID + 512 ] );   __syncthreads();
#  endif

#  if ( RED_NTHREAD >= 512 )
   if ( ID < 256 )   s_Reduction[ID] = RED( s_Reduction[ID], s_Reduction[ ID + 256 ] );   __syncthreads();
#  endif

#  if ( RED_NTHREAD >= 256 )
   if ( ID < 128 )   s_Reduction[ID] = RED( s_Reduction[ID], s_Reduction[ ID + 128 ] );   __syncthreads();
#  endif

#  if ( RED_NTHREAD >= 128 )
   if ( ID <  64 )   s_Reduction[ID] = RED( s_Reduction[ID], s_Reduction[ ID +  64 ] );   __syncthreads();
#  endif


// parallel reduction with the warp-synchronous mechanism (assuming warpSize == 32)
#  if ( WARP_SIZE != 32 )
#     error : ERROR : WARP_SIZE != 32 !!
#  endif

   if ( ID < WARP_SIZE )
   {
//    declare volatile pointer to ensure that the operations are not reordered
      volatile real *s_Reduction_Volatile = s_Reduction;

//    here we have assumed that RED_NTHREAD >= 64
#     if ( RED_NTHREAD < 64 )
#       error : ERROR : RED_NTHREAD < 64 !!
#     endif

      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+32] );
      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+16] );
      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+ 8] );
      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+ 4] );
      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+ 2] );
      s_Reduction_Volatile[ID] = RED( s_Reduction_Volatile[ID], s_Reduction_Volatile[ID+ 1] );
   }
   __syncthreads();


   return s_Reduction[0];

} // FUNCTION : BlockReduction_WarpSync

#undef RED


#endif // #ifdef GPU
