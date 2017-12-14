#include "Macro.h"

#ifdef GPU


// checks
// one must define RED_NTHREAD for the reduction kernel in advance since we use the static shared memory
#ifndef RED_NTHREAD
#  error : ERROR : RED_NTHREAD is not defined in BlockReduction_Shuffle !!
#endif

// WARP_SIZE must be defined to be the same as the CUDA predefined constant "warpSize"
#ifndef WARP_SIZE
#  error : ERROR : WARP_SIZE is not defined in BlockReduction_Shuffle !!
#endif

// RED_NTHREAD must be a multiple of WARP_SIZE
#if ( RED_NTHREAD % WARP_SIZE != 0 )
#  error : ERROR : RED_NTHREAD must be a multiple of WARP_SIZE !!
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
// Function    :  WarpReduction_Shuffle
// Description :  GPU reduction within each warp using the register shuffling
//
// Note        :  1. Invoked by BlockReduction_Shuffle
//                2. Only thread 0 will hold the correct result
//
// Parameter   :  val : Per-thread value for the reduction
//
// Return value:  Reduction of "val" within each warp
//---------------------------------------------------------------------------------------------------
__inline__ __device__
real WarpReduction_Shuffle( real val )
{

   for (int offset=WARP_SIZE/2; offset>0; offset/=2)
   {
//    this line somehow fails on K20X for RED_MAX (and RED_MIN likely)
//    --> perhaps it's because when using "val = (val > __shfl_down(val,offset) ) ? val : __shfl_down(val,offset);"
//        the second " __shfl_down(val,offset)" becomes ill-defined since "val" in other threads might be modified in advance
//        if these threads have "(val > __shfl_down(val,offset) )"
//    val = RED( val, __shfl_down(val,offset,WARP_SIZE) );

//    use this approach instead to invoke "__shfl_down(val,offset, WARP_SIZE)" only once
      const real tmp = __shfl_down( val, offset, WARP_SIZE );
      val = RED( val, tmp );
   }

   return val;

} // FUNCTION : WarpReduction_Shuffle



//-------------------------------------------------------------------------------------------------------
// Function    :  BlockReduction_Shuffle
// Description :  GPU reduction within each thread block using the register shuffling
//
// Note        :  1. Improve reduction performance by register shuffling
//                2. Reference: https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
//                   --> By Justin Luitjens
//                3. Assuming warp size == 32
//                4. Must define RED_NTHREAD in advance since we use the static shared memory
//                   --> RED_NTHREAD must be a multiple of the warp size
//                5. Must define either RED_SUM, RED_MAX, or RED_MIN in advance to determine the reduction operation
//                6. Only thread 0 will hold the correct result after calling this function
//
// Parameter   :  val : Per-thread value for the reduction
//
// Return value:  Reduction of "val" within each thread block
//---------------------------------------------------------------------------------------------------
__inline__ __device__
real BlockReduction_Shuffle( real val )
{

   const uint tid_x   = threadIdx.x;
   const uint tid_y   = threadIdx.y;
   const uint tid_z   = threadIdx.z;
   const uint bdim_x  = blockDim.x;
   const uint bdim_y  = blockDim.y;
   const uint ID      = __umul24( tid_z, __umul24(bdim_x,bdim_y) ) + __umul24( tid_y, bdim_x ) + tid_x;
   const int lane     = ID % WARP_SIZE;         // local lane ID within a warp [0 ... WARP_SIZE-1]
   const int wid      = ID / WARP_SIZE;         // warp ID
   const int MaxNWarp = 32;                     // maximum number of warps allowed == MaxBlockSize/WARP_SIZE == 1024/32 == 32
                                                // --> all current compute capabilities have MaxBlockSize==1024 and  WARP_SIZE==32
   const int NWarp    = RED_NTHREAD/WARP_SIZE;  // actual number of warps (which must be <= WARP_SIZE since we apply the
                                                // final reduction only to the first warp)

   static __shared__ real shared[MaxNWarp];     // maximum shared memory required for 32 partial sums (must be <= WARP_SIZE)

// perform reduction within each warp
   val = WarpReduction_Shuffle( val );

// write reduced value to the shared memory
   if ( lane == 0 )  shared[wid] = val;

// wait for all partial reductions
   __syncthreads();

// here we have assumed that NWarp < WARP_SIZE
   if ( wid == 0 )
   {
//    read from the shared memory only if that warp exists
      val = ( ID < NWarp ) ? shared[lane] :
#     if   defined RED_SUM
                             (real)0.0;
#     elif defined RED_MAX
                             (real)-HUGE_NUMBER;
#     elif defined RED_MIN
                             (real)+HUGE_NUMBER;
#     else
#       error : undefined reduction operation !!
#     endif

//    final reduction within first warp
      val = WarpReduction_Shuffle( val );
   }

   return val;

} // FUNCTION : BlockReduction_Shuffle



#endif // #ifdef GPU
