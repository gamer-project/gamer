#include "Copyright.h"
#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU )




//-------------------------------------------------------------------------------------------------------
// Function    :  WarpReductionWithShuffle
// Description :  GPU reduction using register shuffling. Sum up the elements in one warp.
//
// Note        :  1. Increase reduction speed by using registers
//                2. Reference: https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
//                   --> By Justin Luitjens
//                3. Only thread 0 will hold the result
//
// Parameter   :  val   : Value that each thread holds that will be added together
//
// Return      :  Sum of the values in the same warp (for thread 0 only)
//---------------------------------------------------------------------------------------------------
__inline__ __device__
real WarpReductionWithShuffle( real val )
{

   for (int offset=warpSize/2; offset>0; offset/=2)
      val += __shfl_down( val, offset, warpSize );

   return val;

} // FUNCTION : WarpReductionWithShuffle



//-------------------------------------------------------------------------------------------------------
// Function    :  BlockReductionWithShuffle
// Description :  GPU reduction using register shuffling. Sum up the elements in one block.
//
// Note        :  1. Increasesreduction speed by using registers
//                2. Reference: https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
//                   --> By Justin Luitjens
//                3. Works only if the total number of threads per block (i.e., the CUDA block size) is a multiple
//                   of warpSize (32)
//                   --> Always satisfied for PATCH_SIZE==8 and POT_GHOST_SIZE==5, for which
//                       block size == RHO_NXT/2 * RHO_NXT * POT_BLOCK_SIZE_Z = 8*16*POT_BLOCK_SIZE_Z = 128*POT_BLOCK_SIZE_Z
//
// Parameter   :  val   : Value that each thread holds that will be added together
//                ID    : Thread index
//
// Return value:  sum of the block
//---------------------------------------------------------------------------------------------------
__inline__ __device__
real BlockReductionWithShuffle( real val, const int ID )
{

// const int lane         = ID % warpSize;         // local lane ID within a warp [0 ... warpSize-1]
// const int wid          = ID / warpSize;         // warp ID
   const int lane         = ID & 0x1f;             // optimized for warpSize == 32
   const int wid          = ID >> 5;               // optimized for warpSize == 32
   const int MaxNWarp     = 32;                    // maximum number of warps allowed == MaxBlockSize/warpSize == 1024/32 == 32
                                                   // --> all current compute capabilities have MaxBlockSize==1024 warpSize==32
// const int NWarp        = POT_NTHREAD/warpSize;  // actual number of warps (which must be <= warpSize since we apply the
                                                   // final reduce only to the first warp)
   const int NWarp        = POT_NTHREAD >> 5;      // optimized for warpSize == 32

   static __shared__ real shared[MaxNWarp];        // maximum shared memory required for 32 partial sums (must be <= warpSize)

// perform reduction within each warp
   val = WarpReductionWithShuffle( val );

// write reduced value to shared memory
   if ( lane == 0 )  shared[wid] = val;

// wait for all partial reductions
   __syncthreads();

// here we have assumed that NWarp < warpSize
   if ( wid == 0 )
   {
//    read from shared memory only if that warp existed
      val = ( ID < NWarp ) ? shared[lane] : (real)0.0;

//    final reduce within first warp
      val = WarpReductionWithShuffle( val );
   }

   return val;

} // FUNCTION : BlockReductionWithShuffle



#endif // #if ( defined GRAVITY  &&  defined GPU )
