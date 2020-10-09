#include "CUPOT.h"

#ifdef GRAVITY



// include c_ExtPot_AuxArray_Flt/Int[]
#ifdef __CUDACC__
#  include "CUDA_ConstMemory.h"
#endif




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_ExtPotSolver
// Description :  Add external potential
//
// Note        :  1. External potential is specified by the input function ExtPot_Func()
//                2. Set PotIsInit to false if the input potential g_Pot_Array[] has not been initialized
//                   --> Useful when self-gravity is disabled
//                3. Invoked by Gra_AdvanceDt(), CPU_PoissonGravitySolver(), and CUAPI_Asyn_PoissonGravitySolver()
//
// Parameter   :  g_Pot_Array               : Array storing the input and output potential data of each target patch
//                g_Corner_Array            : Array storing the physical corner coordinates of each patch
//                g_ExtPotTable             : Array storing the external potential 3D table
//                NPatchGroup               : Number of target patch groups (for CPU only)
//                dh                        : Cell size
//                ExtPot_Func               : Function pointer to the external potential routine (for both CPU and GPU)
//                c_ExtPot_AuxArray_Flt/Int : Auxiliary floating-point/integer arrays for adding external potential (for CPU only)
//                                            --> When using GPU, this array is stored in the constant memory header
//                                                CUDA_ConstMemory.h and does not need to be passed as a function argument
//                Time                      : Target physical time
//                PotIsInit                 : Whether the input potential has been initialized
//                                            --> true : **add** external potential to the input data
//                                                false: **overwrite** the input data
//
// Return      :  g_Pot_Array[]
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_ExtPotSolver( real g_Pot_Array[][ CUBE(GRA_NXT) ],
                         const double g_Corner_Array[][3],
                         const real g_ExtPotTable[],
                         const real dh, const ExtPot_t ExtPot_Func,
                         const double Time, const bool PotIsInit )
#else
void CPU_ExtPotSolver  ( real g_Pot_Array[][ CUBE(GRA_NXT) ],
                         const double g_Corner_Array[][3],
                         const real g_ExtPotTable[],
                         const int NPatchGroup,
                         const real dh, const ExtPot_t ExtPot_Func,
                         const double c_ExtPot_AuxArray_Flt[],
                         const int    c_ExtPot_AuxArray_Int[],
                         const double Time, const bool PotIsInit )
#endif
{

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
      const double x0 = g_Corner_Array[P][0] - GRA_GHOST_SIZE*dh;
      const double y0 = g_Corner_Array[P][1] - GRA_GHOST_SIZE*dh;
      const double z0 = g_Corner_Array[P][2] - GRA_GHOST_SIZE*dh;

      double x, y, z;
      real   ExtPot;

//    loop over all cells of the target patch
      CGPU_LOOP( t, CUBE(GRA_NXT) )
      {
         const int i = t % GRA_NXT;
         const int j = t % SQR(GRA_NXT) / GRA_NXT;
         const int k = t / SQR(GRA_NXT);

         x = x0 + double(i*dh);
         y = y0 + double(j*dh);
         z = z0 + double(k*dh);

         ExtPot = ExtPot_Func( x, y, z, Time, c_ExtPot_AuxArray_Flt, c_ExtPot_AuxArray_Int, EXT_POT_USAGE_ADD, g_ExtPotTable );

         if ( PotIsInit )  g_Pot_Array[P][t] += ExtPot;  // add to the input potential
         else              g_Pot_Array[P][t]  = ExtPot;  // overwrite the input potential
      }
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_ExtPotSolver



#endif // #ifdef GRAVITY
