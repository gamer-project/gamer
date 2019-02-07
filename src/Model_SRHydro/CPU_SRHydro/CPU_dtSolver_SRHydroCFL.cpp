#include "CUFLU.h"


#if ( MODEL == SR_HYDRO )

//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) from the CFL condition of the hydro solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//                   --> CFL condition
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  [1] g_dt_Array  : Array to store the minimum dt in each target patch
//                [2] g_Flu_Array : Array storing the prepared fluid data of each target patch
//                [3] NPG         : Number of target patch groups
//                [4] dh          : Cell size
//                [5] Safety      : dt safety factor
//                [6] Gamma       : Ratio of specific heats
//                [7] MinPres     : Minimum allowed pressure
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_SRHydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#else
void CPU_dtSolver_SRHydroCFL  ( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const int NPG,
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#endif
{
   const real dhSafety = Safety*dh;

// loop over all patches
#  ifdef __CUDACC__
   const int p = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<8*NPatch; p++)
#  endif
   {
     g_dt_Array[p] = dhSafety;
   }

} // FUNCTION : CPU_dtSolver_HydroCFL



#endif // #if ( MODEL == SR_HYDRO )
