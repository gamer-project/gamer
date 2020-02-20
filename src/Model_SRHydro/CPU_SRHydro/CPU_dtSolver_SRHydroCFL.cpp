#include "CUFLU.h"

#if ( MODEL == SR_HYDRO )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

// parallel reduction routine
#define RED_NTHREAD  DT_FLU_BLOCK_SIZE
#define RED_MAX

#ifdef DT_FLU_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif

#else // #ifdef __CUDACC__

#  include "../../../include/SRHydroPrototypes.h"

#endif // #ifdef __CUDACC__ ... else ...




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_dtSolver_SRHydroCFL
// Description :  Estimate the evolution time-step (dt) from the CFL condition of the sr-hydro solver
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
//                [3] NPG         : Number of target patch groups (for CPU only)
//                [4] dh          : Cell size
//                [5] Safety      : dt safety factor, defined by (Step==0)?DT__FLUID_INIT:DT__FLUID
//                [6] Gamma       : Ratio of specific heats
//                [7] MinPres     : Minimum allowed pressure
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_SRHydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                              const real dh, const real Safety, const real Gamma, const real MinTemp )
#else
void CPU_dtSolver_SRHydroCFL  ( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const int NPG,
                              const real dh, const real Safety, const real Gamma, const real MinTemp )
#endif
{

   const real dhSafety         = Safety*dh;

// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//                     to work on different patches
#  ifdef __CUDACC__
   const int p = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<8*NPG; p++)
#  endif
   {
      real MaxCFL=(real)0.0;

      CGPU_LOOP( t, CUBE(PS1) )
      {
         real fluid[NCOMP_FLUID], Pri[NCOMP_FLUID], LorentzFactor;
         real Cs, Cs_sq, Us, LorentzFactor_s;
         real Us_Max, LorentzFactor_s_Max;
         real U_Max, LorentzFactor_Max;

         for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = g_Flu_Array[p][v][t];

         LorentzFactor = SRHydro_Con2Pri( fluid, Pri, Gamma, MinTemp );

         Cs_sq = SoundSpeedSquare( Pri[4]/Pri[0], Gamma );
		  
         Cs = SQRT( Cs_sq );

         LorentzFactor_s = (real)1.0 / SQRT( (real)1.0 - Cs_sq );

         Us = Cs * LorentzFactor_s;

#        ifndef USE_3_VELOCITY
		 // 4-velocities -> 3-velocities
		 Pri[1] /= LorentzFactor;
		 Pri[2] /= LorentzFactor;
		 Pri[3] /= LorentzFactor;
#        endif

#        define TIME_TYPE1

//#        if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
#        if ( defined TIME_TYPE1 )
		 U_Max = FABS(Pri[1]) + FABS(Pri[2]) + FABS(Pri[3]);
         Us_Max = (real)3.0 * Us;
#        elif ( defined TIME_TYPE2 )
		 U_Max = FMAX(FABS(Pri[1]), FABS(Pri[2]) );
		 U_Max = FMAX(       U_Max, FABS(Pri[3]) );
         Us_Max = Us;
#        endif

         LorentzFactor_Max   = SQRT( (real)1.0 +  U_Max *  U_Max ); 
         LorentzFactor_s_Max = SQRT( (real)1.0 + Us_Max * Us_Max ); 

         MaxCFL = FMAX( Us_Max * LorentzFactor_Max + LorentzFactor_s_Max * U_Max, MaxCFL );

      } // CGPU_LOOP( t, CUBE(PS1) )

//    perform parallel reduction to get the maximum CFL speed in each thread block
//    --> store in the thread 0
#     ifdef __CUDACC__
#     ifdef DT_FLU_USE_SHUFFLE
      MaxCFL = BlockReduction_Shuffle ( MaxCFL );
#     else
      MaxCFL = BlockReduction_WarpSync( MaxCFL );
#     endif
      if ( threadIdx.x == 0 )
#     endif // #ifdef __CUDACC__
      g_dt_Array[p] = dhSafety / ( MaxCFL / SQRT( (real)1.0 + MaxCFL*MaxCFL ) );

   } // for (int p=0; p<8*NPG; p++)

} // FUNCTION : CPU/CUFLU_dtSolver_SRHydroCFL



#endif // #if ( MODEL == HYDRO )
