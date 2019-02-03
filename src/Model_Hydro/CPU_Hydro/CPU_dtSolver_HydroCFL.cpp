#include "CUFLU.h"

#if ( MODEL == HYDRO )



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

real Hydro_GetPressure( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                        const real Gamma_m1, const bool CheckMinPres, const real MinPres );

#endif // #ifdef __CUDACC__ ... else ...




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
// Parameter   :  g_dt_Array  : Array to store the minimum dt in each target patch
//                g_Flu_Array : Array storing the prepared fluid data of each target patch
//                NPG         : Number of target patch groups (for CPU only)
//                dh          : Cell size
//                Safety      : dt safety factor
//                Gamma       : Ratio of specific heats
//                MinPres     : Minimum allowed pressure
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#else
void CPU_dtSolver_HydroCFL  ( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ], const int NPG,
                              const real dh, const real Safety, const real Gamma, const real MinPres )
#endif
{

   const bool CheckMinPres_Yes = true;
   const real Gamma_m1         = Gamma - (real)1.0;
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
         real fluid[NCOMP_FLUID], _Rho, Vx, Vy, Vz, Pres, Cs, MaxV;

         for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = g_Flu_Array[p][v][t];

        _Rho  = (real)1.0 / fluid[DENS];
         Vx   = FABS( fluid[MOMX] )*_Rho;
         Vy   = FABS( fluid[MOMY] )*_Rho;
         Vz   = FABS( fluid[MOMZ] )*_Rho;
         Pres = Hydro_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                   Gamma_m1, CheckMinPres_Yes, MinPres );
         Cs   = SQRT( Gamma*Pres*_Rho );

#        if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU )
         MaxV   = FMAX( Vx, Vy );
         MaxV   = FMAX( Vz, MaxV );
         MaxCFL = FMAX( MaxV+Cs, MaxCFL );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
         MaxV   = Vx + Vy + Vz;
         MaxCFL = FMAX( MaxV+(real)3.0*Cs, MaxCFL );
#        endif
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
      g_dt_Array[p] = dhSafety/MaxCFL;

   } // for (int p=0; p<8*NPG; p++)

} // FUNCTION : CPU/CUFLU_dtSolver_HydroCFL



#endif // #if ( MODEL == HYDRO )
