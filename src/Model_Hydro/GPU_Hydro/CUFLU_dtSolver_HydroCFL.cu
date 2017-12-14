#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO )



#include "CUFLU_Shared_FluUtility.cu"


// parallel reduction routine
#define RED_NTHREAD  DT_FLU_BLOCK_SIZE
#define RED_MAX

#ifdef DT_FLU_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) required for the hydro solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. time-step is estimated by the stability criterion from the von Neumann stability analysis
//
// Parameter   :  g_dt_Array  : Global memory array to store the minimum dt in each target patch
//                g_Flu_Array : Global memory array storing the prepared fluid data of each target patch
//                dh          : Grid size
//                Safety      : dt safety factor
//                Gamma       : Ratio of specific heats
//                MinPres     : Minimum allowed pressure
//
// Return      :  g_dt_Array
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][NCOMP_FLUID][ CUBE(PS1) ],
                                         const real dh, const real Safety, const real Gamma, const real MinPres )
{

   const uint bx               = blockIdx.x;
   const uint ID               = threadIdx.x;
   const real Gamma_m1         = Gamma - (real)1.0;
   const bool CheckMinPres_Yes = true;

   real fluid[NCOMP_FLUID], _Rho, Vx, Vy, Vz, Pres, Cs, MaxV, MaxCFL;
   uint t;


// get the maximum CFL speed evaluated by each thread
   t      = ID;
   MaxCFL = (real)0.0;

   while ( t < CUBE(PS1) )
   {
      for (int v=0; v<NCOMP_FLUID; v++)   fluid[v] = g_Flu_Array[bx][v][t];

     _Rho  = (real)1.0 / fluid[DENS];
      Vx   = FABS( fluid[MOMX] )*_Rho;
      Vy   = FABS( fluid[MOMY] )*_Rho;
      Vz   = FABS( fluid[MOMZ] )*_Rho;
      Pres = CUFLU_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY],
                                Gamma_m1, CheckMinPres_Yes, MinPres );
      Cs   = SQRT( Gamma*Pres*_Rho );

#     if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU  ||  FLU_SCHEME == WAF )
      MaxV   = FMAX( Vx, Vy );
      MaxV   = FMAX( Vz, MaxV );
      MaxCFL = FMAX( MaxV+Cs, MaxCFL );

#     elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
      MaxV   = Vx + Vy + Vz;
      MaxCFL = FMAX( MaxV+(real)3.0*Cs, MaxCFL );
#     endif

      t += DT_FLU_BLOCK_SIZE;
   } // while ( t < CUBE(PS1) )


// perform parallel reduction to get the maximum CFL speed in each thread block
#  ifdef DT_FLU_USE_SHUFFLE
   MaxCFL = BlockReduction_Shuffle ( MaxCFL );
#  else
   MaxCFL = BlockReduction_WarpSync( MaxCFL );
#  endif


// store the minimum dt in each patch back to the global memory
   if ( ID == 0 )    g_dt_Array[bx] = Safety*dh/MaxCFL;

} // FUNCTION : CUFLU_dtSolver_HydroCFL



#endif // #if ( defined GPU  &&  MODEL == HYDRO )
