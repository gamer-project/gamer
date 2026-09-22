#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

// parallel reduction routine
#define RED_NTHREAD  DT_FLU_BLOCK_SIZE
#define RED_MAX

#ifdef DT_FLU_USE_SHUFFLE
#  include "../../GPU_Utility/CUUTI_BlockReduction_Shuffle.cu"
#else
#  include "../../GPU_Utility/CUUTI_BlockReduction_WarpSync.cu"
#endif

#endif // #ifdef __CUDACC__




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_dtSolver_HydroCFL
// Description :  Estimate the evolution time-step (dt) from the CFL condition of the hydro/MHD solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Time-step is estimated by the stability criterion from the von Neumann stability analysis
//                   --> CFL condition
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_dt_Array   : Array to store the minimum dt in each target patch
//                g_Flu_Array  : Array storing the prepared fluid   data of each target patch
//                g_Mag_Array  : Array storing the prepared B field data of each target patch
//                NPG          : Number of target patch groups (for CPU only)
//                dh           : Cell size
//                Safety       : dt safety factor
//                MinPres      : Minimum allowed pressure
//                PassiveFloor : Bitwise flag to specify the passive scalars to be floored
//                EoS          : EoS object
//                MicroPhy     : Microphysics object
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                              const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ],
                              const real dh, const real Safety, const real MinPres,
                              const long PassiveFloor, const EoS_t EoS, const MicroPhy_t MicroPhy )
#else
void CPU_dtSolver_HydroCFL  ( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                              const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const int NPG,
                              const real dh, const real Safety, const real MinPres,
                              const long PassiveFloor, const EoS_t EoS, const MicroPhy_t MicroPhy )
#endif
{

   const real dhSafety         = Safety*dh;
#  ifdef CR_DIFFUSION
   const real dh2Safety        = MicroPhy.CR_safety*0.5*dh*dh;
#  endif

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
         real fluid[FLU_NIN_T];

         for (int v=0; v<FLU_NIN_T; v++)  fluid[v] = g_Flu_Array[p][v][t];

         real B[3] = { (real)0.0, (real)0.0, (real)0.0 };
#        ifdef MHD
         int  i, j, k;

         i    = t % PS1;
         j    = t % SQR(PS1) / PS1;
         k    = t / SQR(PS1);

         MHD_GetCellCenteredBField( B, g_Mag_Array[p][MAGX], g_Mag_Array[p][MAGY], g_Mag_Array[p][MAGZ], PS1, PS1, PS1, i, j, k );
#        endif // #ifdef MHD

         MaxCFL = FMAX( MaxCFL,
                        Hydro_GetCFL( fluid, B, MinPres, PassiveFloor,
                                      EoS.DensEint2Pres_FuncPtr, EoS.DensPres2Eint_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                      EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                                      EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table ) );

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

#     ifdef SRHD
      g_dt_Array[p] = dhSafety / ( MaxCFL / SQRT( (real)1.0 + MaxCFL*MaxCFL ) );
#     else
      g_dt_Array[p] = dhSafety/MaxCFL;
#     endif

//    The CFL condition determined by the cosmic ray diffusion
#     ifdef CR_DIFFUSION
      real MaxCFL_CRDiffusion=(real)0.0;

      CGPU_LOOP( t, CUBE(PS1) )
      {
        MaxCFL_CRDiffusion = FMAX( MaxCFL_CRDiffusion, Hydro_GetCFL_CRDiffusion( MicroPhy ) );

      } // CGPU_LOOP( t, CUBE(PS1) )

#     ifdef __CUDACC__
#     ifdef DT_FLU_USE_SHUFFLE
      MaxCFL_CRDiffusion = BlockReduction_Shuffle ( MaxCFL_CRDiffusion );
#     else
      MaxCFL_CRDiffusion = BlockReduction_WarpSync( MaxCFL_CRDiffusion );
#     endif
      if ( threadIdx.x == 0 )
#     endif // #ifdef __CUDACC__
      g_dt_Array[p] = FMIN( dh2Safety/MaxCFL_CRDiffusion, g_dt_Array[p] );

#     endif // #ifdef CR_DIFFUSION

   } // for (int p=0; p<8*NPG; p++)

} // FUNCTION : CPU/CUFLU_dtSolver_HydroCFL



#endif // #if ( MODEL == HYDRO )
