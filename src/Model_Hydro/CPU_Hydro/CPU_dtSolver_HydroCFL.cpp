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
// Parameter   :  g_dt_Array             : Array to store the minimum dt in each target patch
//                g_Flu_Array            : Array storing the prepared fluid   data of each target patch
//                g_Mag_Array            : Array storing the prepared B field data of each target patch
//                NPG                    : Number of target patch groups (for CPU only)
//                dh                     : Cell size
//                Safety                 : dt safety factor
//                MinPres                : Minimum allowed pressure
//                EoS_DensEint2Pres_Func : Function pointers to the EoS routines
//                EoS_DensPres2CSqr_Func : ...
//                c_EoS_AuxArray_*       : Auxiliary arrays for the EoS routines (for CPU only)
//                c_EoS_Table            : EoS tables                            (for CPU only)
//                                         --> When using GPU, these CPU-only variables are stored in the constant memory
//                                             header CUDA_ConstMemory.h and do not need to be passed as function arguments
//
// Return      :  g_dt_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                              const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ],
                              const real dh, const real Safety, const real MinPres,
                              const EoS_DE2P_t EoS_DensEint2Pres_Func, const EoS_DP2C_t EoS_DensPres2CSqr_Func )
#else
void CPU_dtSolver_HydroCFL  ( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                              const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const int NPG,
                              const real dh, const real Safety, const real MinPres,
                              const EoS_DE2P_t EoS_DensEint2Pres_Func, const EoS_DP2C_t EoS_DensPres2CSqr_Func,
                              const double c_EoS_AuxArray_Flt[], const int c_EoS_AuxArray_Int[],
                              const real* const c_EoS_Table[EOS_NTABLE_MAX] )
#endif
{

   const bool CheckMinPres_Yes = true;
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
         real fluid[FLU_NIN_T], _Rho, Vx, Vy, Vz, Pres, Emag, a2, CFLx, CFLy, CFLz;
#        ifdef MHD
         int  i, j, k;
         real B[3], Bx2, By2, Bz2, B2, Ca2_plus_a2, Ca2_min_a2, Ca2_min_a2_sqr, four_a2_over_Rho;
#        endif

         for (int v=0; v<FLU_NIN_T; v++)  fluid[v] = g_Flu_Array[p][v][t];

#        ifdef MHD
         i    = t % PS1;
         j    = t % SQR(PS1) / PS1;
         k    = t / SQR(PS1);

         MHD_GetCellCenteredBField( B, g_Mag_Array[p][MAGX], g_Mag_Array[p][MAGY], g_Mag_Array[p][MAGZ], PS1, PS1, PS1, i, j, k );

         Bx2  = SQR( B[MAGX] );
         By2  = SQR( B[MAGY] );
         Bz2  = SQR( B[MAGZ] );
         B2   = Bx2 + By2 + Bz2;
         Emag = (real)0.5*B2;
#        else
         Emag = NULL_REAL;
#        endif

        _Rho   = (real)1.0 / fluid[DENS];
         Vx    = FABS( fluid[MOMX] )*_Rho;
         Vy    = FABS( fluid[MOMY] )*_Rho;
         Vz    = FABS( fluid[MOMZ] )*_Rho;
         Pres  = Hydro_Con2Pres( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], fluid[ENGY], fluid+NCOMP_FLUID,
                                 CheckMinPres_Yes, MinPres, Emag,
                                 EoS_DensEint2Pres_Func, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int, c_EoS_Table, NULL );
         a2    = EoS_DensPres2CSqr_Func( fluid[DENS], Pres, fluid+NCOMP_FLUID, c_EoS_AuxArray_Flt, c_EoS_AuxArray_Int,
                                         c_EoS_Table ); // sound speed squared

//       compute the maximum information propagating speed
//       --> hydro: bulk velocity + sound wave
//           MHD  : bulk velocity +  fast wave
#        ifdef MHD
         Ca2_plus_a2      = B2*_Rho + a2;
         Ca2_min_a2       = B2*_Rho - a2;
         Ca2_min_a2_sqr   = SQR( Ca2_min_a2 );
         four_a2_over_Rho = (real)4.0*a2*_Rho;
         CFLx             = (real)0.5*(  Ca2_plus_a2 + SQRT( Ca2_min_a2_sqr + four_a2_over_Rho*(By2+Bz2) )  );
         CFLy             = (real)0.5*(  Ca2_plus_a2 + SQRT( Ca2_min_a2_sqr + four_a2_over_Rho*(Bx2+Bz2) )  );
         CFLz             = (real)0.5*(  Ca2_plus_a2 + SQRT( Ca2_min_a2_sqr + four_a2_over_Rho*(Bx2+By2) )  );
         CFLx             = SQRT( CFLx );
         CFLy             = SQRT( CFLy );
         CFLz             = SQRT( CFLz );
#        else
         CFLx             = SQRT( a2 );
         CFLy             = CFLx;
         CFLz             = CFLx;
#        endif // #ifdef MHD ... else ...

         CFLx += Vx;
         CFLy += Vy;
         CFLz += Vz;

#        if   ( FLU_SCHEME == RTVD  ||  FLU_SCHEME == CTU )
         MaxCFL = FMAX( CFLx, MaxCFL );
         MaxCFL = FMAX( CFLy, MaxCFL );
         MaxCFL = FMAX( CFLz, MaxCFL );
#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
         MaxCFL = FMAX( CFLx+CFLy+CFLz, MaxCFL );
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
