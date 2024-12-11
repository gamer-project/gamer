// include <mpi.h> BEFORE "Macro.h" to avoid conflicting definition of the symbolic constant "REAL"
// --> for example, it's also defined in "openmpi-1.4.3-qlc/include/openmpi/ompi/mpi/cxx/constants.h(119)"
#if ( !defined __CUDACC__  &&  !defined SERIAL )
#include <mpi.h>
#endif

#include "CUPOT.h"

#if ( MODEL == ELBDM  &&  defined GRAVITY  &&  ELBDM_SCHEME == ELBDM_HYBRID )




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_ELBDMGravitySolver_HamiltonJacobi
// Description :  CPU/GPU ELBDM gravity solver for advancing the phase S by S = S - Eta*(Phi + Lambda*Rho)*dt
//
// Note        :  1. ELBDM gravity solver requires NO potential and fluid ghost zone
//                   --> Optimized performance can be achieved if GRA_GHOST_SIZE==0 (and thus GRA_NXT==PATCH_SIZE)
//                   --> But the code supports GRA_GHOST_SIZE>0 as well
//                       --> Mainly for the STORE_POT_GHOST option
//                2. ELBDM gravity solver does only need the density information if QUARTIC_SELF_INTERACTION is on
//                   --> GRA_NIN==2 (only store the density and phase)
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//                4. No GPU shared memory is used in this kernel since no computational stencil is required
//                   and hence no data needed to be shared
//
// Parameter   :  g_Flu_Array : Array to store the input and output data
//                g_Pot_Array : Array storing the input potential for evaluating the gravitational acceleration
//                NPatchGroup : Number of patch groups to be evaluated (for CPU only)
//                EtaDt       : Particle mass / Planck constant * dt
//                dh          : Cell size
//                Lambda      : Quartic self-interaction coefficient in ELBDM
//
// Return      :  g_Flu_Array[]
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_ELBDMGravitySolver_HamiltonJacobi(       real g_Flu_Array[][GRA_NIN][ CUBE(PS1) ],
                                              const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                              const real EtaDt, const real dh, const real Lambda )
#else
void CPU_ELBDMGravitySolver_HamiltonJacobi  (       real g_Flu_Array[][GRA_NIN][ CUBE(PS1) ],
                                              const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                              const int NPatchGroup,
                                              const real EtaDt, const real dh, const real Lambda )
#endif
{

   const int PS1_sqr = SQR(PS1);

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
//    loop over all cells of the target patch
      CGPU_LOOP( idx_flu, CUBE(PS1) )
      {
         const int i_flu   = idx_flu % PS1;
         const int j_flu   = idx_flu % PS1_sqr / PS1;
         const int k_flu   = idx_flu / PS1_sqr;

         const int i_pot   = i_flu + GRA_GHOST_SIZE;
         const int j_pot   = j_flu + GRA_GHOST_SIZE;
         const int k_pot   = k_flu + GRA_GHOST_SIZE;
         const int idx_pot = IDX321( i_pot, j_pot, k_pot, GRA_NXT, GRA_NXT );

         real Pot = g_Pot_Array[P][idx_pot];

#        ifdef QUARTIC_SELF_INTERACTION
         const real Rho = g_Flu_Array[P][0][idx_flu];

         Pot += Lambda * Rho;
#        endif

         g_Flu_Array[P][PHAS][idx_flu] -= EtaDt * Pot;

      } // CGPU_LOOP( idx_flu, CUBE(PS1) )
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_ELBDMGravitySolver_HJ



#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  ELBDM_SCHEME == ELBDM_HYBRID )
