#ifndef __CUFLU_FULLSTEPUPDATE__
#define __CUFLU_FULLSTEPUPDATE__


#include <assert.h>
#include "CUFLU.h"

#if (  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

# include "CUFLU_Shared_FluUtility.cu"

#endif // #ifdef __CUDACC__




//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Note        :  1. This function is shared by MHM, MHM_RP, and CTU schemes
//                2. Invoke dual-energy check if DualEnergySwitch is on
//
// Parameter   :  g_Input          : Array storing the input fluid data
//                g_Output         : Array to store the updated fluid data
//                g_DE_Status      : Array to store the dual-energy status
//                g_Flux           : Array storing the input face-centered fluxes
//                                   --> Accessed with the array stride N_FL_FLUX even thought its actually
//                                       allocated size is N_FC_FLUX^3
//                dt               : Time interval to advance solution
//                dh               : Cell size
//                Gamma            : Ratio of specific heats
//                MinDens          : Minimum allowed density
//                MinTemp          : Minimum allowed temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                             const real g_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ], const real dt, const real dh,
                             const real Gamma, const real MinDens, const real MinTemp, int *state )
{

   const int  didx_flux[3] = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh        = dt/dh;

   real dFlux[3][NCOMP_TOTAL], Output_1Cell[NCOMP_TOTAL];


   const int size_ij = SQR(PS2);
   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
      const int i_out    = idx_out % PS2;
      const int j_out    = idx_out % size_ij / PS2;
      const int k_out    = idx_out / size_ij;
      const int idx_flux = IDX321( i_out, j_out, k_out, N_FL_FLUX, N_FL_FLUX );

      const int i_in     = i_out + FLU_GHOST_SIZE;
      const int j_in     = j_out + FLU_GHOST_SIZE;
      const int k_in     = k_out + FLU_GHOST_SIZE;
      const int idx_in   = IDX321( i_in, j_in, k_in, FLU_NXT, FLU_NXT );


//    1. calculate flux difference to update the fluid data
      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         dFlux[d][v] = g_Flux[d][v][ idx_flux+didx_flux[d] ] - g_Flux[d][v][idx_flux];

//    2. full update
      for (int v=0; v<NCOMP_TOTAL; v++)
         Output_1Cell[v] = FMA( dFlux[0][v] + dFlux[1][v] + dFlux[2][v], - dt_dh, g_Input[v][idx_in] );

//    3. check unphysical cell
#     ifdef CHECK_MIN_TEMP
      Output_1Cell[ENGY] = SRHydro_CheckMinTempInEngy( Output_1Cell, MinTemp, Gamma );
#     endif

      if( SRHydro_CheckUnphysical(Output_1Cell, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true) )
      {
#       ifdef __CUDACC__
        atomicOr ( (int*)state, 1);
#       else
        *state = *state | 1;
#       endif
      }

//    waiting all threads within a block
#     ifdef __CUDACC__
      __syncthreads();
#     endif

//    return all threads within a block 
      if ( *state == 1 ) return;


//    4. store results to the output array
      for (int v=0; v<NCOMP_TOTAL; v++)   g_Output[v][idx_out] = Output_1Cell[v];

   } // CGPU_LOOP
  
} // FUNCTION : SRHydro_FullStepUpdate



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )



#endif // #ifndef __CUFLU_FULLSTEPUPDATE__
