#ifndef __CUFLU_FULLSTEPUPDATE__
#define __CUFLU_FULLSTEPUPDATE__



#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



// external functions
#ifdef __CUDACC__

#if ( NCOMP_PASSIVE > 0 )
# include "CUFLU_Shared_FluUtility.cu"
#endif

#ifdef DUAL_ENERGY
# include "CUFLU_Shared_DualEnergy.cu"
#endif

#endif // #ifdef __CUDACC__




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_FullStepUpdate
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
//                MinPres          : Minimum allowed pressure
//                DualEnergySwitch : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive      : true --> normalize passive scalars so that the sum of their mass density
//                                            is equal to the gas mass density
//                NNorm            : Number of passive scalars to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx          : Target variable indices to be normalized
//                                   --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ], const real dt, const real dh,
                           const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                           const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const int  didx_flux[3] = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh        = dt/dh;
#  ifdef DUAL_ENERGY
   const real  Gamma_m1    = Gamma - (real)1.0;
   const real _Gamma_m1    = (real)1.0 / Gamma_m1;
#  endif

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

      for (int v=0; v<NCOMP_TOTAL; v++)
         Output_1Cell[v] = g_Input[v][idx_in] - dt_dh*( dFlux[0][v] + dFlux[1][v] + dFlux[2][v] );


//    we no longer ensure positive density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//        because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
//    --> this consideration holds even when DUAL_ENERGY is adopted (e.g., when density is negative,
//        even when DUAL_ENERGY is on, we still want to try the 1st-order-flux correction before setting a floor value)
      /*
      Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );
      Output_1Cell[ENGY] = Hydro_CheckMinPresInEngy( Output_1Cell[DENS], Output_1Cell[MOMX],
                                                     Output_1Cell[MOMY], Output_1Cell[MOMZ],
                                                     Output_1Cell[ENGY], Gamma_m1, _Gamma_m1, MinPres );
      */


//    2. floor and normalize passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Output_1Cell[v] = FMAX( Output_1Cell[v], TINY_NUMBER );

      if ( NormPassive )
         Hydro_NormalizePassive( Output_1Cell[DENS], Output_1Cell+NCOMP_FLUID, NNorm, NormIdx );
#     endif


//    3. apply the dual-energy formalism to correct the internal energy
//    --> currently, even when UNSPLIT_GRAVITY is on (which would update the internal energy), we still invoke
//        Hydro_DualEnergyFix() here and will fix the internal energy in the gravity solver for cells updated
//        by the dual-energy formalism (i.e., for cells with their dual-energy status marked as DE_UPDATED_BY_DUAL)
//    --> this feature might be modified in the future
#     ifdef DUAL_ENERGY
//    we no longer apply the minimum density and pressure checks here since we want to enable 1st-order-flux correction for that
      const bool CheckMinPres_No = false;
//    Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );

      Hydro_DualEnergyFix( Output_1Cell[DENS], Output_1Cell[MOMX], Output_1Cell[MOMY], Output_1Cell[MOMZ],
                           Output_1Cell[ENGY], Output_1Cell[ENPY], g_DE_Status[idx_out],
                           Gamma_m1, _Gamma_m1, CheckMinPres_No, NULL_REAL, DualEnergySwitch );
#     endif // #ifdef DUAL_ENERGY


//    4. store results to the output array
      for (int v=0; v<NCOMP_TOTAL; v++)   g_Output[v][idx_out] = Output_1Cell[v];


//    5. check the negative density and energy
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(Output_1Cell[DENS]) )
         printf( "WARNING : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Output_1Cell[DENS], __FILE__, __LINE__, __FUNCTION__ );

      if ( Hydro_CheckNegative(Output_1Cell[ENGY]) )
         printf( "WARNING : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Output_1Cell[ENGY], __FILE__, __LINE__, __FUNCTION__ );
#     endif

   } // i,j,k

} // FUNCTION : Hydro_FullStepUpdate



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )



#endif // #ifndef __CUFLU_FULLSTEPUPDATE__
