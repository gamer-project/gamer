#include "GAMER.h"
#include "CUFLU.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp
// Description :  1. Use the corrected coarse-fine boundary fluxes to fix the data at level "lv"
//                2. Use the average data at level "lv+1" to replace the data at level "lv"
//
// Note        :  1. Boundary fluxes from the neighboring ranks must be received in advance by invoking
//                   Buf_GetBufferData()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp( const int lv )
{

   const real Const[6]   = { real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]),
                             real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]),
                             real(-1.0/amr->dh[lv]), real(+1.0/amr->dh[lv]) };
   const int  FluSg      = amr->FluSg[lv];
   const int  Offset[6]  = { 0, PS1-1, 0, (PS1-1)*PS1, 0, (PS1-1)*SQR(PS1) }; // x=0/PS1-1, y=0/PS1-1, z=0/PS1-1 faces
   const int  didx[3][2] = { PS1, SQR(PS1), 1, SQR(PS1), 1, PS1 };

//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
// when enabling cooling, we want to fix the specific internal energy so that it won't be corrected by the flux fix-up operation
// --> it doesn't make much sense to correct the total energy density, pressure, or dual-energy variable when cooling is
//     adopted since the total energy is not conserved anymore
// --> experiments show that fixing the specific internal energy instead of pressure works better
//     --> fixing pressure is found to lead to extremely small dt
//     --> fixing specific internal energy works better since it is directly proportional to the sound speed square
#  ifdef SUPPORT_GRACKLE
   const bool FixSEint   = ( GRACKLE_MODE != GRACKLE_MODE_NONE );
#  else
   const bool FixSEint   = false;
#  endif
*/

   real CorrVal[NFLUX_TOTAL];    // values after applying the flux correction
   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;
   real *FluidPtr1D0[NCOMP_TOTAL], *FluidPtr1D[NCOMP_TOTAL];
   int  didx_m, didx_n;

#  if   ( MODEL == HYDRO  ||  MODEL == MHD )
   const real  Gamma_m1       = GAMMA - (real)1.0;
   const real _Gamma_m1       = (real)1.0 / Gamma_m1;
   const bool CheckMinPres_No = false;
#  ifdef DUAL_ENERGY
   char *DE_StatusPtr1D0, *DE_StatusPtr1D;
#  endif
#  endif // MODEL


// a. correct the coarse-fine boundary fluxes
   if ( OPT__FIXUP_FLUX )
   {
//    check
#     ifdef GAMER_DEBUG

      if ( !amr->WithFlux )
         Aux_Error( ERROR_INFO, "amr->WithFlux is off -> no flux array is allocated for OPT__FIXUP_FLUX !!\n" );

#     if ( MODEL == ELBDM )

#     ifndef CONSERVE_MASS
      Aux_Error( ERROR_INFO, "CONSERVE_MASS is not turned on in the Makefile for the option OPT__FIXUP_FLUX !!\n" );
#     endif

#     if ( NFLUX_TOTAL != 1 )
      Aux_Error( ERROR_INFO, "NFLUX_TOTAL (%d) != 1 for the option OPT__FIXUP_FLUX !!\n", NFLUX_TOTAL );
#     endif

#     if ( DENS != 0 )
      Aux_Error( ERROR_INFO, "DENS (%d) != 0 for the option OPT__FIXUP_FLUX !!\n", DENS );
#     endif

#     if ( FLUX_DENS != 0 )
      Aux_Error( ERROR_INFO, "FLUX_DENS (%d) != 0 for the option OPT__FIXUP_FLUX !!\n", FLUX_DENS );
#     endif


//    if "NCOMP_TOTAL != NFLUX_TOTAL", one must specify how to correct cell data from the flux arrays
//    --> specifically, how to map different flux variables to fluid active/passive variables
//    --> for ELBDM, we have assumed that FLUX_DENS == DENS == 0 and NFLUX_TOTAL == 1
#     elif ( NCOMP_TOTAL != NFLUX_TOTAL )
#        error : NCOMP_TOTAL != NFLUX_TOTAL (one must specify how to map flux variables to fluid active/passive variables) !!

#     endif // #if ( MODEL == ELBDM ) ... elif ...

#     endif // #ifdef GAMER_DEBUG


#     if (  ( MODEL == HYDRO || MODEL == MHD )  &&  defined DUAL_ENERGY  )
#     pragma omp parallel for private( CorrVal, FluxPtr, FluidPtr1D0, FluidPtr1D, didx_m, didx_n, \
                                       DE_StatusPtr1D0, DE_StatusPtr1D ) schedule( runtime )
#     else
#     pragma omp parallel for private( CorrVal, FluxPtr, FluidPtr1D0, FluidPtr1D, didx_m, didx_n ) schedule( runtime )
#     endif
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       a1. sum up the coarse-grid and fine-grid fluxes for bitwise reproducibility
#        ifdef BITWISE_REPRODUCIBILITY
         for (int s=0; s<6; s++)
         {
            FluxPtr = amr->patch[0][lv][PID]->flux[s];

            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] += amr->patch[0][lv][PID]->flux_bitrep[s][v][m][n];
            }
         }
#        endif


//       a2. correct fluid variables by the difference between the coarse-grid and fine-grid fluxes
//       loop over all six faces of a given patch
         for (int s=0; s<6; s++)
         {
//          skip the faces not adjacent to the coarse-fine boundaries
            if ( NULL == (FluxPtr = amr->patch[0][lv][PID]->flux[s]) )  continue;


//          set the pointers to the target face
            for (int v=0; v<NCOMP_TOTAL; v++)
            FluidPtr1D0[v]  = amr->patch[FluSg][lv][PID]->fluid [v][0][0] + Offset[s];
#           ifdef DUAL_ENERGY
            DE_StatusPtr1D0 = amr->patch[    0][lv][PID]->de_status[0][0] + Offset[s];
#           endif


//          set the array index strides
            didx_m = didx[s/2][1];
            didx_n = didx[s/2][0];


//          loop over all cells on a given face
            for (int m=0; m<PS1; m++)
            {
               for (int v=0; v<NCOMP_TOTAL; v++)
               FluidPtr1D[v]  = FluidPtr1D0[v]  + m*didx_m;
#              ifdef DUAL_ENERGY
               DE_StatusPtr1D = DE_StatusPtr1D0 + m*didx_m;
#              endif

               for (int n=0; n<PS1; n++)
               {
//                from now on we also correct cells updated by either the minimum pressure threshold or the 1st-order-flux correction
//                --> note that we have stored the 1st-order fluxes across the coarse-fine boundaries in Flu_Close()
                  /*
#                 ifdef DUAL_ENERGY
                  if ( *DE_StatusPtr1D == DE_UPDATED_BY_MIN_PRES  ||  *DE_StatusPtr1D == DE_UPDATED_BY_1ST_FLUX )    continue;
#                 endif
                  */


//                calculate the corrected results
//                --> do NOT **store** these results yet since we want to skip the cells with unphysical results
                  for (int v=0; v<NFLUX_TOTAL; v++)   CorrVal[v] = *FluidPtr1D[v] + FluxPtr[v][m][n]*Const[s];


//                calculate the pressure
#                 if ( MODEL == HYDRO  ||  MODEL == MHD )
                  real Pres;
                  real *ForPres = CorrVal;

//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
                  real ForPres[NCOMP_TOTAL];
//                when FixSEint is on, use FluidPtr1D to calculate the original pressure
                  if ( FixSEint )
                     for (int v=0; v<NCOMP_TOTAL; v++)   ForPres[v] = *FluidPtr1D[v];
                  else
                     for (int v=0; v<NCOMP_TOTAL; v++)   ForPres[v] = CorrVal    [v];
*/

#                 if   ( DUAL_ENERGY == DE_ENPY )
//                must determine to use CPU_GetPressure() or CPU_DensEntropy2Pres() since the fluid variables stored
//                in CorrVal[] may not be fully consistent (as it's not corrected by CPU_DualEnergyFix())
//                --> note that currently we adopt CPU_DensEntropy2Pres() for DE_UPDATED_BY_MIN_PRES
                  Pres = ( *DE_StatusPtr1D == DE_UPDATED_BY_ETOT  ||  *DE_StatusPtr1D == DE_UPDATED_BY_ETOT_GRA ) ?
                         CPU_GetPressure( ForPres[DENS], ForPres[MOMX], ForPres[MOMY], ForPres[MOMZ], ForPres[ENGY],
                                          Gamma_m1, CheckMinPres_No, NULL_REAL )
                       : CPU_DensEntropy2Pres( ForPres[DENS], ForPres[ENPY], Gamma_m1, CheckMinPres_No, NULL_REAL );

#                 elif ( DUAL_ENERGY == DE_EINT )
#                 error : DE_EINT is NOT supported yet !!

#                 else
                  Pres = CPU_GetPressure( ForPres[DENS], ForPres[MOMX], ForPres[MOMY], ForPres[MOMZ], ForPres[ENGY],
                                          Gamma_m1, CheckMinPres_No, NULL_REAL );
#                 endif // DUAL_ENERGY
#                 endif // MODEL


//###EXPERIMENTAL: (does not work well and thus has been disabled for now)
/*
//                correct pressure to restore the original specific internal energy
                  if ( FixSEint )   Pres *= CorrVal[DENS] / *FluidPtr1D[DENS];
*/


//                do not apply the flux correction if there are any unphysical results
#                 if   ( MODEL == HYDRO  ||  MODEL == MHD )
                  if ( CorrVal[DENS] <= MIN_DENS  ||  Pres <= MIN_PRES  ||  !isfinite(Pres)
#                      if   ( DUAL_ENERGY == DE_ENPY )
                       ||  ( (*DE_StatusPtr1D == DE_UPDATED_BY_DUAL || *DE_StatusPtr1D == DE_UPDATED_BY_MIN_PRES)
                              && CorrVal[ENPY] <= (real)2.0*TINY_NUMBER )

#                      elif ( DUAL_ENERGY == DE_EINT )
#                      error : DE_EINT is NOT supported yet !!
#                      endif
                     )

#                 elif ( MODEL == ELBDM  &&  defined CONSERVE_MASS )
                  if ( CorrVal[DENS] <= MIN_DENS )
#                 endif
                     continue;


//                floor and normalize the passive scalars
#                 if ( NCOMP_PASSIVE > 0 )
                  for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  CorrVal[v] = FMAX( CorrVal[v], TINY_NUMBER );

                  if ( OPT__NORMALIZE_PASSIVE )
                     CPU_NormalizePassive( CorrVal[DENS], CorrVal+NCOMP_FLUID, PassiveNorm_NVar, PassiveNorm_VarIdx );
#                 endif


//                ensure the consistency between pressure, total energy density, and dual-energy variable
//                --> assuming the variable "Pres" is correct
//                --> no need to check the minimum pressure here since we have skipped those cells already
#                 if   ( MODEL == HYDRO  ||  MODEL == MHD )
                  CorrVal[ENGY] = (real)0.5*( SQR(CorrVal[MOMX]) + SQR(CorrVal[MOMY]) + SQR(CorrVal[MOMZ]) ) / CorrVal[DENS]
                                  + Pres*_Gamma_m1;

#                 ifdef DUAL_ENERGY
#                 if   ( DUAL_ENERGY == DE_ENPY )
                  CorrVal[ENPY] = CPU_DensPres2Entropy( CorrVal[DENS], Pres, Gamma_m1 );

#                 elif ( DUAL_ENERGY == DE_EINT )
#                 error : DE_EINT is NOT supported yet !!
#                 endif
#                 endif // #ifdef DUAL_ENERGY
#                 endif // HYDRO/MHD


//                store the corrected results
                  for (int v=0; v<NFLUX_TOTAL; v++)   *FluidPtr1D[v] = CorrVal[v];


//                rescale the real and imaginary parts to be consistent with the corrected amplitude
//                --> must NOT use CorrVal[REAL] and CorrVal[IMAG] below since NFLUX_TOTAL == 1 for ELBDM
#                 if ( MODEL == ELBDM  &&  defined CONSERVE_MASS )
                  real Re, Im, Rho_Corr, Rho_Wrong, Rescale;

                  Re        = *FluidPtr1D[REAL];
                  Im        = *FluidPtr1D[IMAG];
                  Rho_Corr  = *FluidPtr1D[DENS];
                  Rho_Wrong = SQR(Re) + SQR(Im);

//                be careful about the negative density introduced from the round-off errors
                  if ( Rho_Wrong <= (real)0.0  ||  Rho_Corr <= (real)0.0 )
                  {
                     *FluidPtr1D[DENS] = (real)0.0;
                     Rescale           = (real)0.0;
                  }
                  else
                     Rescale = SQRT( Rho_Corr/Rho_Wrong );

                  *FluidPtr1D[REAL] *= Rescale;
                  *FluidPtr1D[IMAG] *= Rescale;
#                 endif


//                update the fluid pointers
                  for (int v=0; v<NCOMP_TOTAL; v++)
                  FluidPtr1D[v]  += didx_n;
#                 ifdef DUAL_ENERGY
                  DE_StatusPtr1D += didx_n;
#                 endif

               } // for (int n=0; n<PS1; n++}
            } // for (int m=0; m<PS1; m++}
         } // for (int s=0; s<6; s++)
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


//    a3. reset all flux arrays (in both real and buffer patches) to zero for bitwise reproducibility
#     ifdef BITWISE_REPRODUCIBILITY
#     pragma omp parallel for private( FluxPtr ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
      {
         for (int s=0; s<6; s++)
         {
            FluxPtr = amr->patch[0][lv][PID]->flux[s];
            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] = 0.0;
            }

            FluxPtr = amr->patch[0][lv][PID]->flux_bitrep[s];
            if ( FluxPtr != NULL )
            {
               for (int v=0; v<NFLUX_TOTAL; v++)
               for (int m=0; m<PS1; m++)
               for (int n=0; n<PS1; n++)
                  FluxPtr[v][m][n] = 0.0;
            }
         }
      }
#     endif

   } // if ( OPT__FIXUP_FLUX )


// b. average over the data at level "lv+1" to correct the data at level "lv"
   if ( OPT__FIXUP_RESTRICT )    Flu_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], NULL_INT, NULL_INT, _TOTAL );

} // FUNCTION : Flu_FixUp
