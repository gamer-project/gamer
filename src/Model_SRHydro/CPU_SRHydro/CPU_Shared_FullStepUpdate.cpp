#include "GAMER.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == SR_HYDRO  && ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )

bool CPU_CheckUnphysical( const real Con[], const real Pri[], const char s[], const int line);
real CPU_CheckMinTempInEngy (const real Con[]);
static bool boolean;
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Parameter   :  [ 1] Input            : Array storing the input initial data
//                [ 2] Output           : Array to store the ouptut updated data
//                [ 3] DE_Status        : Array to store the dual-energy status
//                [ 4] Flux             : Array storing the input face-centered flux
//                                        --> Size is assumed to be N_FL_FLUX^3
//                [ 5] dt               : Time interval to advance solution
//                [ 6] dh               : Grid size
//                [ 7] Gamma            : Ratio of specific heats
//                [ 8] MinDens          : Minimum allowed density
//                [ 9] MinPres          : Minimum allowed pressure
//                [10] DualEnergySwitch : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                [11] NormPassive      : true --> normalize passive scalars so that the sum of their mass density
//                                                 is equal to the gas mass density
//                [12] NNorm            : Number of passive scalars to be normalized
//                                        --> Should be set to the global variable "PassiveNorm_NVar"
//                [13] NormIdx          : Target variable indices to be normalized
//                                        --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                         const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                         const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                         const bool NormPassive, const int NNorm, const int NormIdx[] )
{

#  ifdef DUAL_ENERGY
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif
   const int  dID1[3]   = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh     = dt/dh;

   int  ID1, ID2, ID3;
   real dF[3][NCOMP_TOTAL];

#  if ( NCOMP_PASSIVE > 0 )
   real Passive[NCOMP_PASSIVE];
#  endif


   for (int k1=0, k2=FLU_GHOST_SIZE;  k1<PS2;  k1++, k2++)
   for (int j1=0, j2=FLU_GHOST_SIZE;  j1<PS2;  j1++, j2++)
   for (int i1=0, i2=FLU_GHOST_SIZE;  i1<PS2;  i1++, i2++)
   {

      ID1 = (k1*N_FL_FLUX + j1)*N_FL_FLUX + i1;
      ID2 = (k1*PS2       + j1)*PS2       + i1;
      ID3 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

#     ifdef CHECK_NEGATIVE_IN_FLUID
      real Con[NCOMP_FLUID];
      for (int v = 0;v<NCOMP_FLUID;v++) Con[v] = Input[v][ID3];
      boolean = CPU_CheckUnphysical(Con, NULL, __FUNCTION__, __LINE__);
#     endif

      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)   dF[d][v] = Flux[ ID1+dID1[d] ][d][v] - Flux[ID1][d][v];

      for (int v=0; v<NCOMP_TOTAL; v++)
         Output[v][ID2] = Input[v][ID3] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );


//    we no longer ensure positive density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//        because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
//    --> this consideration holds even when DUAL_ENERGY is adopted (e.g., when density is negative, even when DUAL_ENERGY is on,
//        we still want to try the 1st-order-flux correction before setting a floor value)
/*      
      Output[DENS][ID2] = FMAX( Output[DENS][ID2], MinDens );

      Output[ENGY][ID2] = CPU_CheckMinPresInEngy( Output[DENS][ID2], Output[MOMX][ID2], Output[MOMY][ID2], Output[MOMZ][ID2],
                                                  Output[ENGY][ID2], Gamma_m1, _Gamma_m1, MinPres );
*/
      

//    floor and normalize passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Output[v][ID2] = FMAX( Output[v][ID2], TINY_NUMBER );

      if ( NormPassive )
      {
         for (int v=0; v<NCOMP_PASSIVE; v++)    Passive[v] = Output[ NCOMP_FLUID + v ][ID2];

         CPU_NormalizePassive( Output[DENS][ID2], Passive, NNorm, NormIdx );

         for (int v=0; v<NCOMP_PASSIVE; v++)    Output[ NCOMP_FLUID + v ][ID2] = Passive[v];
      }
#     endif


//    apply the dual-energy formalism to correct the internal energy
//    --> currently, even when UNSPLIT_GRAVITY is on (which would update the internal energy), we still invoke
//        CPU_DualEnergyFix() here and will fix the internal energy in the gravity solver for cells updated
//        by the dual-energy formalism (i.e., for cells with their dual-energy status marked as DE_UPDATED_BY_DUAL)
//    --> this feature might be modified in the future
#     ifdef DUAL_ENERGY
//    we no longer apply the minimum density and pressure checks here since we want to enable 1st-order-flux correction for that
      const bool CheckMinPres_No = false;
//    Output[DENS][ID2] = FMAX( Output[DENS][ID2], MinDens );

      CPU_DualEnergyFix( Output[DENS][ID2], Output[MOMX][ID2], Output[MOMY][ID2], Output[MOMZ][ID2],
                         Output[ENGY][ID2], Output[ENPY][ID2], DE_Status[ID2],
                         Gamma_m1, _Gamma_m1, CheckMinPres_No, NULL_REAL, DualEnergySwitch );
#     endif // #ifdef DUAL_ENERGY

      real Cons[NCOMP_FLUID];
      for (int v = 0;v<NCOMP_FLUID;v++) Cons[v] = Output[v][ID2];
#     ifdef CHECK_MIN_TEMP
      Output[ENGY][ID2] = CPU_CheckMinTempInEngy( Cons );
      Cons[ENGY] = Output[ENGY][ID2];
#     endif

#     ifdef CHECK_NEGATIVE_IN_FLUID
      boolean = CPU_CheckUnphysical(Cons, NULL, __FUNCTION__, __LINE__);
#     endif

   } // i,j,k

} // FUNCTION : CPU_FullStepUpdate



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_StoreFlux
// Description :  Store the inter-patch fluxes for the AMR fix-up operation
//
// Parameter   :  Output   : Array to store the inter-patch fluxes
//                FC_Flux  : Array storing the face-centered fluxes
//                           --> Size is assumed to be N_FL_FLUX^3
//-------------------------------------------------------------------------------------------------------
void CPU_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL]  )
{

   int Face, ID1, ID2[9];

   for (int m=0; m<PS2; m++)
   for (int n=0; n<PS2; n++)
   {
      ID1    = m*PS2 + n;

      ID2[0] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + 0;
      ID2[1] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS1;
      ID2[2] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS2;
      ID2[3] = (  m*N_FL_FLUX +   0)*N_FL_FLUX + n;
      ID2[4] = (  m*N_FL_FLUX + PS1)*N_FL_FLUX + n;
      ID2[5] = (  m*N_FL_FLUX + PS2)*N_FL_FLUX + n;
      ID2[6] = (  0*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[7] = (PS1*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[8] = (PS2*N_FL_FLUX +   m)*N_FL_FLUX + n;

      for (int t=0; t<9; t++)
      {
         Face = t/3;

         for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Array[t][v][ID1] = FC_Flux[ ID2[t] ][Face][v];
      }
   }

} // FUNCTION : CPU_StoreFlux



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

