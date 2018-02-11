#ifdef __CUDACC__
#include "Macro.h"
#else
#include "GAMER.h"
#endif
#include "CUFLU.h"

// some functions in this file need to be defined even when using GPU
#if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )

#if ( DUAL_ENERGY == DE_ENPY  &&  defined __CUDACC__ )
static __device__ real CUFLU_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 );
static __device__ real CUFLU_DensEntropy2Pres( const real Dens, const real Enpy, const real Gamma_m1,
                                               const bool CheckMinPres, const real MinPres );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_DualEnergyFix
// Description :  Correct the internal and total energies using the dual-energy formalism
//
// Note        :  1. Invoked by the functions "CPU_FullStepUpdate, CUFLU_FullStepUpdate"
//                   --> This function is invoked by both CPU and GPU codes
//                2. A floor value "MinPres" is applied to the corrected pressure if CheckMinPres is on
//                3  A floor value "TINY_NUMBER" is applied to the inputted entropy as well
//                4. Call-by-reference for "Etot, Enpy, and DE_Status"
//                5. This function always ensures that the returned fluid variables are consistent with each other
//                   --> They must satisfy "entropy = pressure / density^(Gamma-1)", where pressure is calculated
//                       by (Etot - Ekin)*(Gamma-1.0)
//                   --> It doesn't matter we use entropy to correct Eint or vice versa, and it also holds even when
//                       the floor value is applied to pressure
//
// Parameter   :  Dens             : Mass density
//                MomX/Y/Z         : Momentum density
//                Etot             : Total energy density
//                Enpy             : Entropy
//                DE_Status        : Assigned to (DE_UPDATED_BY_ETOT / DE_UPDATED_BY_DUAL / DE_UPDATED_BY_MIN_PRES)
//                                   to indicate whether this cell is updated by the total energy, dual energy variable,
//                                   or minimum allowed pressure (MinPres)
//                Gamma_m1         : Adiabatic index - 1.0
//                _Gamma_m1        : 1.0/Gamma_m1
//                CheckMinPres     : Return CPU/CUFLU_CheckMinPres()
//                                   --> In some cases we actually want to check if pressure becomes unphysical,
//                                       for which we don't want to enable this option
//                MinPres          : Minimum allowed pressure
//                DualEnergySwitch : if ( Eint/Ekin < DualEnergySwitch ) ==> correct Eint and Etot
//                                   else                                ==> correct Enpy
//
// Return      :  Etot, Enpy, DE_Status
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
void CUFLU_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Enpy, char &DE_Status, const real Gamma_m1, const real _Gamma_m1,
                          const bool CheckMinPres, const real MinPres, const real DualEnergySwitch )
#else
void   CPU_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Enpy, char &DE_Status, const real Gamma_m1, const real _Gamma_m1,
                          const bool CheckMinPres, const real MinPres, const real DualEnergySwitch )
#endif
{

   const bool CheckMinPres_No = false;

// apply the minimum entropy check
   Enpy = FMAX( Enpy, TINY_NUMBER );


// calculate energies --> note that here Eint can even be negative due to numerical errors
   real Ekin, Eint, Pres;

   Ekin = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) )/Dens;
   Eint = Etot - Ekin;


// determine whether or not to use the dual-energy variable (entropy or internal energy) to correct the total energy density
   if ( Eint/Ekin < DualEnergySwitch )
   {
//    correct total energy
//    --> we will check the minimum pressure later
#     if   ( DUAL_ENERGY == DE_ENPY )
#     ifdef __CUDACC__
      Pres = CUFLU_DensEntropy2Pres( Dens, Enpy, Gamma_m1, CheckMinPres_No, NULL_REAL );
#     else
      Pres =   CPU_DensEntropy2Pres( Dens, Enpy, Gamma_m1, CheckMinPres_No, NULL_REAL );
#     endif
      Eint = Pres*_Gamma_m1;

#     elif ( DUAL_ENERGY == DE_EINT )
#     error : DE_EINT is NOT supported yet !!
#     endif

      Etot      = Ekin + Eint;
      DE_Status = DE_UPDATED_BY_DUAL;
   }

   else
   {
//    correct entropy
      Pres      = Eint*Gamma_m1;
#     ifdef __CUDACC__
      Enpy      = CUFLU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
#     else
      Enpy      =   CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
#     endif
      DE_Status = DE_UPDATED_BY_ETOT;
   } // if ( Eint/Ekin < DualEnergySwitch ) ... else ...


// apply the minimum pressure check
   if ( CheckMinPres  &&  Pres < MinPres )
   {
      Pres = MinPres;
      Eint = Pres*_Gamma_m1;

//    ensure that both energy and entropy are consistent with the pressure floor
      Etot      = Ekin + Eint;
#     ifdef __CUDACC__
      Enpy      = CUFLU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
#     else
      Enpy      =   CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
#     endif
      DE_Status = DE_UPDATED_BY_MIN_PRES;
   }

} // FUNCTION : CPU/CUFLU_DualEnergyFix



#if ( DUAL_ENERGY == DE_ENPY )

// CPU_Fluid2Entropy() does not need to have the GPU version
#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Fluid2Entropy
// Description :  Evaluate the gas entropy from the input fluid variables
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "Hydro_Init_ByFunction_AssignData, Gra_Close"
//                3. Currently this function does NOT apply the minimum pressure check when calling CPU_GetPressure()
//                   --> However, note that CPU_DensPres2Entropy() does apply a floor value (TINY_NUMBER) for entropy
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Engy     : Energy density
//                Gamma_m1 : Adiabatic index - 1.0
//
// Return      :  Enpy
//-------------------------------------------------------------------------------------------------------
real CPU_Fluid2Entropy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy, const real Gamma_m1 )
{

// currently this function does NOT apply the minimum pressure check when calling CPU_GetPressure()
   const bool CheckMinPres_No = false;

   real Pres, Enpy;

// calculate pressure and convert it to entropy
   Pres = CPU_GetPressure( Dens, MomX, MomY, MomZ, Engy, Gamma_m1, CheckMinPres_No, NULL_REAL );
   Enpy = CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );

   return Enpy;

} // FUNCTION : CPU_Fluid2Entropy
#endif // ifndef __CUDACC__



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_DensPres2Entropy
// Description :  Evaluate the gas entropy from the input density and pressure
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "CPU_Fluid2Entropy, CPU_DualEnergyFix, CUFLU_DualEnergyFix"
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value (TINY_NUMBER) is applied to the returned value
//
// Parameter   :  Dens     : Mass density
//                Pres     : Pressure
//                Gamma_m1 : Adiabatic index - 1.0
//
// Return      :  Enpy
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
real CUFLU_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 )
#else
real   CPU_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 )
#endif
{

   real Enpy;

// calculate entropy
   Enpy = Pres*POW( Dens, -Gamma_m1 );

// apply a floor value
   Enpy = FMAX( Enpy, TINY_NUMBER );

   return Enpy;

} // FUNCTION : CPU/CUFLU_DensPres2Entropy



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_DensEntropy2Pres
// Description :  Evaluate the gas pressure from the input density and entropy
//                --> Here entropy is defined as "pressure / density^(Gamma-1)" (i.e., entropy per volume)
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "CPU_DualEnergyFix, Flu_Close, Hydro_Aux_Check_Negative,
//                   Flu_FixUp, CUFLU_DualEnergyFix"
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value "MinPres" is applied to the returned pressure if CheckMinPres is on
//
// Parameter   :  Dens         : Mass density
//                Enpy         : Enpy
//                Gamma_m1     : Adiabatic index - 1.0
//                CheckMinPres : Return CPU/CUFLU_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//
// Return      :  Pres
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__forceinline__ __device__
real CUFLU_DensEntropy2Pres( const real Dens, const real Enpy, const real Gamma_m1,
                             const bool CheckMinPres, const real MinPres )
#else
real   CPU_DensEntropy2Pres( const real Dens, const real Enpy, const real Gamma_m1,
                             const bool CheckMinPres, const real MinPres )
#endif
{

   real Pres;

// calculate pressure
   Pres = Enpy*POW( Dens, Gamma_m1 );

// apply a floor value
   if ( CheckMinPres )
   {
#     ifdef __CUDACC__
      Pres = CUFLU_CheckMinPres( Pres, MinPres );
#     else
      Pres =   CPU_CheckMinPres( Pres, MinPres );
#     endif
   }

   return Pres;

} // FUNCTION : CPU/CUFLU_DensEntropy2Pres

#endif // #if ( DUAL_ENERGY == DE_ENPY )



#endif // #if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
