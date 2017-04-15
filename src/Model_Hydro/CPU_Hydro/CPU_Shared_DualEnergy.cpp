#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

// some functions in this file need to be defined even when using GPU
#if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )

#if ( DUAL_ENERGY == DE_ENTROPY )
real CPU_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DualEnergyFix
// Description :  Correct the internal and total energies using the dual-energy formalism
//
// Note        :  1. Invoked by the functions "CPU_FullStepUpdate, ..."
//                2. This function checks the minimum pressure since currently DUAL_ENERGY does NOT work
//                   with OPT__1ST_FLUX_CORR
//                3. Call-by-reference for "Etot and Entropy"
//                4. This function always ensures that the returned fluid variables are consistent with each other
//                   --> They must satisfy "entropy = pressure / density^(Gamma-1)", where pressure is calculated
//                       by (Etot - Ekin)*(Gamma-1.0)
//                   --> It doesn't matter we use entropy to correct Eint or vice versa, and it also holds even when
//                       the floor value is applied to pressure
//
// Parameter   :  Dens             : Mass density
//                MomX/Y/Z         : Momentum density
//                Etot             : Total energy density
//                Entropy          : Entropy
//                Gamma_m1         : Adiabatic index - 1.0
//                _Gamma_m1        : 1.0/Gamma_m1
//                MinPres          : Minimum allowed pressure
//                DualEnergySwitch : if ( Eint/Ekin < DualEnergySwitch ) ==> correct Eint and Etot
//                                   else                                ==> correct Entropy
//
// Return      :  Etot, Entropy
//-------------------------------------------------------------------------------------------------------
void CPU_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ, real &Etot, real &Entropy,
                        const real Gamma_m1, const real _Gamma_m1, const real MinPres, const real DualEnergySwitch )
{

// calculate energies --> note that here Eint can even be negative due to numerical errors
   real Ekin, Eint, Pres;

   Ekin = (real)0.5*( SQR(MomX) + SQR(MomY) + SQR(MomZ) )/Dens;
   Eint = Etot - Ekin;


// determine whether or not to use the dual-energy variable (entropy or internal energy) to correct the total energy density
   if ( Eint/Ekin < DualEnergySwitch )
   {
//    correct total energy
#     if   ( DUAL_ENERGY == DE_ENTROPY )
      Pres = CPU_DensEntropy2Pres( Dens, Entropy, Gamma_m1, MinPres );
      Eint = Pres*_Gamma_m1;
#     elif ( DUAL_ENERGY == DE_EINT )
#     error : DE_EINT is NOT supported yet !!
#     endif

      Etot = Ekin + Eint;
   }

   else
   {
//    correct entropy
      Pres    = Eint*Gamma_m1;
      Entropy = CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
   } // if ( Eint/Ekin < DualEnergySwitch ) ... else ...


// apply the minimum pressure check
// --> must include "=" since we already apply MinPres check when calling CPU_DensEntropy2Pres() above
   if ( Pres <= MinPres )
   {
      Pres = MinPres;
      Eint = Pres*_Gamma_m1;

//    ensure that both energy and entropy are consistent with the pressure floor
      Etot    = Ekin + Eint;
      Entropy = CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );
   }

} // FUNCTION : CPU_DualEnergyFix



#if ( DUAL_ENERGY == DE_ENTROPY )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_Fluid2Entropy
// Description :  Evaluate the gas entropy from the input fluid variables
//                --> Here entropy is defined as "pressure / density^(Gamma-1)"
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "Hydro_Init_StartOver_AssignData, ..."
//                3. Currently this function does NOT apply the minimum pressure check when calling CPU_GetPressure()
//                   --> However, note that CPU_DensPres2Entropy() does apply a floor value (TINY_VALUE) for entropy
//
// Parameter   :  Dens     : Mass density
//                MomX/Y/Z : Momentum density
//                Engy     : Energy density
//                Gamma_m1 : Adiabatic index - 1.0
//
// Return      :  Entropy
//-------------------------------------------------------------------------------------------------------
real CPU_Fluid2Entropy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy, const real Gamma_m1 )
{

// currently this function does NOT apply the minimum pressure check when calling CPU_GetPressure()
   const bool CheckMinPres_No = false;

   real Pres, Entropy;

// calculate pressure and convert it to entropy
   Pres    = CPU_GetPressure( Dens, MomX, MomY, MomZ, Engy, Gamma_m1, CheckMinPres_No, NULL_REAL );
   Entropy = CPU_DensPres2Entropy( Dens, Pres, Gamma_m1 );

   return Entropy;

} // FUNCTION : CPU_Fluid2Entropy



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DensPres2Entropy
// Description :  Evaluate the gas entropy from the input density and pressure
//                --> Here entropy is defined as "pressure / density^(Gamma-1)"
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "CPU_Fluid2Entropy, ..."
//                3. A floor value (TINY_VALUE) is applied to the returned value
//
// Parameter   :  Dens     : Mass density
//                Pres     : Pressure
//                Gamma_m1 : Adiabatic index - 1.0
//
// Return      :  Entropy
//-------------------------------------------------------------------------------------------------------
real CPU_DensPres2Entropy( const real Dens, const real Pres, const real Gamma_m1 )
{

   real Entropy;

// calculate entropy
   Entropy = Pres*POW( Dens, -Gamma_m1 );

// apply a floor value
   Entropy = FMAX( Entropy, TINY_NUMBER );

   return Entropy;

} // FUNCTION : CPU_DensPres2Entropy



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DensEntropy2Pres
// Description :  Evaluate the gas pressure from the input density and entropy
//                --> Here entropy is defined as "pressure / density^(Gamma-1)"
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by the functions "CPU_Shared_FullStepUpdate, ..."
//                3. A floor value "MinPres" is applied to the returned pressure
//
// Parameter   :  Dens     : Mass density
//                Entropy  : Entropy
//                Gamma_m1 : Adiabatic index - 1.0
//                MinPres  : Minimum allowed pressure
//
// Return      :  Pres
//-------------------------------------------------------------------------------------------------------
real CPU_DensEntropy2Pres( const real Dens, const real Entropy, const real Gamma_m1, const real MinPres )
{

   real Pres;

// calculate pressure
   Pres = Entropy*POW( Dens, Gamma_m1 );

// apply a floor value
   Pres = CPU_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : CPU_DensEntropy2Pres
#endif // #if ( DUAL_ENERGY == DE_ENTROPY )



#endif // #if ( MODEL == HYDRO  &&  defined DUAL_ENERGY )
