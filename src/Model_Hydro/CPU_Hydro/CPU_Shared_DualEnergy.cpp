#ifndef __CUFLU_DUALENERGY__
#define __CUFLU_DUALENERGY__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined DUAL_ENERGY  &&  !defined SRHD )



// internal functions
#ifdef __CUDACC__
GPU_DEVICE
static real Hydro_DensEint2Dual( const real Dens, const real Eint, const real Passive[],
                                 const EoS_DE2S_t EoS_DensEint2Entr, const double EoS_AuxArray_Flt[],
                                 const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] );
GPU_DEVICE
static real Hydro_DensDual2Pres( const real Dens, const real Dual, const real Passive[],
                                 const bool CheckMinPres, const real MinPres,
                                 const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[],
                                 const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DualEnergyFix
// Description :  Correct the internal and total energies using the dual-energy formalism
//
// Note        :  1. Invoked by Hydro_FullStepUpdate(), InterpolateGhostZon(), ...
//                2. A floor value "MinPres" is applied to the corrected pressure if CheckMinPres is on
//                3  A floor value "TINY_NUMBER" is applied to the input dual-energy variable as well
//                4. Call-by-reference for "Etot, Dual, and DE_Status"
//                5. The dual-energy variable is determined by DUAL_ENERGY, which can be either
//                   DE_ENPY (entropy) or DE_EINT (internal energy)
//                   --> DE_ENPY: entropy = pressure / density^(Gamma-1)
//                       DE_EINT: gas internal energy density
//                   --> Note that the entropy here is a monotonic function of entropy per volume
//                       instead of the real thermodynamic entropy (see Eqs. 48 and 49 in the Arepo code paper)
//                6. Fluid variables returned by this function are guaranteed to be consistent with each other
//                   --> It doesn't matter we use the dual-energy variable to correct Eint or vice versa,
//                       and it also holds even when the floor value is applied to pressure
//                7. DE_ENPY only supports the Gamma-law EoS;
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Etot              : Total energy density
//                Dual              : Dual-energy variable
//                DE_Status         : Assigned to (DE_UPDATED_BY_ETOT / DE_UPDATED_BY_DUAL / DE_UPDATED_BY_MIN_PRES)
//                                    to indicate whether this cell is updated by the total energy, dual-energy variable,
//                                    or pressure floor (MinPres)
//                Passive           : Passive scalars
//                CheckMinPres      : Return Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which we don't want to enable this option
//                MinPres           : Minimum allowed pressure
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                DualEnergySwitch  : if ( Eint/(Ekin+Emag) < DualEnergySwitch ) ==> correct Eint and Etot
//                                    else                                       ==> correct Dual
//                Emag              : Magnetic energy density (0.5*B^2) --> for MHD only
//                EoS_DensEint2Pres : EoS routine to compute pressure
//                EoS_DensPres2Eint : EoS routine to compute internal energy density
//                EoS_DensEint2Entr : EoS routine to compute entropy
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//
// Return      :  Etot, Dual, DE_Status
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Dual, char &DE_Status, const real Passive[],
                          const bool CheckMinPres, const real MinPres, const long PassiveFloor, const real DualEnergySwitch,
                          const real Emag, const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                          const EoS_DE2S_t EoS_DensEint2Entr,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{

#  if ( DUAL_ENERGY == DE_ENPY )
   const bool CheckMinPres_No = false;
#  endif
   const bool CheckMinEint_No = false;

// apply the dual-energy floor
   Dual = FMAX( Dual, TINY_NUMBER );


// calculate energies
// --> note that here Eint can even be negative due to numerical errors
// --> Enth (i.e., non-thermal energy) includes both kinetic and magnetic energies
   real Enth, Eint, Pres;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Etot, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
   Enth = Etot - Eint;


// determine whether or not to use the dual-energy variable to correct the total energy density
   if ( Eint/Enth < DualEnergySwitch )
   {
//    correct total energy
//    --> we will apply pressure floor later
#     if   ( DUAL_ENERGY == DE_ENPY )
      Pres      = Hydro_DensDual2Pres( Dens, Dual, Passive, CheckMinPres_No, NULL_REAL,
                                       EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
      Eint      = EoS_DensPres2Eint( Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     elif ( DUAL_ENERGY == DE_EINT )
      Eint      = Dual;
      Pres      = EoS_DensEint2Pres( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     endif
      Etot      = Enth + Eint;
      DE_Status = DE_UPDATED_BY_DUAL;
   }

   else
   {
//    correct dual-energy variable
      Pres      = EoS_DensEint2Pres( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     if   ( DUAL_ENERGY == DE_ENPY )
      Dual      = Hydro_DensEint2Dual( Dens, Eint, Passive, EoS_DensEint2Entr,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     elif ( DUAL_ENERGY == DE_EINT )
      Dual      = Eint;
#     endif
      DE_Status = DE_UPDATED_BY_ETOT;
   } // if ( Eint/Enth < DualEnergySwitch ) ... else ...


// apply pressure floor
   if ( CheckMinPres  &&  Pres < MinPres )
   {
      Pres = MinPres;
      Eint = EoS_DensPres2Eint( Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

//    ensure that both energy and dual-energy variable are consistent with the pressure floor
      Etot      = Enth + Eint;

#     if   ( DUAL_ENERGY == DE_ENPY )
      Dual      = Hydro_DensEint2Dual( Dens, Eint, Passive, EoS_DensEint2Entr,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     elif ( DUAL_ENERGY == DE_EINT )
      Dual      = Eint;
#     endif
      DE_Status = DE_UPDATED_BY_MIN_PRES;
   }

} // FUNCTION : Hydro_DualEnergyFix



// Hydro_Con2Dual() is used by CPU only
#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Con2Dual
// Description :  Evaluate the dual-energy variable from the input fluid variables
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_Init_ByFunction_AssignData(), Gra_Close(), Init_ByFile(), ...
//                3. Currently this function does NOT apply a pressure or dual-energy floor
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Total energy density
//                Emag              : Magnetic energy density (0.5*B^2) --> for MHD only
//                Passive           : Passive scalars
//                EoS_DensEint2Entr : EoS routine to compute the gas entropy
//                EoS_AuxArray_*    : Auxiliary arrays for EoS_DensEint2Entr()
//                EoS_Table         : EoS tables
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//
// Return      :  Dual
//-------------------------------------------------------------------------------------------------------
real Hydro_Con2Dual( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Emag, const real Passive[], const EoS_DE2S_t EoS_DensEint2Entr,  const double EoS_AuxArray_Flt[],
                     const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX],
                     const long PassiveFloor )
{

// currently this function does NOT apply pressure floor when calling Hydro_Con2Pres()
   const bool CheckMin_No = false;
   real Dual;

#  if   ( DUAL_ENERGY == DE_ENPY )
   Dual = Hydro_Con2Entr( Dens, MomX, MomY, MomZ, Engy, Passive, CheckMin_No, NULL_REAL, PassiveFloor,
                          Emag, EoS_DensEint2Entr, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  elif ( DUAL_ENERGY == DE_EINT )
   Dual = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMin_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
#  endif

   return Dual;

} // FUNCTION : Hydro_Con2Dual
#endif // ifndef __CUDACC__



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensEint2Dual
// Description :  Evaluate the dual-energy variable from the input density and internal energy
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_DualEnergyFix(), Grackle_Close(), and Flu_FixUp_Flux()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value (TINY_NUMBER) is applied to the returned value
//
// Parameter   :  Dens              : Mass density
//                Eint              : Internal energy density
//                Passive           : Passive scalars
//                EoS_DensEint2Entr : EoS routine to compute entropy (DE_ENPY only)
//                EoS_AuxArray_*    : EoS auxiliary arrays
//                EoS_Table         : EoS tables
//
// Return      :  Dual
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_DensEint2Dual( const real Dens, const real Eint, const real Passive[],
                          const EoS_DE2S_t EoS_DensEint2Entr, const double EoS_AuxArray_Flt[],
                          const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real Dual;

// calculate the dual-energy variable
#  if   ( DUAL_ENERGY == DE_ENPY )
   Dual = EoS_DensEint2Entr( Dens, Eint, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  elif ( DUAL_ENERGY == DE_EINT )
   Dual = Eint;
#  endif

// apply a floor value
   Dual = FMAX( Dual, TINY_NUMBER );

   return Dual;

} // FUNCTION : Hydro_DensEint2Dual



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensDual2Pres
// Description :  Evaluate the gas pressure from the input density and dual-energy variable
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_DualEnergyFix(), Flu_Close(), Hydro_Aux_Check_Negative(), and Flu_FixUp()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value "MinPres" is applied to the returned pressure if CheckMinPres is on
//
// Parameter   :  Dens         : Mass density
//                Dual         : Dual-energy variable
//                Passive      : Passive scalars
//                CheckMinPres : Return Hydro_CheckMinPres()
//                               --> In some cases we actually want to check if pressure becomes unphysical,
//                                   for which we don't want to enable this option
//                MinPres      : Minimum allowed pressure
//                EoS_DensEint2Pres : EoS routine to compute pressure (DE_EINT only)
//                EoS_AuxArray_*    : EoS auxiliary arrays (DE_ENPY uses Gamma-1 stored at index 1)
//                EoS_Table         : EoS tables
//
// Return      :  Pres
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_DensDual2Pres( const real Dens, const real Dual, const real Passive[],
                          const bool CheckMinPres, const real MinPres,
                          const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[],
                          const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real Pres;

// calculate pressure
#  if   ( DUAL_ENERGY == DE_ENPY )
   const real Gamma_m1 = EoS_AuxArray_Flt[1];
   Pres = Dual*POW( Dens, Gamma_m1 );
#  elif ( DUAL_ENERGY == DE_EINT )
   Pres = EoS_DensEint2Pres( Dens, Dual, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  endif

// apply a floor value
   if ( CheckMinPres )  Pres = Hydro_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : Hydro_DensDual2Pres



#endif // ( MODEL == HYDRO  &&  defined DUAL_ENERGY  &&  !defined SRHD )



#endif // #ifndef __CUFLU_DUALENERGY__
