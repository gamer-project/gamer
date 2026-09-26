#ifndef __CUFLU_DUALENERGY__
#define __CUFLU_DUALENERGY__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined DUAL_ENERGY  &&  !defined SRHD )



// internal functions
#ifdef __CUDACC__
GPU_DEVICE
static real Hydro_DensPres2Dual( const real Dens, const real Pres, const real Passive[],
                                 const EoS_DP2E_t EoS_DensPres2Eint,
                                 const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                                 const real *const EoS_Table[EOS_NTABLE_MAX] );
GPU_DEVICE
static real Hydro_DensDual2Pres( const real Dens, const real Dual, const real Passive[],
                                 const bool CheckMinPres, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                 const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                                 const real *const EoS_Table[EOS_NTABLE_MAX] );
#endif


// ****************************************************************
// Unless otherwise specified, the internal energy and pressure
// in all dual-energy routines EXCLUDE cosmic rays
// --> This is different from the cosmic-ray EoS driver, where the
//     internal energy and pressure always include cosmic rays
// ****************************************************************


//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DualEnergyFix
// Description :  Correct the internal and total energies using the dual-energy formalism
//
// Note        :  1. Invoked by Hydro_FullStepUpdate(), InterpolateGhostZone(), ...
//                2. A floor value "MinPres" is applied to the corrected pressure if CheckMinPres is on
//                3. A floor value "TINY_NUMBER" is applied to the input dual-energy variable as well
//                4. Call-by-reference for "Etot, Dual, and DE_Status"
//                5. The dual-energy variable is determined by DUAL_ENERGY, which can be either
//                   DE_ENPY (entropy) or DE_EINT (internal energy)
//                   --> DE_ENPY: entropy
//                       DE_EINT: internal_energy
//                   --> Note that the entropy here is a monotonic function of entropy per volume
//                       instead of the real thermodynamic entropy (see Eqs. 48 and 49 in the Arepo code paper)
//                6. Fluid variables returned by this function are guaranteed to be consistent with each other
//                   --> It doesn't matter we use the dual-energy variable to correct Eint or vice versa,
//                       and it also holds even when the floor value is applied to pressure
//                7. Optionally including cosmic rays
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Etot              : Total energy density
//                Dual              : Dual-energy variable
//                Passive           : Passive scalars
//                DE_Status         : Assigned to (DE_UPDATED_BY_ETOT / DE_UPDATED_BY_DUAL / DE_UPDATED_BY_MIN_PRES)
//                                    to indicate whether this cell is updated by the total energy, dual-energy variable,
//                                    or pressure floor (MinPres)
//                CheckMinPres      : Return Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which we don't want to enable this option
//                MinPres           : Minimum allowed pressure
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                                    --> It's actually useless here since it's only relevant for SRHD,
//                                        which does not support dual energy
//                DualEnergySwitch  : if ( Eint/(Ekin+Emag) < DualEnergySwitch ) ==> correct Eint and Etot
//                                    else                                       ==> correct Dual
//                Emag              : Magnetic energy density (0.5*B^2) --> for MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                                    --> EoS_AuxArray_Flt[1/2] store Gamma-1 and 1/(Gamma-1) for EOS_GAMMA
//                EoS_Table         : EoS tables
//
// Return      :  Etot, Dual, DE_Status
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Dual, const real Passive[], char &DE_Status,
                          const bool CheckMinPres, const real MinPres, const long PassiveFloor, const real DualEnergySwitch,
                          const real Emag, const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef COSMIC_RAY
   if ( Passive == NULL )
      printf( "ERROR : Passive == NULL for COSMIC_RAY at file <%s>, line <%d>, function <%s> !!\n", ERROR_INFO );
#  endif
#  endif // GAMER_DEBUG


#  if ( DUAL_ENERGY == DE_ENPY )
   const real  Gamma_m1 = EoS_AuxArray_Flt[1];
   const real _Gamma_m1 = EoS_AuxArray_Flt[2];
#  endif
   const bool CheckMinPres_No = false;
   const bool CheckMinEint_No = false;

// apply the dual-energy floor
   Dual = FMAX( Dual, TINY_NUMBER );


// calculate energies
// --> note that here Eint can even be negative due to numerical errors
// --> Enth (i.e., non-thermal energy) includes kinetic, magnetic, and cosmic-ray energies
   real Enth, Eint, Pres;

   Eint = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Etot, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );
// exclude cosmic-ray energy
#  ifdef COSMIC_RAY
   Eint -= Passive[ CRAY-NCOMP_FLUID ];
#  endif
   Enth = Etot - Eint;


// determine whether or not to use the dual-energy variable to correct the total energy density
   if ( Eint/Enth < DualEnergySwitch )
   {
//    correct total energy
//    --> we will apply pressure floor later
      Pres      = Hydro_DensDual2Pres( Dens, Dual, Passive, CheckMinPres_No, NULL_REAL,
                                       EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     if   ( DUAL_ENERGY == DE_ENPY )
      Eint      = Pres*_Gamma_m1;
#     elif ( DUAL_ENERGY == DE_EINT )
      Eint      = Dual;
#     endif
      Etot      = Enth + Eint;
      DE_Status = DE_UPDATED_BY_DUAL;
   }

   else
   {
//    correct dual-energy variable
#     if   ( DUAL_ENERGY == DE_ENPY )
      Pres      = Eint*Gamma_m1;
      Dual      = Hydro_DensPres2Dual( Dens, Pres, Passive, EoS_DensPres2Eint,
                                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     elif ( DUAL_ENERGY == DE_EINT )
      Pres      = Hydro_DensDual2Pres( Dens, Eint, Passive, CheckMinPres_No, NULL_REAL,
                                       EoS_DensEint2Pres, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
      Dual      = Eint;
#     endif
      DE_Status = DE_UPDATED_BY_ETOT;
   } // if ( Eint/Enth < DualEnergySwitch ) ... else ...


// apply pressure floor
   if ( CheckMinPres  &&  Pres < MinPres )
   {
      Pres = MinPres;

//    ensure that both energy and dual-energy variable are consistent with the pressure floor
      Dual = Hydro_DensPres2Dual( Dens, Pres, Passive, EoS_DensPres2Eint,
                                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#     if   ( DUAL_ENERGY == DE_ENPY )
      Eint      = Pres*_Gamma_m1;
#     elif ( DUAL_ENERGY == DE_EINT )
      Eint      = Dual;
#     endif
      Etot      = Enth + Eint;
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
//                3. Currently this function does NOT apply pressure floor when calling Hydro_Con2Pres()
//                   --> However, note that Hydro_DensPres2Dual() does apply a floor value (TINY_NUMBER) to the
//                       dual-energy variable
//
// Parameter   :  Dens              : Mass density
//                MomX/Y/Z          : Momentum density
//                Engy              : Total energy density
//                Passive           : Passive scalars
//                Emag              : Magnetic energy density (0.5*B^2) --> for MHD only
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_CREint2CRPres : EoS routine to compute the cosmic-ray pressure
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                                    --> EoS_AuxArray_Flt[1/2] store Gamma-1 and 1/(Gamma-1) for EOS_GAMMA
//                EoS_Table         : EoS tables
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                                    --> It's actually useless here since it's only relevant for SRHD,
//                                        which does not support dual energy
//
// Return      :  Dual
//-------------------------------------------------------------------------------------------------------
real Hydro_Con2Dual( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_CRE2CRP_t EoS_CREint2CRPres,
                     const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const long PassiveFloor )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef COSMIC_RAY
   if ( Passive == NULL )
      printf( "ERROR : Passive == NULL for COSMIC_RAY at file <%s>, line <%d>, function <%s> !!\n", ERROR_INFO );
#  endif
#  endif // GAMER_DEBUG


   real Dual;

#  if   ( DUAL_ENERGY == DE_ENPY )
// currently this function does NOT apply pressure floor when calling Hydro_Con2Pres()
   const bool CheckMinPres_No = false;
   real Pres;

// calculate pressure
   Pres = Hydro_Con2Pres( Dens, MomX, MomY, MomZ, Engy, Passive, CheckMinPres_No, NULL_REAL, PassiveFloor, Emag,
                          EoS_DensEint2Pres, NULL, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                          EoS_Table, NULL );

// exclude cosmic-ray pressure since Hydro_Con2Pres() returns gas+cosmic-ray pressures
#  ifdef COSMIC_RAY
   const real E_CR = Passive[ CRAY-NCOMP_FLUID ];
   Pres -= EoS_CREint2CRPres( E_CR, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  endif

// convert gas pressure to the dual-energy variable
   Dual = Hydro_DensPres2Dual( Dens, Pres, NULL, NULL,
                               EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

#  elif ( DUAL_ENERGY == DE_EINT )
// currently this function does NOT apply internal energy floor when calling Hydro_Con2Eint()
   const bool CheckMinEint_No = false;

// calculate internal energy
   Dual = Hydro_Con2Eint( Dens, MomX, MomY, MomZ, Engy, CheckMinEint_No, NULL_REAL, PassiveFloor, Emag,
                          NULL, NULL, NULL, NULL, NULL );

// exclude cosmic-ray energy since Hydro_Con2Eint() returns gas+cosmic-ray energies
#  ifdef COSMIC_RAY
   Dual -= Passive[ CRAY-NCOMP_FLUID ];
#  endif

#  endif // DUAL_ENERGY

   return Dual;

} // FUNCTION : Hydro_Con2Dual
#endif // ifndef __CUDACC__



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensPres2Dual
// Description :  Evaluate the dual-energy variable from the input density and pressure
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_Con2Dual() and Hydro_DualEnergyFix()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value (TINY_NUMBER) is applied to the returned value
//                4. Both the input pressure and returned entropy/internal energy exclude cosmic rays
//                   --> Convert to/from the total quantities required by EOS_COSMIC_RAY internally
//
// Parameter   :  Dens              : Mass density
//                Pres              : Pressure
//                Passive           : Passive scalars
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                                    --> EoS_AuxArray_Flt[1/2] store Gamma-1 and 1/(Gamma-1) for EOS_GAMMA
//                EoS_Table         : EoS tables
//
// Return      :  Dual
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
real Hydro_DensPres2Dual( const real Dens, const real Pres, const real Passive[],
                          const EoS_DP2E_t EoS_DensPres2Eint,
                          const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{

   real Dual;

// calculate the dual-energy variable
#  if   ( DUAL_ENERGY == DE_ENPY )
   const real Gamma_m1 = EoS_AuxArray_Flt[1];
   Dual = Pres*POW( Dens, -Gamma_m1 );

#  elif ( DUAL_ENERGY == DE_EINT )
#  ifdef COSMIC_RAY
   const real GammaCR_m1 = EoS_AuxArray_Flt[5];
   const real E_CR       = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR    = GammaCR_m1*E_CR;

   Dual  = EoS_DensPres2Eint( Dens, Pres+Pres_CR, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Dual -= E_CR;
#  else
   Dual  = EoS_DensPres2Eint( Dens, Pres,         Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  endif

#  endif // DUAL_ENERGY

// apply a floor value
   Dual = FMAX( Dual, TINY_NUMBER );

   return Dual;

} // FUNCTION : Hydro_DensPres2Dual



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DensDual2Pres
// Description :  Evaluate the gas pressure from the input density and dual-energy variable
//
// Note        :  1. Used by the dual-energy formalism
//                2. Invoked by Hydro_DualEnergyFix(), Flu_Close(), Hydro_Aux_Check_Negative(), and Flu_FixUp()
//                   --> This function is invoked by both CPU and GPU codes
//                3. A floor value "MinPres" is applied to the returned pressure if CheckMinPres is on
//                4. Both the input entropy/internal energy and returned pressure exclude cosmic rays
//                   --> Convert to/from the total quantities required by EOS_COSMIC_RAY internally
//
// Parameter   :  Dens              : Mass density
//                Dual              : Dual-energy variable
//                Passive           : Passive scalars
//                CheckMinPres      : Return Hydro_CheckMinPres()
//                                    --> In some cases we actually want to check if pressure becomes unphysical,
//                                        for which we don't want to enable this option
//                MinPres           : Minimum allowed pressure
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                                    --> EoS_AuxArray_Flt[1/2] store Gamma-1 and 1/(Gamma-1) for EOS_GAMMA
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
#  ifdef COSMIC_RAY
   const real GammaCR_m1 = (real)EoS_AuxArray_Flt[5];
   const real E_CR       = Passive[ CRAY-NCOMP_FLUID ];
   const real Pres_CR    = GammaCR_m1*E_CR;

   Pres  = EoS_DensEint2Pres( Dens, Dual+E_CR, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   Pres -= Pres_CR;
#  else
   Pres  = EoS_DensEint2Pres( Dens, Dual,      Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
#  endif

#  endif // DUAL_ENERGY

// apply a floor value
   if ( CheckMinPres )  Pres = Hydro_CheckMinPres( Pres, MinPres );

   return Pres;

} // FUNCTION : Hydro_DensDual2Pres



#endif // #if ( MODEL == HYDRO  &&  defined DUAL_ENERGY  &&  !defined SRHD )



#endif // #ifndef __CUFLU_DUALENERGY__
