#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Microphysics_Init
// Description :  Initialize the microphysics routines
//
// Note        :  1. Invoked by Init_GAMER()
//             :  2. Must be called AFTER Init_Load_Parameter() and Init_Unit()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Microphysics_Init()
{

// check if microphysics has been initialized already
   static bool MicroPhy_Initialized = false;

   if ( MicroPhy_Initialized )  return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   MicroPhy.useless            = NULL_BOOL;  // to avoid the empty structure issue
#  ifdef CR_DIFFUSION
   MicroPhy.CR_safety          = DT__CR_DIFFUSION;
   MicroPhy.CR_diff_coeff_para = CR_DIFF_PARA;
   MicroPhy.CR_diff_coeff_perp = CR_DIFF_PERP;
   MicroPhy.CR_diff_min_b      = CR_DIFF_MIN_B;
#  endif // #ifdef CR_DIFFUSION

#  ifdef CONDUCTION
   // Runtime parameters
   MicroPhy.CondType = CONDUCTION_TYPE;
   MicroPhy.CondFluxType = CONDUCTION_FLUX_TYPE;
   MicroPhy.Cond_safety = DT__CONDUCTION;
   MicroPhy.CondConstCoeff = CONDUCTION_CONSTANT_COEFF;
   MicroPhy.CondSpitzerFraction = CONDUCTION_SPITZER_FRACTION;
   MicroPhy.CondCoulombLog = CONDUCTION_COULOMB_LOG;
   MicroPhy.CondMaxDiffusivity = CONDUCTION_MAX_DIFFUSIVITY;

   // Compute the prefactor for the Spitzer conducitivity
   // This prefactor is in CGS (erg/s/cm/K), assuming input temperature in
   // units of 10^7 K, and we must convert it to code units per K. To avoid precision 
   // errors, we do this one step at a time. Still need to correct for differences 
   // in mean molecular weight
   MicroPhy.CondPrefactor = (real)5818590894709.818/MicroPhy.CondCoulombLog;
   MicroPhy.CondPrefactor /= UNIT_E;
   MicroPhy.CondPrefactor *= UNIT_L;
   MicroPhy.CondPrefactor *= UNIT_T;
#  endif // #ifdef CONDUCTION

#  ifdef VISCOSITY
   // Runtime parameters
   MicroPhy.ViscType = VISCOSITY_TYPE;
   MicroPhy.ViscDiffType = VISCOSITY_DIFF_TYPE;
   MicroPhy.ViscCoeffType = VISCOSITY_COEFF_TYPE;
   MicroPhy.Visc_safety = DT__VISCOSITY;
   MicroPhy.ViscKineticCoeff = VISCOSITY_KINETIC_COEFF;
   MicroPhy.ViscDynamicCoeff = VISCOSITY_DYNAMIC_COEFF;
   MicroPhy.ViscSpitzerFraction = VISCOSITY_SPITZER_FRACTION;
   MicroPhy.ViscCoulombLog = VISCOSITY_COULOMB_LOG;
   MicroPhy.ViscMaxDiffusivity = VISCOSITY_MAX_DIFFUSIVITY;

   // Compute the prefactor for the Spitzer viscosity
   // This prefactor is in CGS (g/cm/s), assuming input temperature in units 
   // of 10^7 K, and we must convert it to code units. To avoid precision errors, 
   // we do this one step at a time. Still need to correct for differences in mean 
   // molecular weight
   MicroPhy.ViscPrefactor = (real)695.7010852370435/MicroPhy.ViscCoulombLog;
   MicroPhy.ViscPrefactor /= UNIT_M;
   MicroPhy.ViscPrefactor *= UNIT_L;
   MicroPhy.ViscPrefactor *= UNIT_T;
#  endif // #ifdef VISCOSITY

   MicroPhy_Initialized = true;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Microphysics_Init
