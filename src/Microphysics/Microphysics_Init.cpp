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
   MicroPhy.CondType = CONDUCTION_TYPE;
   MicroPhy.CondFluxType = CONDUCTION_FLUX_TYPE;
   MicroPhy.Cond_safety = DT__CONDUCTION;
   MicroPhy.CondConstCoeff = CONDUCTION_CONSTANT_COEFF;
   MicroPhy.CondSpitzerFraction = CONDUCTION_SPITZER_FRACTION;
   MicroPhy.CondCoulombLog = CONDUCTION_COULOMB_LOG;
   MicroPhy.CondMaxDiffusivity = CONDUCTION_MAX_DIFFUSIVITY;
   MicroPhy.CondSaturation = CONDUCTION_SATURATION;
   MicroPhy.CondMue = CONDUCTION_MUE;
   if ( OPT_UNIT )
      MicroPhy.CondSpecificHeat = Const_kB / ( MOLECULAR_WEIGHT * MU_NORM ) * (UNIT_M/UNIT_E);
   else
      MicroPhy.CondSpecificHeat = (real)1.0/MOLECULAR_WEIGHT;
   MicroPhy.CondSpecificHeat /= ( GAMMA - (real)1.0 );

   if ( MicroPhy.CondSaturation ) 
   {
      // This calculates the prefactor for the electron MFP
      // Note that this assumes CGS units for the electron charge
      MicroPhy.CondMFPConst = (real)0.7329037678543799 * Const_kB * Const_kB / ( UNIT_E * UNIT_E );
      MicroPhy.CondMFPConst /= POW( Const_e / SQRT( UNIT_E*UNIT_L ), (real)4.0 ) * MicroPhy.CondCoulombLog;
      MicroPhy.CondMFPConst *= CONDUCTION_MUE * Const_amu / UNIT_M;
   }      

   if ( MicroPhy.CondType == CONSTANT_CONDUCTION ) 
   {
      // This coefficient is in CGS (erg/s/cm/K), and we must convert it to code 
      // units per K. To avoid precision errors, we do this one step at a time.
      MicroPhy.CondConstCoeff /= UNIT_E;
      MicroPhy.CondConstCoeff *= UNIT_L;
      MicroPhy.CondConstCoeff *= UNIT_T;
   }
   else if ( MicroPhy.CondType == SPITZER_CONDUCTION ) 
   {
      // Compute the prefactor for the Spitzer conducitivity
      // This prefactor is in CGS (erg/s/cm/K), assuming input temperature in
      // units of 10^7 K, and we must convert it to code units per K. To avoid precision 
      // errors, we do this one step at a time. Still need to correct for differences 
      // in mean molecular weight
      MicroPhy.CondPrefactor = (real)5818590894709.818/MicroPhy.CondCoulombLog;
      MicroPhy.CondPrefactor /= UNIT_E;
      MicroPhy.CondPrefactor *= UNIT_L;
      MicroPhy.CondPrefactor *= UNIT_T;
      MicroPhy.CondPrefactor *= MicroPhy.CondSpitzerFraction;
   } 
   else 
   {
      Aux_Error( ERROR_INFO, "unsupported conduction type (%d) !!\n", MicroPhy.CondType );
   }

#  endif // #ifdef CONDUCTION

#  ifdef VISCOSITY
   MicroPhy.ViscType = VISCOSITY_TYPE;
   MicroPhy.ViscFluxType = VISCOSITY_FLUX_TYPE;
   MicroPhy.ViscCoeffType = VISCOSITY_COEFF_TYPE;
   MicroPhy.Visc_safety = DT__VISCOSITY;
   MicroPhy.ViscConstCoeff = VISCOSITY_CONST_COEFF;
   MicroPhy.ViscSpitzerFraction = VISCOSITY_SPITZER_FRACTION;
   MicroPhy.ViscCoulombLog = VISCOSITY_COULOMB_LOG;
   MicroPhy.ViscMaxDiffusivity = VISCOSITY_MAX_DIFFUSIVITY;

   if ( MicroPhy.CondType == CONSTANT_VISCOSITY ) 
   {
      // We must convert the viscosity coefficient to code units. To avoid 
      // precision errors, we do this one step at a time. 
      if ( MicroPhy.ViscCoeffType == VISCOSITY_KINETIC_COEFF ) {
         // This coefficient is in units of cm^2/s
         MicroPhy.ViscConstCoeff /= UNIT_L;
         MicroPhy.ViscConstCoeff /= UNIT_L;
         MicroPhy.ViscConstCoeff *= UNIT_T;
      } else if ( MicroPhy.ViscCoeffType == VISCOSITY_DYNAMIC_COEFF ) {
         // This coefficient is in units of g/cm/s
         MicroPhy.ViscConstCoeff /= UNIT_M;
         MicroPhy.ViscConstCoeff *= UNIT_L;
         MicroPhy.ViscConstCoeff *= UNIT_T;
      }
   } 
   else if ( MicroPhy.ViscType == SPITZER_VISCOSITY ) 
   {
      // Compute the prefactor for the Spitzer viscosity
      // This prefactor is in CGS (g/cm/s), assuming input temperature in units 
      // of 10^7 K, and we must convert it to code units. To avoid precision errors, 
      // we do this one step at a time. Still need to correct for differences in mean 
      // molecular weight
      MicroPhy.ViscPrefactor = (real)695.7010852370435/MicroPhy.ViscCoulombLog;
      MicroPhy.ViscPrefactor /= UNIT_M;
      MicroPhy.ViscPrefactor *= UNIT_L;
      MicroPhy.ViscPrefactor *= UNIT_T;
      MicroPhy.ViscPrefactor *= MicroPhy.ViscSpitzerFraction;
   }
   else 
   {
      Aux_Error( ERROR_INFO, "unsupported viscosity type (%d) !!\n", MicroPhy.ViscType );
   }

#  endif // #ifdef VISCOSITY

   MicroPhy_Initialized = true;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Microphysics_Init
