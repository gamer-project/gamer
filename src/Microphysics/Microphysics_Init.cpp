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

   // This is the numerical constant prefactor for the mean free path,
   // whether electron or ion: (3**1.5)*(k_B**2)*amu/(4*(pi**0.5)*e**4)
   // (see Sarazin 1988, etc.). Note that this assumes CGS units for
   // the elementary charge e
   const double mfp_prefactor = 4.35851e-19;

#  ifdef CONDUCTION
   MicroPhy.CondType = CONDUCTION_TYPE;
   MicroPhy.CondFluxType = CONDUCTION_FLUX_TYPE;
   MicroPhy.Cond_safety = DT__CONDUCTION;
   MicroPhy.CondConstCoeff = CONDUCTION_CONSTANT_COEFF;
   MicroPhy.CondSpitzerFraction = CONDUCTION_SPITZER_FRAC;
   MicroPhy.CondCoulombLog = CONDUCTION_COULOMB_LOG;
   MicroPhy.CondMaxDiffusivity = CONDUCTION_MAX_DIFFUSIVITY;
   MicroPhy.CondSaturation = CONDUCTION_SATURATION;
   MicroPhy.CondSatWhistler = CONDUCTION_SAT_WHISTLER;
   MicroPhy.CondMue = CONDUCTION_MUE;
   if ( OPT__UNIT )
      MicroPhy.CondSpecificHeat = (real)( Const_kB / ( MOLECULAR_WEIGHT * MU_NORM ) * ( UNIT_M/UNIT_E ) );
   else
      MicroPhy.CondSpecificHeat = (real)1.0/MOLECULAR_WEIGHT;
   MicroPhy.CondPresConv = MicroPhy.CondSpecificHeat;
   MicroPhy.CondSpecificHeat /= ( GAMMA - (real)1.0 );

   if ( MicroPhy.CondSaturation )
   {
      // This calculates the prefactor for the electron MFP
      MicroPhy.CondMFPConst = (real)( mfp_prefactor * CONDUCTION_MUE / MicroPhy.CondCoulombLog );
      MicroPhy.CondMFPConst *= (real)( UNIT_L * UNIT_L / UNIT_M );
   }

   if ( MicroPhy.CondType == CONSTANT_CONDUCTION )
   {
      // This coefficient is in CGS (erg/s/cm/K), and we must convert it to code
      // units per K. To avoid precision errors, we do this one step at a time.
      MicroPhy.CondConstCoeff *= (real)( UNIT_T*UNIT_L/UNIT_E );

   }
   else if ( MicroPhy.CondType == SPITZER_CONDUCTION )
   {
      // Compute the prefactor for the Spitzer conducitivity
      // This prefactor is in CGS (erg/s/cm/K), assuming input temperature in
      // units of 10^7 K, and we must convert it to code units per K. To avoid precision
      // errors, we do this one step at a time. Still need to correct for differences
      // in mean molecular weight
      MicroPhy.CondPrefactor = (real)5818590894709.818/MicroPhy.CondCoulombLog;
      MicroPhy.CondPrefactor *= (real)( UNIT_T*UNIT_L/UNIT_E );
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
   MicroPhy.ViscConstCoeff = VISCOSITY_CONSTANT_COEFF;
   MicroPhy.ViscSpitzerFraction = VISCOSITY_SPITZER_FRAC;
   MicroPhy.ViscCoulombLog = VISCOSITY_COULOMB_LOG;
   MicroPhy.ViscMaxDiffusivity = VISCOSITY_MAX_DIFFUSIVITY;
   MicroPhy.ViscBounds = VISCOSITY_BOUNDS;
   MicroPhy.ViscSaturation = VISCOSITY_SATURATION;
   MicroPhy.ViscMui = VISCOSITY_MUI;

   if ( MicroPhy.ViscSaturation )
   {
      // This calculates the prefactor for the ion MFP
      MicroPhy.ViscMFPConst = (real)( mfp_prefactor * VISCOSITY_MUI / MicroPhy.ViscCoulombLog );
      MicroPhy.ViscMFPConst *= (real)( UNIT_L * UNIT_L / UNIT_M );
      // This prefactor converts temperature to thermal speed squared: v_th^2 = 2kT/m_i
      MicroPhy.ViscThermalSpeedConv = (real)( 2.0 * Const_kB / ( VISCOSITY_MUI * Const_amu ) * ( UNIT_M/UNIT_E ) );
   }

   if ( MicroPhy.ViscType == CONSTANT_VISCOSITY )
   {
      // We must convert the viscosity coefficient to code units. To avoid
      // precision errors, we do this one step at a time.
      if ( MicroPhy.ViscCoeffType == VISCOSITY_KINETIC_COEFF ) {
         // This coefficient is in units of cm^2/s
         MicroPhy.ViscConstCoeff *= (real)( UNIT_T/UNIT_L/UNIT_L );
      } else if ( MicroPhy.ViscCoeffType == VISCOSITY_DYNAMIC_COEFF ) {
         // This coefficient is in units of g/cm/s
         MicroPhy.ViscConstCoeff *= (real)( UNIT_T*UNIT_L/UNIT_M );
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
      MicroPhy.ViscPrefactor *= (real)( UNIT_T*UNIT_L/UNIT_M );
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
