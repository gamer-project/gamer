#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )

static bool SF_CreateStar_Check_CellMassDepletion( const real GasMass );
static bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real CosmoScaleFactor, const real Threshold );
static bool SF_CreateStar_Check_GasOverdensity( const real GasDensity );
static bool SF_CreateStar_Check_GasTemperature( const real GasTemperature, const real CosmoScaleFactor, const real Threshold );
static bool SF_CreateStar_Check_GasJeansLength( const real GasDensity, const real GasCs2, const real CosmoScaleFactor, const real Threshold );




//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check
// Description :  Check if the target cell (i,j,k) satisfies the star formation criteria
//
// Note        :  1. Useless input arrays are set to NULL
//                2. Each scheme can have a combination of multiple star formation criteria
//                   --> A star can only form when all of the criteria are met
//                   --> Different combinations of choices are represented by different schemes
//
// Parameter   :  lv               : Target refinement level
//                PID              : Target patch ID
//                i,j,k            : Indices of the target cell
//                dh               : Cell size at the target level
//                CosmoScaleFactor : Scale factor "a" in cosmology
//                                   --> Must be set to unity when COMOVING is disabled
//                fluid            : Input fluid array (with NCOMP_TOTAL components)
//                Temp             : Input temperature array
//                Pres             : Input pressure array
//                Cs2              : Input squared sound speed array
//
// Return      :  "true"  if the specified star formation criteria are satisfied
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check( const int lv, const int PID, const int i, const int j, const int k, const double dh, const real CosmoScaleFactor,
                          const real fluid[][PS1][PS1][PS1], const real Temp[][PS1][PS1], const real Pres[][PS1][PS1], const real Cs2[][PS1][PS1] )
{
#  ifdef GAMER_DEBUG
   const long SupportedCriteria = ( SF_CREATE_STAR_CRITERIA_HIGH_GAS_DENSITY | SF_CREATE_STAR_CRITERIA_LOW_GAS_TEMPERATURE | SF_CREATE_STAR_CRITERIA_UNRESOLVED_JEANS_LENGTH );
   if ( SF_CREATE_STAR_CRITERIA & ~SupportedCriteria )
      Aux_Error( ERROR_INFO, "unsupported SF_CREATE_STAR_CRITERIA = %ld !!\n", SF_CREATE_STAR_CRITERIA );
#  endif

   bool AllowSF = true;

   if ( SF_CREATE_STAR_CRITERIA == SF_CREATE_STAR_CRITERIA_NONE )
   {
      AllowSF &= false;
      if ( !AllowSF )    return AllowSF;
   }

   if ( SF_CREATE_STAR_CRITERIA & SF_CREATE_STAR_CRITERIA_HIGH_GAS_DENSITY )
   {
//    create star particles only if the gas density is higher than the given threshold
      AllowSF &= SF_CreateStar_Check_GasDensity( fluid[DENS][k][j][i], CosmoScaleFactor, SF_CREATE_STAR_MIN_GAS_DENS );
      if ( !AllowSF )    return AllowSF;
   }

   if ( SF_CREATE_STAR_CRITERIA & SF_CREATE_STAR_CRITERIA_LOW_GAS_TEMPERATURE )
   {
//    create star particles only if the gas temperature is lower than the given threshold
      AllowSF &= SF_CreateStar_Check_GasTemperature( Temp[k][j][i], CosmoScaleFactor, SF_CREATE_STAR_MAX_GAS_TEMP );
      if ( !AllowSF )    return AllowSF;
   }

   if ( SF_CREATE_STAR_CRITERIA & SF_CREATE_STAR_CRITERIA_UNRESOLVED_JEANS_LENGTH )
   {
//    create star particles only if the gas Jeans length is less than the given threshold
      AllowSF &= SF_CreateStar_Check_GasJeansLength( fluid[DENS][k][j][i], Cs2[k][j][i], CosmoScaleFactor, dh*SF_CREATE_STAR_MAX_GAS_JEANSL );
      if ( !AllowSF )    return AllowSF;
   }


   return AllowSF;

} // FUNCTION : SF_CreateStar_Check



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_CellMassDepletion
// Description :  Check if the gas mass is sufficient to spawn a stellar particle
//
// Note        :  1. The gas mass should be greater than the stellar particle mass
//                   by a factor of 1/(allowed depletion fraction),
//                   even in the stochastic star formation model,
//                   to avoid spawning a stellar particle of a capped mass
//                2. It can be seen as an effective physical density threshold of
//                   MinStarMass / ( MaxStarMFrac * dh^3 * a^3 ),
//                   where dh is the comoving cell size
//                3. This function is not used currently
//                   --> It will be used in cosmological simulations in the future
//
// Parameter   :  GasMass : Gas mass in the cell
//
// Return      :  "true"  if the gas mass is sufficient to spawn a stellar particle
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_CellMassDepletion( const real GasMass )
{

   bool AllowSF = false;

   const real MinStarFormationCellMass = SF_CREATE_STAR_MIN_STAR_MASS / SF_CREATE_STAR_MAX_STAR_MFRAC;

   if ( GasMass >= MinStarFormationCellMass )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_CellMassDepletion



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasDensity
// Description :  Check if the gas density exceeds the given threshold
//
// Note        :  1. The density threshold is defined in the physical frame even when COMOVING is enabled
//                2. When COMOVING is enabled, "GasDensity" should be the comoving density, a^3*\rho,
//                   where \rho is the density in the physical frame
//
// Parameter   :  GasDensity       : Gas density
//                CosmoScaleFactor : Scale factor "a" in cosmology
//                                   --> Must be set to unity when COMOVING is disabled
//                Threshold        : Threshold for the star formation
//
// Return      :  "true"  if the gas density is larger than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real CosmoScaleFactor, const real Threshold )
{

   const real a3inv = (real)1.0 / CUBE( CosmoScaleFactor );  // a^-3

   bool AllowSF = false;

// converted into the density in the physical frame
   if ( GasDensity * a3inv >= Threshold )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasDensity



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasOverdensity
// Description :  Check if the gas overdensity exceeds the given threshold
//
// Note        :  1. COMOVING must be enabled, so the input "GasDensity" is
//                   the comoving density in units of the background matter density,
//                   whose value equals a^3*\rho_{gas} / (\Omega_{m,0}*\rho_{crit,0}),
//                   where \rho_{gas} is the local gas density in the physical frame,
//                         \Omega_{m,0} is the current density parameter of matter, and
//                         \rho_{crit,0} is the current critical density
//                2. The overdensity is defined as a dimensionless quantity:
//                   \delta_{gas} = \rho_{gas} / (\langle \rho_{gas} \rangle)
//                                = \rho_{gas} / (\Omega_b * \rho_{crit})
//                                = \rho_{gas} / (\Omega_b/\Omega_m * \Omega_m*\rho_{crit})
//                                = \rho_{gas} / (\Omega_b/\Omega_m * a^{-3}*\Omega_{m,0}*\rho_{crit,0})
//                                = "GasDensity" / (\Omega_b/\Omega_m),
//                   where \langle \rho_{gas} \rangle is the background gas density,
//                         \rho_{crit} is the critical density,
//                         \Omega_b is the density parameter of baryons, and
//                         \Omega_m is the density parameter of matter
//                3. This function is not used currently
//                   --> It will be used in cosmological simulations in the future
//
// Parameter   :  GasDensity : Gas density
//
// Return      :  "true"  if the gas overdensity is larger than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasOverdensity( const real GasDensity )
{

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

   bool AllowSF = false;

#  ifdef COMOVING
   const double Omega_B0          = 0.04;                  // Omega_baryon at the present time
   const real   BaryonMatterRatio = Omega_B0 / OMEGA_M0;   // Omega_baryon / Omega_matter
   const real   GasOverdensity    = GasDensity / BaryonMatterRatio;

// The threshold value is hard-coded for now
// --> It should be a runtime parameter in the future
// --> The value is taken from GADGET-4's default; see the references
//     1. The parameter "CritOverDensity" in the documentation
//        https://wwwmpa.mpa-garching.mpg.de/gadget4/#documentation
//     2. The variables "All.OverDensThresh" and "All.CritOverDensity"
//        in GADGET-4's src/cooling_sfr/sfr_eos.cc
   const real CritOverDensity = 57.7;    // overdensity at R200 of an NFW halo
   const real Threshold       = CritOverDensity;

   if ( GasOverdensity >= Threshold )    AllowSF = true;
#  endif // #ifdef COMOVING

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasOverdensity



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasTemperature
// Description :  Check if the gas temperature falls below the given threshold
//
// Note        :  1. The temperature threshold is defined in the physical frame even when COMOVING is enabled.
//                2. When COMOVING is enabled, "GasTemperature" should be the comoving temperature, a^2*T,
//                   where T is the density in the physical frame
//
// Parameter   :  GasTemperature   : Gas temperature
//                CosmoScaleFactor : Scale factor "a" in cosmology
//                                   --> Must be set to unity when COMOVING is disabled
//                Threshold        : Threshold for the star formation
//
// Return      :  "true"  if the gas temperature is lower than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasTemperature( const real GasTemperature, const real CosmoScaleFactor, const real Threshold )
{

   const real a2inv = (real)1.0 / SQR( CosmoScaleFactor );  // a^-2

   bool AllowSF = false;

// converted into the temperature in the physical frame
   if ( GasTemperature * a2inv <= Threshold )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasTemperature



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasJeansLength
// Description :  Check if the gas Jeans length is below the given threshold
//
// Note        :  1. Gas Jeans length = \sqrt{ \frac{ \pi Cs^2 }{ G \rho } }
//                2. When COMOVING is enabled,
//                   (1) Newton G in the above formula should be replaced by G*a
//                   (2) The Jeans length threshold is defined in the comoving frame, a^{-1}*\lambda_J,
//                       where \lambda_J is the Jeans length in the physical frame
//                   (3) "GasDensity" should be the comoving density, a^3*\rho,
//                       where \rho is the density in the physical frame
//                   (4) "GasCs2" should be the comoving sound speed squared, a^2*Cs^2,
//                       where Cs is the sound speed in the physical frame
//
// Parameter   :  GasDensity       : Gas density
//                GasCs2           : Gas squared sound speed
//                CosmoScaleFactor : Scale factor "a" in cosmology
//                                   --> Must be set to unity when COMOVING is disabled
//                Threshold        : Threshold for the star formation
//
// Return      :  "true"  if the gas Jeans length is smaller than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasJeansLength( const real GasDensity, const real GasCs2, const real CosmoScaleFactor, const real Threshold )
{

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

   bool AllowSF = false;

#  ifdef GRAVITY
   const real GasJeansL2 = ( M_PI * GasCs2 ) / ( NEWTON_G * CosmoScaleFactor * GasDensity );

// the threshold is defined in the comoving frame when COMOVING is enabled
   if ( GasJeansL2 <= SQR(Threshold) )    AllowSF = true;
#  endif // #ifdef GRAVITY

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasJeansLength



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
