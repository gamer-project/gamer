#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )

static bool SF_CreateStar_Check_CellMassDepletion( const real GasMass );
static bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real CosmoScaleFactor, const real Threshold );
static bool SF_CreateStar_Check_GasOverDensity( const real GasDensity );
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
//                fluid            : Input fluid array (with NCOMP_TOTAL components)
//                Pres             : Input pressure array
//                Cs2              : Input squared sound speed array
//
// Return      :  "true"  if the specified star formation criteria are satisfied
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check( const int lv, const int PID, const int i, const int j, const int k, const double dh, const real CosmoScaleFactor,
                          const real fluid[][PS1][PS1][PS1], const real Pres[][PS1][PS1], const real Cs2[][PS1][PS1] )
{
   bool AllowSF = true;

   switch ( SF_CREATE_STAR_SCHEME )
   {
      case SF_CREATE_STAR_SCHEME_AGORA:
//       create star particles only if the gas density exceeds the given threshold
//       Ref: (1) Nathan Goldbaum, et al., 2015, ApJ, 814, 131 (arXiv: 1510.08458), sec. 2.4
//            (2) Ji-hoon Kim, et al., 2016, ApJ, 833, 202 (arXiv: 1610.03066), sec. 3.2
//       Currently this function does not check whether the cell mass exceeds the Jeans mass
//       --> Ref: "jeanmass" in star_maker_ssn.F of Enzo
         AllowSF &= SF_CreateStar_Check_GasDensity( fluid[DENS][k][j][i], CosmoScaleFactor, SF_CREATE_STAR_MIN_GAS_DENS );
         if ( !AllowSF )    return AllowSF;
         break;

      case SF_CREATE_STAR_SCHEME_DWARFGALAXY:
//       create star particles only if the gas Jeans length is less than the given threshold
         AllowSF &= SF_CreateStar_Check_GasJeansLength( fluid[DENS][k][j][i], Cs2[k][j][i], CosmoScaleFactor, dh*SF_CREATE_STAR_MAX_GAS_JEANSL );
         if ( !AllowSF )    return AllowSF;
         break;

      case SF_CREATE_STAR_SCHEME_NONE:
         AllowSF &= false;
         if ( !AllowSF )    return AllowSF;
         break;

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "SF_CREATE_STAR_SCHEME", SF_CREATE_STAR_SCHEME );
   } // switch ( SF_CREATE_STAR_SCHEME )


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
//                   MinStarMass / MaxStarMFrac / dh^3 * a^-3,
//                   where dh is the comoving cell size
//                3. This function is not used currently
//                   --> It will be used in the cosmological simulations in the future
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
// Note        :  1. The density threshold is in physical frame
//
// Parameter   :  GasDensity       : Gas density
//                CosmoScaleFactor : Scale factor "a" in cosmology
//                Threshold        : Threshold for the star formation
//
// Return      :  "true"  if the gas density is larger than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real CosmoScaleFactor, const real Threshold )
{
   const real a3inv = (real)1.0 / CUBE( CosmoScaleFactor );  // a^-3

   bool AllowSF = false;

   if ( GasDensity * a3inv >= Threshold )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasDensity



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasOverDensity
// Description :  Check if the gas overdensity exceeds the given threshold
//
// Note        :  1. The threshold value is hard-coded for now
//                2. This function is not used currently
//                   --> It will be used in the cosmological simulations in the future
//
// Parameter   :  GasDensity : Gas density
//
// Return      :  "true"  if the gas overdensity is larger than or equal to the given threshold
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasOverDensity( const real GasDensity )
{

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

   bool AllowSF = false;

#  ifdef COMOVING
   const real CritOverDensity = 57.7;    // overdensity at R200 of an NFW halo
   const real BaryonRatio     = 0.167;   // Omega_Baryon / Omega_Matter
   const real OverDensThres   = CritOverDensity * BaryonRatio;

   if ( GasDensity >= OverDensThres )    AllowSF = true;
#  endif // #ifdef COMOVING

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasOverDensity



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasJeansLength
// Description :  Check if the gas Jeans length is below the given threshold
//
// Note        :  1. Gas Jeans length = \sqrt{ \frac{ \pi Cs^2 }{ G \rho } }
//                2. In comoving coordinates, Newton G in the formula should be replaced by G*a
//
// Parameter   :  GasDensity       : Gas density
//                GasCs2           : Gas squared sound speed
//                CosmoScaleFactor : Scale factor "a" in cosmology
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

   if ( GasJeansL2 <= SQR(Threshold) )    AllowSF = true;
#  endif // #ifdef GRAVITY

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasJeansLength



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
