
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )

static bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real Threshold );
# ifdef GRAVITY
static bool SF_CreateStar_Check_GasJeansLength( const real GasDensity, const real GasCs2, const real Threshold );
# endif



//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check
// Description :  Check if the target cell (i,j,k) satisfies the star formation criteria
//
// Note        :  1. Useless input arrays are set to NULL (e.g, Pot[] if GRAVITY is off)
//
// Parameter   :  lv            : Target refinement level
//                PID           : Target patch ID
//                i,j,k         : Indices of the target cell
//                dh            : Cell size at the target level
//                fluid         : Input fluid array (with NCOMP_TOTAL components)
//                Pres          : Input pressure array
//                Cs2           : Input squared sound speed array
//
// Return      :  "true"  if all  of the star formation criteria is satisfied
//                "false" if none of the star formation criteria is satisfied
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check( const int lv, const int PID, const int i, const int j, const int k, const double dh,
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
         AllowSF &= SF_CreateStar_Check_GasDensity( fluid[DENS][k][j][i], SF_CREATE_STAR_MIN_GAS_DENS );
         if ( !AllowSF )    return AllowSF;
         break;

      case SF_CREATE_STAR_SCHEME_DWARFGALAXY:
#        ifdef GRAVITY
//       create star particles only if the gas Jeans length is less than the given threshold
         AllowSF &= SF_CreateStar_Check_GasJeansLength( fluid[DENS][k][j][i], Cs2[k][j][i], dh*SF_CREATE_STAR_MAX_GAS_JEANSL );
#        endif
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
// Function    :  SF_CreateStar_Check_GasDensity
// Description :  Check if the gas desnity exceeds the given threshold
//
// Note        :  1.
//
// Parameter   :  GasDensity : Gas density
//                Threshold  : Threshold for the star formation
//
// Return      :  "true"  if the gas density is larger           than the given threshold
//                "false" if the gas density is equal or smaller than the given threshold
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasDensity( const real GasDensity, const real Threshold )
{
   bool AllowSF = false;

   if ( GasDensity >= Threshold )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasDensity



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  SF_CreateStar_Check_GasJeansLength
// Description :  Check if the gas Jeans length is below the given threshold
//
// Note        :  1. Gas Jeans length = \sqrt{ \frac{ \pi Cs^2 }{ G \rho } }
//
// Parameter   :  GasDensity : Gas density
//                GasCs2     : Gas squared sound speed
//                Threshold  : Threshold for the star formation
//
// Return      :  "true"  if the gas Jeans length is smaller         than the given threshold
//                "false" if the gas Jeans length is equal or larger than the given threshold
//-------------------------------------------------------------------------------------------------------
bool SF_CreateStar_Check_GasJeansLength( const real GasDensity, const real GasCs2, const real Threshold )
{
   bool AllowSF = false;

   real GasJeansL = SQRT( ( M_PI * GasCs2 ) / ( NEWTON_G * GasDensity ) );

   if ( GasJeansL <= Threshold )    AllowSF = true;

   return AllowSF;

} // FUNCTION : SF_CreateStar_Check_GasDensity
#endif // #ifdef GRAVITY



#endif // #if ( defined PARTICLE  &&  defined STAR_FORMATION  &&  MODEL == HYDRO )
