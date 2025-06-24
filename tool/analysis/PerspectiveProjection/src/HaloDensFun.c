#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  haloDensFun
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real haloDensFun( const real x, const real y, const real z )
{
  const real r = SQRT( x*x + y*y + z*z );

  return PEAK_DENS * POW( (real)1.0 + SQR( r / CORE_RADIUS ), -(real)1.5 * BETA );
} // FUNCTION : haloDensFun



//-------------------------------------------------------------------------------------------------------
// Function    :  xRayHalo
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real xRayHalo( const real x, const real y, const real z, int numRow, real *tempTable, real *lambdaTable )
{
  real dens   = haloDensFun(x, y, z) * MASS_PROTON_GRAM * MU_ELECTRON;

  real lambda = Lambda( HALO_TEMP, numRow, tempTable, lambdaTable );

  return Xray_emissivity( dens, lambda );
} // FUNCTION : xRayHalo
