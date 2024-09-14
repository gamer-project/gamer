#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  TrilinearInterpolation
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real TrilinearInterpolation( real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz )
{
   real x1, y1, z1, x0, y0, z0, xd, yd, zd, x, y, z;
   real c000, c001, c010, c100, c011, c101, c110, c111, c00, c01, c10, c11, c0, c1, c;

   x0 = xyz000[0];
   y0 = xyz000[1];
   z0 = xyz000[2];

   x1 = xyz000[0] + dxyz[0];
   y1 = xyz000[1] + dxyz[1];
   z1 = xyz000[2] + dxyz[2];

   x = xyz[0];
   y = xyz[1];
   z = xyz[2];

   c000 = FieldAtVertices[0];
   c001 = FieldAtVertices[1];
   c010 = FieldAtVertices[2];
   c100 = FieldAtVertices[3];
   c011 = FieldAtVertices[4];
   c101 = FieldAtVertices[5];
   c110 = FieldAtVertices[6];
   c111 = FieldAtVertices[7];

   xd = (x-x0)/(x1-x0);
   yd = (y-y0)/(y1-y0);
   zd = (z-z0)/(z1-z0);

   c00 = c000*(1.0-xd) + c100*xd;
   c01 = c001*(1.0-xd) + c101*xd;
   c10 = c010*(1.0-xd) + c110*xd;
   c11 = c011*(1.0-xd) + c111*xd;

   c0  = c00*(1.0-yd) + c10*yd;
   c1  = c01*(1.0-yd) + c11*yd;

   c = c0*(1.0-zd) + c1*zd;

   return c;
} // FUNCTION : TrilinearInterpolation



//-------------------------------------------------------------------------------------------------------
// Function    :  LinearInterpolation
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real LinearInterpolation( real x1, real y1, real x2, real y2, real x )
{
  real dx = x2 - x1;
  return ( y2 * ( x - x1 ) + y1 * ( x2 - x ) ) / dx;
} //FUNCTION : LinearInterpolation
