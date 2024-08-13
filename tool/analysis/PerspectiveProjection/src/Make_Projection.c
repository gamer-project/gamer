#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  PerspectiveProject
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
real PerspectiveProject( real b, real l, real dt, real *XYZ[], int numCellXYZ[], real dxyz[], real azimuthalAngle,
                         real BoxSize[], real ***TargetQuantity, int numRow, real *tempTable, real *lambdaTable )
{
   if ( FABS(COS(l)) < EPSILON )
   {
      printf("l=%e, COS(l)=%e\n", l, COS(l));
      fflush(stdout);
      exit(0);
   }

   real x, y, z, r2, xp, yp, zp;
   real t                  = 0.0;
   double ProjectedValue   = 0.0;

   b *= M_PI/180.0;
   l *= M_PI/180.0;

   do
   {
      if ( -0.5*M_PI < l  &&  l < +0.5*M_PI ) t += dt;
      else                                    t -= dt;

      x = (real)R_SUN - t;
      y = TAN(l) * t;
      z = TAN(b) / COS(l) *  t;

      r2 = x*x + y*y + z*z;

      if( r2 > (real)SQR(HALO_RADIUS) )   break;

      rotationMatrix( x, y, z, &xp, &yp, &zp, azimuthalAngle );

      real xyz[3] = { xp, yp, zp };

#     if ( defined XRAY_ROSAT  ||  defined XRAY_EROSITA )
//    project the X-ray from the Galactic halo
      if ( numRow > 0 )   ProjectedValue += xRayHalo( x, y, z, numRow, tempTable, lambdaTable );
#     endif

      const int Idx = BinarySearch( XYZ[0], 0, numCellXYZ[0]-1, xp );
      const int Jdx = BinarySearch( XYZ[1], 0, numCellXYZ[1]-1, yp );
      const int Kdx = BinarySearch( XYZ[2], 0, numCellXYZ[2]-1, zp );

      if ( (Idx < 0 || Idx > numCellXYZ[0]-2)  ||
           (Jdx < 0 || Jdx > numCellXYZ[1]-2)  ||
           (Kdx < 0 || Kdx > numCellXYZ[2]-2) )
         continue;

      real xyz000[3]    = { XYZ[0][Idx], XYZ[1][Jdx], XYZ[2][Kdx] };

      real Vertex000[1] = { TargetQuantity[Idx  ][Jdx  ][Kdx  ] };
      real Vertex001[1] = { TargetQuantity[Idx  ][Jdx  ][Kdx+1] };
      real Vertex010[1] = { TargetQuantity[Idx  ][Jdx+1][Kdx  ] };
      real Vertex100[1] = { TargetQuantity[Idx+1][Jdx  ][Kdx  ] };
      real Vertex011[1] = { TargetQuantity[Idx  ][Jdx+1][Kdx+1] };
      real Vertex101[1] = { TargetQuantity[Idx+1][Jdx  ][Kdx+1] };
      real Vertex110[1] = { TargetQuantity[Idx+1][Jdx+1][Kdx  ] };
      real Vertex111[1] = { TargetQuantity[Idx+1][Jdx+1][Kdx+1] };

      real FieldAtVertices[8] = { Vertex000[0], Vertex001[0], Vertex010[0], Vertex100[0],
                                  Vertex011[0], Vertex101[0], Vertex110[0], Vertex111[0] };

      ProjectedValue += TrilinearInterpolation( FieldAtVertices, xyz000, dxyz, xyz );
   } while( 1 );

   return (real)ProjectedValue;
} // FUNCTION : PerspectiveProject



//-------------------------------------------------------------------------------------------------------
// Function    :  rotationMatrix
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void rotationMatrix( real x, real y, real z, real *xp, real *yp, real *zp, real angle )
{
   *xp = x * COS(angle) - y * SIN(angle);
   *yp = x * SIN(angle) + y * COS(angle);
   *zp = z;
} // FUNCTION : rotationMatrix
