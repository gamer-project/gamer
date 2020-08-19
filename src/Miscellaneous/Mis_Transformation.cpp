# include "GAMER.h"

// Reference: https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
//
// r[0] = r
// r[1] = theta
// r[2] = phi

void Cartesian2Spherical( const double x[], double r[] )
{

  r[0] = sqrt( SQR(x[0]) + SQR(x[1]) + SQR(x[2]) );

  if ( SQR(x[0]) + SQR(x[1]) == 0.0 && x[2] == 0.0 )
  {
    printf("error! : %s: %d\n", __FUNCTION__, __LINE__);
	exit(1);
  }

  r[1] = atan2( sqrt( SQR(x[0]) + SQR(x[1]) ), x[2] );

  if ( x[1] == 0.0 && x[0] == 0.0 )
  {
    printf("error! : %s: %d\n", __FUNCTION__, __LINE__);
	exit(1);
  }

  r[2] = atan2( x[1], x[0] );
}


void Cartesian2Cylindrical( const double x[], double R[] )
{
  R[0] = sqrt( SQR(x[0]) + SQR(x[1]) );

  if ( x[1] == 0.0 && x[0] == 0.0 )
  {
    printf("error! : %s: %d\n", __FUNCTION__, __LINE__);
	exit(1);
  }

  R[1] = atan2( x[1], x[0] );

  R[2] = x[2];
}

void Spherical2Cartesian( const double r[], double x[] )
{
  x[0] = r[0]*sin(r[1])*cos(r[2]);

  x[1] = r[0]*sin(r[1])*sin(r[2]);

  x[2] = r[0]*cos(r[1]);
}

void CartesianRotate( double x[], double theta, double phi, bool inverse )
{
  double xp[3];

  if ( inverse )
  {
     xp[0] = -            sin(phi)*x[0] - cos(theta)*cos(phi)*x[1] + sin(theta)*cos(phi)*x[2];
	 xp[1] = +            cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*sin(phi)*x[2];
	 xp[2] =                            + sin(theta)*         x[1] + cos(theta)*         x[2];
  }
  else
  {
     xp[0] = -            sin(phi)*x[0] +            cos(phi)*x[1];
	 xp[1] = - cos(theta)*cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*         x[2];
	 xp[2] = + sin(theta)*cos(phi)*x[0] + sin(theta)*sin(phi)*x[1] + cos(theta)*         x[2];
  }

  for (int i=0;i<3;i++) x[i] = xp[i];
}
