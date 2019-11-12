#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"
#include "SRHydroPrototypes.h"

#ifdef FLOAT8
#define real double
#else
#define real float
#endif

int
main ()
{
  real Gamma = 1.333333;

  real Con[5] = { 0 };

  Con[0] = 5.000250e+01;
  Con[1] = 8.011784e+04;
  Con[2] = 0.0;
  Con[3] = 0.0;
  Con[4] = 8.036986e+04;

  real discriminant, Msqr;

  Msqr = VectorDotProduct( Con[MOMX], Con[MOMY], Con[MOMZ] );

  discriminant = SQR(Con[ENGY]/Con[DENS]) + (real)2*(Con[ENGY]/Con[DENS]) - Msqr/SQR(Con[DENS]);

  if ( discriminant <= 0.0 )
  {
    printf("discriminant = %e < 0.0 !!\n", discriminant);
    exit(0);
  }

  real Pri[5] = { 0 };

  printf ("Initial cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con[0], Con[1], Con[2], Con[3], Con[4]);

  SRHydro_Con2Pri (Con, Pri, Gamma, 0.0);

  printf ("Transform to prim:\n\n d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);

  real Vx, Vy, Vz;

  real G = sqrt(1.0 + Pri[1]*Pri[1] + Pri[2]*Pri[2] + Pri[3]*Pri[3]);

  Vx = Pri[1] / G;
  Vy = Pri[2] / G;
  Vz = Pri[3] / G;

  //printf ("Transform to prim:\n\n d=%e, Vx=%e, Vy=%e, Vz=%e, V=%e, P=%e\n\n", Pri[0], Vx, Vy, Vz, sqrt(Vx*Vx+Vy*Vy+Vz*Vz), Pri[4]);

  real Con_re[5] = { 0 };
  SRHydro_Pri2Con (Pri, Con_re, Gamma);

  printf ("Transform back to cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con_re[0], Con_re[1], Con_re[2], Con_re[3], Con_re[4]);
  if ((fabs (Con_re[1]) > TINY_NUMBER) 
   && (fabs (Con_re[2]) > TINY_NUMBER)
   && (fabs (Con_re[3]) > TINY_NUMBER))
    {
      real err_d  = (Con_re[0] - Con[0]) / Con[0];
      real err_M1 = (Con_re[1] - Con[1]) / Con[1];
      real err_M2 = (Con_re[2] - Con[2]) / Con[2];
      real err_M3 = (Con_re[3] - Con[3]) / Con[3];
      real err_E  = (Con_re[4] - Con[4]) / Con[4];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E err_M3=%E, err_E=%E\n", err_d,  err_M1, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");

    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M3=%E, err_E=%E\n", err_d, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M2 = (Con_re[2] - Con[2]) / Con[2];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_E=%E\n", err_d, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M1 = (Con_re[1] - Con[1]) / Con[1];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_E=%E\n", err_d, err_M1, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M1 = (Con_re[1] - Con[1]) / Con[1];
      real err_M2 = (Con_re[2] - Con[2]) / Con[2];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E, err_E=%E\n", err_d, err_M1, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M1 = (Con_re[1] - Con[1]) / Con[1];
      real err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M3=%E, err_E=%E\n", err_d, err_M1, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];
      real err_M2 = (Con_re[2] - Con[2]) / Con[2];
      real err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_M3=%E, err_E=%E\n", err_d, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      real err_d = (Con_re[0] - Con[0]) / Con[0];
      real err_E = (Con_re[4] - Con[4]) / Con[4];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_E=%E\n", err_d, err_E);
      printf ("===========================================\n\n");
    }



  real Pri_re[5] = { 0 };
  SRHydro_Con2Pri (Con_re, Pri_re, Gamma, 0.0);
  printf ("Transform to prim:\n\n d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri_re[0], Pri_re[1], Pri_re[2], Pri_re[3], Pri_re[4]);


  if ((fabs (Pri_re[1]) > TINY_NUMBER) 
   && (fabs (Pri_re[2]) > TINY_NUMBER)
   && (fabs (Pri_re[3]) > TINY_NUMBER))
    {
      real err_n  = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_U1 = (Pri_re[1] - Pri[1]) / Pri[1];
      real err_U2 = (Pri_re[2] - Pri[2]) / Pri[2];
      real err_U3 = (Pri_re[3] - Pri[3]) / Pri[3];
      real err_P  = (Pri_re[4] - Pri[4]) / Pri[4];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U1=%E, err_U2=%E err_U3=%E, err_P=%E\n", err_n,  err_U1, err_U2, err_U3, err_P);
      printf ("===========================================\n\n");

    }
  else if ((fabs (Pri_re[1]) < TINY_NUMBER) 
	&& (fabs (Pri_re[2]) < TINY_NUMBER)
	&& (fabs (Pri_re[3]) > TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U3 = (Pri_re[3] - Pri[3]) / Pri[3];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U3=%E, err_P=%E\n", err_n, err_U3, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) < TINY_NUMBER) 
	&& (fabs (Pri_re[2]) > TINY_NUMBER)
	&& (fabs (Pri_re[3]) < TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U2 = (Pri_re[2] - Pri[2]) / Pri[2];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U2=%E, err_P=%E\n", err_n, err_U2, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) > TINY_NUMBER) 
	&& (fabs (Pri_re[2]) < TINY_NUMBER)
	&& (fabs (Pri_re[3]) < TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U1 = (Pri_re[1] - Pri[1]) / Pri[1];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U1=%E, err_P=%E\n", err_n, err_U1, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) > TINY_NUMBER) 
	&& (fabs (Pri_re[2]) > TINY_NUMBER)
	&& (fabs (Pri_re[3]) < TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U1 = (Pri_re[1] - Pri[1]) / Pri[1];
      real err_U2 = (Pri_re[2] - Pri[2]) / Pri[2];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U1=%E, err_U2=%E, err_P=%E\n", err_n, err_U1, err_U2, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) > TINY_NUMBER) 
	&& (fabs (Pri_re[2]) < TINY_NUMBER)
	&& (fabs (Pri_re[3]) > TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U1 = (Pri_re[1] - Pri[1]) / Pri[1];
      real err_U3 = (Pri_re[3] - Pri[3]) / Pri[3];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U1=%E, err_U3=%E, err_P=%E\n", err_n, err_U1, err_U3, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) < TINY_NUMBER) 
	&& (fabs (Pri_re[2]) > TINY_NUMBER)
	&& (fabs (Pri_re[3]) > TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];
      real err_U2 = (Pri_re[2] - Pri[2]) / Pri[2];
      real err_U3 = (Pri_re[3] - Pri[3]) / Pri[3];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_U2=%E, err_U3=%E, err_P=%E\n", err_n, err_U2, err_U3, err_P);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Pri_re[1]) < TINY_NUMBER) 
	&& (fabs (Pri_re[2]) < TINY_NUMBER)
	&& (fabs (Pri_re[3]) < TINY_NUMBER))
    {
      real err_n = (Pri_re[0] - Pri[0]) / Pri[0];
      real err_P = (Pri_re[4] - Pri[4]) / Pri[4];

      printf ("relative error:\n\n");
      printf ("err_n=%E, err_P=%E\n", err_n, err_P);
      printf ("===========================================\n\n");
    }



  return 0;
}
