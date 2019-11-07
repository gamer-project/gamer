#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"
#include "SRHydroPrototypes.h"

int
main ()
{
  double Gamma = 1.333333;

  double Con[5] = { 0 };

  Con[0] = 0.5;
  Con[1] = 1e-15;
  Con[2] = 1e-18;
  Con[3] = -1e-16;
  Con[4] = 1e-16;

  real discriminant, Msqr;

  Msqr = VectorDotProduct( Con[MOMX], Con[MOMY], Con[MOMZ] );

  discriminant = SQR(Con[ENGY]/Con[DENS]) + (real)2*(Con[ENGY]/Con[DENS]) - Msqr/SQR(Con[DENS]);

  if ( discriminant <= 0.0 )
  {
    printf("discriminant = %e < 0.0 !!\n", discriminant);
    exit(0);
  }

  double Pri[5] = { 0 };

  printf ("Initial cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con[0], Con[1], Con[2], Con[3], Con[4]);

  SRHydro_Con2Pri (Con, Pri, Gamma, 0.0);

  printf ("Transform to prim:\n\n d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);

  double Vx, Vy, Vz;

  real G = sqrt(1.0 + Pri[1]*Pri[1] + Pri[2]*Pri[2] + Pri[3]*Pri[3]);

  Vx = Pri[1] / G;
  Vy = Pri[2] / G;
  Vz = Pri[3] / G;

  printf ("Transform to prim:\n\n d=%e, Vx=%e, Vy=%e, Vz=%e, V=%e, P=%e\n\n", Pri[0], Vx, Vy, Vz, sqrt(Vx*Vx+Vy*Vy+Vz*Vz), Pri[4]);

  double Con_re[5] = { 0 };
  SRHydro_Pri2Con (Pri, Con_re, Gamma);

  printf ("Transform back to cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con_re[0], Con_re[1], Con_re[2], Con_re[3], Con_re[4]);

  if ((fabs (Con_re[1]) > TINY_NUMBER) 
   && (fabs (Con_re[2]) > TINY_NUMBER)
   && (fabs (Con_re[3]) > TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_M1 = (Con_re[1] - Con[1]) / Con[1];
      double err_M2 = (Con_re[2] - Con[2]) / Con[2];
      double err_M3 = (Con_re[3] - Con[3]) / Con[3];
      double err_E = (Con_re[4] - Con[4]) / Con[4];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E err_M3=%E, err_E=%E\n", err_d,  err_M1, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");

    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M3=%E, err_E=%E\n", err_d, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M2 = (Con_re[2] - Con[2]) / Con[2];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_E=%E\n", err_d, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M1 = (Con_re[1] - Con[1]) / Con[1];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_E=%E\n", err_d, err_M1, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M1 = (Con_re[1] - Con[1]) / Con[1];
      double err_M2 = (Con_re[2] - Con[2]) / Con[2];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M2=%E, err_E=%E\n", err_d, err_M1, err_M2, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) > TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M1 = (Con_re[1] - Con[1]) / Con[1];
      double err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M1=%E, err_M3=%E, err_E=%E\n", err_d, err_M1, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) > TINY_NUMBER)
	&& (fabs (Con_re[3]) > TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];
      double err_M2 = (Con_re[2] - Con[2]) / Con[2];
      double err_M3 = (Con_re[3] - Con[3]) / Con[3];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_M2=%E, err_M3=%E, err_E=%E\n", err_d, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_E=%E\n", err_d, err_E);
      printf ("===========================================\n\n");
    }



  return 0;
}
