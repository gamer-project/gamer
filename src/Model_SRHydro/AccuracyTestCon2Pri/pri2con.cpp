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
  real Pri[5] = { 0 };

  Pri[0] = 0.5;
  Pri[1] = 1024;
  Pri[2] = 0.0;
  Pri[3] = 0.0;
  Pri[4] = 2.0;


  printf ("Initial prim:\n\n n=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);

  SRHydro_Pri2Con (Pri, Con, Gamma);

  printf ("Transform to cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con[0], Con[1], Con[2], Con[3], Con[4]);

  real M_Dsqr = VectorDotProduct( Con[1], Con[2], Con[3] ) / SQR(Con[0]);
  real E_D    = Con[4]/Con[0];

  real Discriminant = SQR(E_D) + (real)2.0*E_D - M_Dsqr;

  if ( Discriminant < (real)TINY_NUMBER )
  {
		  printf("ERROR!! Discriminant=%e\n", Discriminant);
		  exit(0);
  }

  real Pri_re[5] = { 0 };
  real LorentzFactor;

  LorentzFactor = SRHydro_Con2Pri (Con, Pri_re, Gamma, 0.0);

  printf ("Transform back to prim:\n\n n=%e, Ux=%e, Uy=%e, Uz=%e, P=%e, gamma=%e\n\n", Pri_re[0], Pri_re[1], Pri_re[2], Pri_re[3], Pri_re[4], LorentzFactor);

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
