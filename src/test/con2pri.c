#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define SQR(X) ((X)*(X))
#define TINY_NUMBER 1.0e-100

struct FUN_Q_params
{
  double d;
  double M1;
  double M2;
  double M3;
  double E;
  double Gamma;
};

// some functions in this file need to be defined even when using GPU

static double FUN_Q (double Q, void *ptr);
static double D_FUN_Q (double Q, void *ptr);
static void FDF_FUN_Q (double Q, void *ptr, double *f, double *df);
void CPU_Pri2Con (const double In[], double Out[], const double Gamma);
void CPU_Con2Pri (const double In[], double Out[], const double Gamma);

int
main ()
{
  double Gamma = 1.3333333;

  double Con[5] = { 1.0, 0.5, 0.1, 0.8, 5 };
  double Pri[5] = { 0 };

  printf ("Initial cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con[0], Con[1], Con[2], Con[3], Con[4]);

  CPU_Con2Pri (Con, Pri, Gamma);

  printf ("Transform to prim:\n\n d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);

  double Con_re[5] = { 0 };
  CPU_Pri2Con (Pri, Con_re, Gamma);

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
      printf ("err_d=%E, err_M1=%E, err_M2=%E err_M3=%E, err_P=%E\n", err_d,  err_M1, err_M2, err_M3, err_E);
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
      printf ("err_d=%E, err_M3=%E, err_P=%E\n", err_d, err_M3, err_E);
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
      printf ("err_d=%E, err_M2=%E, err_P=%E\n", err_d, err_M2, err_E);
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
      printf ("err_d=%E, err_M1=%E, err_P=%E\n", err_d, err_M1, err_E);
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
      printf ("err_d=%E, err_M1=%E, err_M2=%E, err_P=%E\n", err_d, err_M1, err_M2, err_E);
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
      printf ("err_d=%E, err_M1=%E, err_M3=%E, err_P=%E\n", err_d, err_M1, err_M3, err_E);
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
      printf ("err_d=%E, err_M2=%E, err_M3=%E, err_P=%E\n", err_d, err_M2, err_M3, err_E);
      printf ("===========================================\n\n");
    }
  else if ((fabs (Con_re[1]) < TINY_NUMBER) 
	&& (fabs (Con_re[2]) < TINY_NUMBER)
	&& (fabs (Con_re[3]) < TINY_NUMBER))
    {
      double err_d = (Con_re[0] - Con[0]) / Con[0];
      double err_E = (Con_re[4] - Con[4]) / Con[4];

      printf ("relative error:\n\n");
      printf ("err_d=%E, err_P=%E\n", err_d, err_E);
      printf ("===========================================\n\n");
    }



  return 0;
}



void
CPU_Pri2Con (const double In[], double Out[], const double Gamma)
{
  double Gamma_m1 = (double) Gamma - 1.0;

  double Factor0 = 1 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
  double Factor1 = sqrt (Factor0);
  double Factor2 = (Gamma / Gamma_m1);
  double Factor3 = In[0] + Factor2 * In[4];
  double Factor4 = Factor3 * Factor1;

  Out[0] = In[0] * Factor1;
  Out[1] = Factor4 * In[1];
  Out[2] = Factor4 * In[2];
  Out[3] = Factor4 * In[3];
  Out[4] = Factor3 * Factor0 - In[4];

}

			       // FUNCTION : CPU_Pri2Con
void
CPU_Con2Pri (const double In[], double Out[], const double Gamma)
{
#ifdef SR_DEBUG
  double In_temp[5] = { In[0], In[1], In[2], In[3], In[4] };
  if ((In[0] < 0) || (In[4] < 0))
    {
      printf ("\n\nerror: D < 0 or E < 0!\n");
      printf ("file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      printf ("line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In[0], In[1], In[2], In[3], In[4]);
      abort ();
    }
#endif
  double Gamma_m1 = Gamma - (double) 1.0;
  double M = sqrt (SQR (In[1]) + SQR (In[2]) + SQR (In[3]));
#ifdef SR_DEBUG
  if (M > In[4])
    {
      printf ("\n\nerror: |M| > E!\n");
      printf ("file: %s\nfunction: %s\n", __FILE__, __FUNCTION__);
      printf ("line:%d\nD=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
      printf ("|M|=%e, E=%e, |M|-E=%e\n\n", M, In[4], M - In[4]);
      abort ();
    }
#endif
  if (fabs (M) > TINY_NUMBER)
    {
      int status;

      int iter = 0;
      int max_iter = 200;

      const gsl_root_fdfsolver_type *T;

      gsl_root_fdfsolver *s;

      double Q, Q0;

/* initial guess Q  */
      if (In[0] > M / Gamma_m1)
	{
	  Q = M * (In[4] - M) / ((1 - 1 / Gamma) * In[0]);
	}
      else
	{
	  Q = In[4] * Gamma;
	}

      printf("Q = %14.7e\n", Q);

      gsl_function_fdf F;

      struct FUN_Q_params params = { In[0], In[1], In[2], In[3], In[4], Gamma };

      F.f = &FUN_Q;
      F.df = &D_FUN_Q;
      F.fdf = &FDF_FUN_Q;
      F.params = &params;

      T = gsl_root_fdfsolver_newton;
      s = gsl_root_fdfsolver_alloc (T);
      gsl_root_fdfsolver_set (s, &F, Q);

      //printf ("status = %s\n", gsl_strerror (status)); 
      do
	{
	  iter++;
	  status = gsl_root_fdfsolver_iterate (s);
	  //printf ("status = %s\n", gsl_strerror (status));
	  Q0 = Q;
	  Q = gsl_root_fdfsolver_root (s);
	  status = gsl_root_test_delta (Q, Q0, 0, 1e-9);
	  //printf ("status = %s\n", gsl_strerror (status));
	}
      while (status == GSL_CONTINUE && iter < max_iter);
      //printf ("status = %s\n", gsl_strerror (status));

      Out[1] = In[1] / Q;	/*Ux */
      Out[2] = In[2] / Q;	/*Uy */
      Out[3] = In[3] / Q;	/*Uz */

      double Factor = sqrt (1 + SQR (Out[1]) + SQR (Out[2]) + SQR (Out[3]));

      Out[0] = In[0] / Factor;	/*primitive density */
      Out[4] = (Gamma_m1 / Gamma) * fabs (Q / Factor - Out[0]);	/*pressure */

#ifdef SR_DEBUG
      if (Out[4] < 0)
	{
	  printf ("\n\nerror: P < 0!\n");
	  printf ("file: %s\nfunction: %s\nline:%d \n", __FILE__, __FUNCTION__, __LINE__);
	  printf ("Q = %e\n", Q);
	  printf ("D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
	  printf ("|M| > TINY_NUMBER!!\n");
	  abort ();
	}
#endif
      gsl_root_fdfsolver_free (s);
    }
  else if (In[4] >= In[0])
    {
      Out[1] = 0.0;
      Out[2] = 0.0;
      Out[3] = 0.0;
      Out[0] = In[0];
      Out[4] = Gamma_m1 * (In[4] - In[0]);
    }
  else
    {
      printf ("Too criticle to solve! Somthing went wrong!\nPlease turn debug mode on!\n");
      printf ("line:%d, D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n", __LINE__, In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4]);
      abort ();
    }
}

				// FUNCTION : CPU_Con2Pri
static double
FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);
  double Gamma = (params->Gamma);
  double Gamma_m1 = (double) Gamma - 1.0;
/*4-Velocity*/
  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double rho = d / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3));

  double pres = (Gamma_m1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  double f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;

  return f;
}

static double
D_FUN_Q (double Q, void *ptr)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double Gamma = (params->Gamma);
  double Gamma_m1 = (double) Gamma - 1.0;

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double dd = -d * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_m1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3))
				    - Q * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - dd);

  double df = sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;

  return df;
}

static void
FDF_FUN_Q (double Q, void *ptr, double *f, double *df)
{
  struct FUN_Q_params *params = (struct FUN_Q_params *) ptr;

  double d = (params->d);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E = (params->E);
  double Gamma = (params->Gamma);
  double Gamma_m1 = (double) Gamma - 1.0;

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double rho = d / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3));
  double pres = (Gamma_m1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double drho = -d * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_m1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3))
				    - Q * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - drho);

  *f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;
  *df = sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;
}
