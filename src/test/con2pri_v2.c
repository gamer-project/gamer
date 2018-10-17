#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define SQR(X) ((X)*(X))
#define TINY_NUMBER 1.0e-100
#define NCOMP_FLUID 5

struct FUN_params
{
  double D;
  double M1;
  double M2;
  double M3;
  double E;
};

  double Gamma = 1.3333333;
// some functions in this file need to be defined even when using GPU

static double FUN (double Q, void *ptr);
static double D_FUN (double Q, void *ptr);
static void FDF_FUN (double Q, void *ptr, double *f, double *df);
void CPU_Pri2Con (const double In[], double Out[], const double Gamma);
void CPU_Con2Pri (const double In[], double Out[], const double Gamma);

int
main ()
{

  double Con[5] = { 1.0, -6.0, 2.0, 10.0, 500 };
  double Pri[5] = { 0 };

  printf ("Initial cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con[0], Con[1], Con[2], Con[3], Con[4]);

  CPU_Con2Pri (Con, Pri, Gamma);

  printf ("Transform to prim:\n\n d=%e, Ux=%e, Uy=%e, Uz=%e, P=%e\n\n", Pri[0], Pri[1], Pri[2], Pri[3], Pri[4]);

  double Con_re[5] = { 0 };
  CPU_Pri2Con (Pri, Con_re, Gamma);

  printf ("Transform back to cons:\n\n D=%e, Mx=%e, My=%e, Mz=%e, E=%e\n\n", Con_re[0], Con_re[1], Con_re[2], Con_re[3], Con_re[4]);

  double err_d  = (Con_re[0] - Con[0]) / Con[0];
  double err_M1 = (Con_re[1] - Con[1]) / Con[1];
  double err_M2 = (Con_re[2] - Con[2]) / Con[2];
  double err_M3 = (Con_re[3] - Con[3]) / Con[3];
  double err_E  = (Con_re[4] - Con[4]) / Con[4];

  printf ("relative error:\n\n");
  printf ("err_D=%E, err_M1=%E, err_M2=%E err_M3=%E, err_E=%E\n", err_d,  err_M1, err_M2, err_M3, err_E);

  return 0;
}
void
CPU_Con2Pri (const double In[], double Out[], const double Gamma)
{
  double In_temp[5] = { In[0], In[1], In[2], In[3], In[4] };
  double Msqr = SQR (In_temp[1]) + SQR (In_temp[2]) + SQR (In_temp[3]);
  double M = sqrt (Msqr); // magnitude of momentum

      int status;

      int iter = 0;
      int max_iter = 200;

      const gsl_root_fdfsolver_type *T;

      gsl_root_fdfsolver *s;

      double Q, Q0;
# ifdef RELATIVISTIC_EOS
# if   ( CONSERVED_ENERGY == 1 )
/* initial guess Q  */
	   double Constant = SQR(In[4]/In[0]) - Msqr/SQR(In[0]) - 1.0;

	   if ( 1.0 - Msqr/(16*In[0]*In[0]) >= 0 )
             {
	       if ( Constant > 1.5 ) 
                 {
		    Q = sqrt( SQR(In[4]/In[0]) - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else Q = Constant / 3.0;
	     }
	   else // 1 - (M/D)**2 < 0
             {
	       if ( Constant >  0.5 + sqrt( Msqr/(16*SQR(In[0])) - 3.0/4.0 )) 
                 {
	            Q = sqrt( SQR(In[4]/In[0]) - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else Q = Constant / 3.0;
	     }
# elif ( CONSERVED_ENERGY == 2 )
/* initial guess Q  */
	   double Constant = SQR(In[4]/In[0]) + 2*In[4]/In[0] - Msqr/SQR(In[0]);

	   if ( 1.0 - Msqr/(16*In[0]*In[0]) >= 0 )
             {
	       if ( Constant > 1.5 ) 
                 {
	            Q = sqrt( SQR(In[4]/In[0]) + 2*In[4]/In[0] - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else Q = Constant / 3.0;
	     }
	   else // 1 - (M/D)**2 < 0
             {
	       if ( Constant >  0.5 + sqrt( Msqr/(16*SQR(In[0])) - 3.0/4.0 )) 
                 {
		    Q = sqrt( SQR(In[4]/In[0]) + 2*In[4]/In[0] - 0.9375*Msqr/SQR(In[0]) - 1.0 ) / 3.0;
		 }
	       else Q = Constant / 3.0;
	     }
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif

# elif defined IDEAL_GAS_EOS
  double Gamma_m1 = Gamma - (double) 1.0;

/* initial guess Q  */
    if ( M > TINY_NUMBER ) 
    {
	if (In_temp[0] > M / Gamma_m1)  Q = M * (In_temp[4] - M) / ((1 - 1 / Gamma) * In_temp[0]);
        else                            Q = In_temp[4] * Gamma;
    } 
    else                                Q = In_temp[4] * Gamma - In_temp[0] * (Gamma_m1);
#else
#error: unsupported EoS!
#endif
      gsl_function_fdf F;

      struct FUN_params params = { In_temp[0], In_temp[1], In_temp[2], In_temp[3], In_temp[4] };

      F.f = &FUN;
      F.df = &D_FUN;
      F.fdf = &FDF_FUN;
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
#         ifdef FLOAT8
	  status = gsl_root_test_delta (Q, Q0, 0, 1e-16);
#         else
	  status = gsl_root_test_delta (Q, Q0, 0, 1e-7);
#         endif
	  //printf ("status = %s\n", gsl_strerror (status));
	}
      while (status == GSL_CONTINUE && iter < max_iter);
      //printf ("status = %s\n", gsl_strerror (status));
#ifdef RELATIVISTIC_EOS
      double h = 2.5*Q + sqrt(2.25*Q*Q + 1.0);

      Out[1] = In_temp[1]/(In_temp[0]*h);
      Out[2] = In_temp[2]/(In_temp[0]*h);
      Out[3] = In_temp[3]/(In_temp[0]*h);

      double Factor = sqrt(1 + SQR (Out[1]) + SQR (Out[2]) + SQR (Out[3]));

      Out[0] = In_temp[0] / Factor;

#     if     ( CONSERVED_ENERGY == 1 )
      Out[4] = In_temp[0] * h * Factor - In_temp[4];
#     elif   ( CONSERVED_ENERGY == 2 )
      Out[4] = In_temp[0] * h * Factor - In_temp[4] - In_temp[0];
#     else
#     error: CONSERVED_ENERGY must be 1 or 2!
#     endif

#elif defined IDEAL_GAS_EOS 
      Out[1] = In_temp[1] / Q;	/*Ux */
      Out[2] = In_temp[2] / Q;	/*Uy */
      Out[3] = In_temp[3] / Q;	/*Uz */

      double Factor = sqrt (1 + SQR (Out[1]) + SQR (Out[2]) + SQR (Out[3]));

      Out[0] = In_temp[0] / Factor;	/*number density in local rest frame*/
      Out[4] = (Gamma_m1 / Gamma) * fabs (Q / Factor - Out[0]);	/*pressure */

#else
#error: unsupported EoS!
#endif
      gsl_root_fdfsolver_free (s);
}				// FUNCTION : CPU_Con2Pri
void
CPU_Pri2Con (const double In[], double Out[], const double Gamma)
{
#ifdef RELATIVISTIC_EOS
  double Temperature;
  if ( In[4] > TINY_NUMBER )  {  
      Temperature = In[4]/In[0]; // T = number_density/pressure
   } else {
      Temperature = 0.0;         // pressure = number_density * T
   }

  double Enthalpy = 2.5*Temperature + sqrt(2.25*SQR(Temperature)+1.0); // approximate enthalpy
  double Factor0 = sqrt(1.0 + SQR (In[1]) + SQR (In[2]) + SQR (In[3])); // Lorentz factor
  double Factor1 = In[0] * Factor0 * Enthalpy;
  
  Out[0] = In[0] * Factor0; // number density in inertial frame
  Out[1] = Factor1 * In[1]; // MomX
  Out[2] = Factor1 * In[2]; // MomX
  Out[3] = Factor1 * In[3]; // MomX
# if   ( CONSERVED_ENERGY == 1 )
  Out[4] = Factor1 * Factor0 - In[4]; // total_energy
# elif ( CONSERVED_ENERGY == 2 )
  Out[4] = Factor1 * Factor0 - In[4] - Out[0]; // ( total_energy ) - ( rest_mass_energy )
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif

#elif defined IDEAL_GAS_EOS
  double Gamma_m1 = (double) Gamma - 1.0;
  double U = sqrt(SQR(In[1])+SQR(In[2])+SQR(In[3]));

    double Factor0 = 1 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
    double Factor1 = sqrt (Factor0);
    double Factor2 = Gamma / Gamma_m1;
    double Factor3 = In[0] + Factor2 * In[4]; // enthalpy * rho
    double Factor4 = Factor3 * Factor1;

    Out[0] = In[0] * Factor1;
    Out[1] = Factor4 * In[1];
    Out[2] = Factor4 * In[2];
    Out[3] = Factor4 * In[3];
    Out[4] = Factor3 * Factor0 - In[4];

#else
#error: unsupported EoS!
#endif
}				// FUNCTION : CPU_Pri2Con
static double
FUN (double Q, void *ptr)
{
  struct FUN_params *params = (struct FUN_params *) ptr;

  double D  = (params->D);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E  = (params->E);

#ifdef RELATIVISTIC_EOS
  double h = 2.5*Q+sqrt(2.25*SQR(Q)+1.0); // approximate enthalpy
  double Msqr = SQR(M1) + SQR(M2) + SQR(M3);

#if   (CONSERVED_ENERGY == 1)
  double Constant = SQR(E/D) - Msqr/(D*D) - 1.0;
#elif (CONSERVED_ENERGY == 2)
  double Constant = SQR(E/D) + 2*(E/D) - Msqr/(D*D);
# else
# error: CONSERVED_ENERGY must be 1 or 2!
#endif
  double f = SQR(h) - 1.0 - 2*h*Q + (SQR(h*Q)) / (SQR(h)+Msqr/SQR(D)) - Constant;
#elif defined IDEAL_GAS_EOS

  double Gamma_m1 = (double) Gamma - 1.0;
/*4-Velocity*/
  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double rho = D / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3));

  double pres = (Gamma_m1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  double f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;
#else
#error: unsupported EoS!
#endif // #ifdef RELATIVISTIC_EOS

  return f;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

static double
D_FUN (double Q, void *ptr)
{
  struct FUN_params *params = (struct FUN_params *) ptr;

  double D  = (params->D);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);

#ifdef RELATIVISTIC_EOS
  double h = 2.5*Q+sqrt(2.25*SQR(Q)+1.0); // approximate enthalpy
  double dh = 2.5 + 9.0*Q / sqrt(36*Q*Q+16);
  double Msqr = SQR(M1) + SQR(M2) + SQR(M3);

  double df = 2*dh*(h-Q) - 2*h + 2*h*Q*((h*h*h + Msqr*h/(D*D) + Q*dh*Msqr/(D*D)) / SQR(h*h+Msqr/(D*D)) );

#elif defined IDEAL_GAS_EOS

  double Gamma_m1 = (double) Gamma - 1.0;

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double dd = -D * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_m1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3))
				  - Q * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - dd);

  double df = sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;
#else
#error: unsupported EoS!
#endif // #ifdef RELATIVISTIC_EOS

  return df;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  
// Description :  
//-------------------------------------------------------------------------------------------------------

static void
FDF_FUN (double Q, void *ptr, double * f, double * df)
{
  struct FUN_params *params = (struct FUN_params *) ptr;

  double D  = (params->D);
  double M1 = (params->M1);
  double M2 = (params->M2);
  double M3 = (params->M3);
  double E  = (params->E);

#ifdef RELATIVISTIC_EOS
  double h = 2.5*Q+sqrt(2.25*SQR(Q)+1.0); // approximate enthalpy
  double dh = 2.5 + 9.0*Q / sqrt(36*Q*Q+16);
  double Msqr = SQR(M1) + SQR(M2) + SQR(M3);
#if   (CONSERVED_ENERGY == 1)
  double Constant = SQR(E/D) - Msqr/(D*D) - 1.0;
#elif (CONSERVED_ENERGY == 2)
  double Constant = SQR(E/D) + 2*(E/D) - Msqr/(D*D);
# else
# error: CONSERVED_ENERGY must be 1 or 2!
#endif
  *f = SQR(h) - 1.0 - 2*h*Q + (SQR(h*Q)) / (SQR(h)+Msqr/SQR(D)) - Constant;
  *df = 2*dh*(h-Q) - 2*h + 2*h*Q*((h*h*h + Msqr*h/(D*D) + Q*dh*Msqr/(D*D)) / SQR(h*h+Msqr/(D*D)) );

#elif defined IDEAL_GAS_EOS

  double Gamma_m1 = (double) Gamma - 1.0;

  double U1 = M1 / Q;
  double U2 = M2 / Q;
  double U3 = M3 / Q;

  double rho = D / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3));
  double pres = (Gamma_m1 / Gamma) * (Q / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - rho);

  double dU1 = -M1 / (Q * Q);
  double dU2 = -M2 / (Q * Q);
  double dU3 = -M3 / (Q * Q);

  double drho = -D * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3);

  double dp = (Gamma_m1 / Gamma) * (1 / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3))
				  - Q * pow (1 + SQR (U1) + SQR (U2) + SQR (U3), -1.5) * (U1 * dU1 + U2 * dU2 + U3 * dU3) - drho);

  *f = Q * sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - pres - E;
  *df = sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) + Q * (U1 * dU1 + U2 * dU2 + U3 * dU3) / sqrt (1 + SQR (U1) + SQR (U2) + SQR (U3)) - dp;
#else
#error: unsupported EoS!
#endif // #ifdef RELATIVISTIC_EOS
}
