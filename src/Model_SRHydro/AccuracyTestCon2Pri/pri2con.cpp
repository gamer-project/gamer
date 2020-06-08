#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"
#include "SRHydroPrototypes.h"


#ifdef FLOAT8
typedef double real;
#else
typedef float real;
#endif

int
main ()
{
  real Gamma = 1.333333;

  real LorentzFactor_re,Temperature_re;
  real LorentzFactor   ,Temperature;

  real Con[5] = { 0 };
  real Pri[5] = { 0 };

  Pri[0] = 1.0;
  Pri[1] = 10.0;
  Pri[2] = 0.0;
  Pri[3] = 0.0;
  Pri[4] = 1e-4;

  LorentzFactor = SQRT( (real)1.0 + SQR(Pri[1]) + SQR(Pri[2]) + SQR(Pri[3]) );

  printf ("\ninitial primitive vars:\n");
  printf ("=========================================\n\n");

  printf ("n     = %30.25e\n"  ,  Pri[0]);
  printf ("Ux    = %30.25e\n"  ,  Pri[1]);
  printf ("Uy    = %30.25e\n"  ,  Pri[2]);
  printf ("Uz    = %30.25e\n"  ,  Pri[3]);
  printf ("P     = %30.25e\n"  ,  Pri[4]);
  printf ("T     = %30.25e\n"  ,  Pri[4]/Pri[0]);
  printf ("gamma = %30.25e\n"  ,  LorentzFactor);

  SRHydro_Pri2Con (Pri, Con, Gamma);

  printf ("\nprimitive --> conservative\n");
  printf ("=========================================\n\n");

  printf ("D     = %30.25e\n"  ,  Con[0]);
  printf ("Mx    = %30.25e\n"  ,  Con[1]);
  printf ("My    = %30.25e\n"  ,  Con[2]);
  printf ("Mz    = %30.25e\n"  ,  Con[3]);
  printf ("E     = %30.25e\n"  ,  Con[4]);


  real M_Dsqr = VectorDotProduct( Con[1], Con[2], Con[3] ) / SQR(Con[0]);
  real M_D    = SQRT( M_Dsqr );
  real E_D    = Con[4]/Con[0];

  real Discriminant = ( E_D + M_D )*( E_D - M_D ) + (real)2.0*E_D;

  if ( Discriminant < (real)TINY_NUMBER )
  {
		  printf("\nERROR!! Discriminant = %30.25e\n", Discriminant);
		  exit(0);
  }
  else
  {
		  printf("\nDiscriminant = %30.25e\n", Discriminant);
  }

  real Pri_re[5] = { 0 };

  LorentzFactor_re = SRHydro_Con2Pri (Con, Pri_re, Gamma, 0.0);
  Temperature      = Pri[4]   / Pri[0];
  Temperature_re   = Pri_re[4]/ Pri_re[0];

  printf ("\nconservative --> primitive\n");
  printf ("=========================================\n\n");
  
  printf ("n     = %30.25e\n"  ,  Pri_re[0]);
  printf ("Ux    = %30.25e\n"  ,  Pri_re[1]);
  printf ("Uy    = %30.25e\n"  ,  Pri_re[2]);
  printf ("Uz    = %30.25e\n"  ,  Pri_re[3]);
  printf ("P     = %30.25e\n"  ,  Pri_re[4]);
  printf ("T     = %30.25e\n"  ,  Temperature_re);
  printf ("gamma = %30.25e\n"  ,  LorentzFactor_re);

  real err_n  = (Pri_re[0] - Pri[0]) / Pri[0];
  real err_U1 = (Pri_re[1] - Pri[1]) / Pri[1];
  real err_U2 = (Pri_re[2] - Pri[2]) / Pri[2];
  real err_U3 = (Pri_re[3] - Pri[3]) / Pri[3];
  real err_P  = (Pri_re[4] - Pri[4]) / Pri[4];
  real err_T  = ( Temperature_re - Temperature ) / Temperature;
  real err_g  = ( LorentzFactor_re - LorentzFactor ) / LorentzFactor;

  printf ("\nrelative error:\n");
  printf ("=========================================\n\n");

  printf ("err_n   = %+30.25e\n", err_n );

  if ( fabs (Pri_re[1]) > TINY_NUMBER ) printf ("err_Ux  = %+30.25e\n", err_U1) ;
  if ( fabs (Pri_re[2]) > TINY_NUMBER ) printf ("err_Uy  = %+30.25e\n", err_U2) ;
  if ( fabs (Pri_re[3]) > TINY_NUMBER ) printf ("err_Uz  = %+30.25e\n", err_U3) ;

  printf ("err_P   = %+30.25e\n", err_P );
  printf ("err_T   = %+30.25e\n", err_T );
  printf ("err_g   = %+30.25e\n", err_g );

//===============================================================
  printf ("\nprimitive --> conservative\n");
  printf ("=========================================\n\n");

  real Con_re[5] = { 0.0 };

  SRHydro_Pri2Con( Pri_re, Con_re, Gamma );

  printf ("D     = %30.25e\n"  ,  Con_re[0]);
  printf ("Mx    = %30.25e\n"  ,  Con_re[1]);
  printf ("My    = %30.25e\n"  ,  Con_re[2]);
  printf ("Mz    = %30.25e\n"  ,  Con_re[3]);
  printf ("E     = %30.25e\n"  ,  Con_re[4]);


  real err_D  = (Con_re[0] - Con[0]) / Con[0];
  real err_Mx = (Con_re[1] - Con[1]) / Con[1];
  real err_My = (Con_re[2] - Con[2]) / Con[2];
  real err_Mz = (Con_re[3] - Con[3]) / Con[3];
  real err_E  = (Con_re[4] - Con[4]) / Con[4];

  printf ("\nrelative error:\n");
  printf ("=========================================\n\n");

  printf ("err_D   = %+30.25e\n", err_D );

  if ( fabs (Con_re[1]) > TINY_NUMBER ) printf ("err_Mx  = %+30.25e\n", err_Mx) ;
  if ( fabs (Con_re[2]) > TINY_NUMBER ) printf ("err_My  = %+30.25e\n", err_My) ;
  if ( fabs (Con_re[3]) > TINY_NUMBER ) printf ("err_Mz  = %+30.25e\n", err_Mz) ;

  printf ("err_E   = %+30.25e\n", err_E );
    

//===============================================================


  return 0;
}
