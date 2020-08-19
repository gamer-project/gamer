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
  real LorentzFactor;

  real Pri[5] = { 0 };
  real Con[5] = { 0 };

  //Con[0] = 1.0000499987500625138636678e+04;
  //Con[1] = 1.0251637480159439146518707e+06;
  //Con[2] = 0.0000000000000000000000000e+00;
  //Con[3] = 0.0000000000000000000000000e+00;
  //Con[4] = 1.0152135049344534054398537e+06;
  Con[0] = 1.0;
  Con[1] = 1000.0;
  Con[2] = 0.0;
  Con[3] = 0.0;
  Con[4] = 2.9970014999996247e+03;


  printf ("\ninitial conservative vars:\n");
  printf ("=========================================\n\n");

  printf ("D     = %30.25e\n"  ,  Con[0]);
  printf ("Mx    = %30.25e\n"  ,  Con[1]);
  printf ("My    = %30.25e\n"  ,  Con[2]);
  printf ("Mz    = %30.25e\n"  ,  Con[3]);
  printf ("E     = %30.25e\n"  ,  Con[4]);


  real M_Dsqr = VectorDotProduct( Con[1]/Con[0], Con[2]/Con[0], Con[3]/Con[0] );
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

  LorentzFactor = SRHydro_Con2Pri (Con, Pri, Gamma, (real)0.0);

  printf ("\nconservative --> primitive\n");
  printf ("=========================================\n\n");

  printf ("n     = %30.25e\n"  ,  Pri[0]);
  printf ("Ux    = %30.25e\n"  ,  Pri[1]);
  printf ("Uy    = %30.25e\n"  ,  Pri[2]);
  printf ("Uz    = %30.25e\n"  ,  Pri[3]);
  printf ("P     = %30.25e\n"  ,  Pri[4]);
  printf ("T     = %30.25e\n"  ,  Pri[4]/Pri[0]);
  printf ("gamma = %30.25e\n"  ,  LorentzFactor);



  real Con_re[5] = { 0 };

  SRHydro_Pri2Con (Pri, Con_re, Gamma);


  printf ("\nprimitive --> conservative\n");
  printf ("=========================================\n\n");
  
  printf ("D     = %30.25e\n"  ,  Con_re[0]);
  printf ("Mx    = %30.25e\n"  ,  Con_re[1]);
  printf ("My    = %30.25e\n"  ,  Con_re[2]);
  printf ("Mz    = %30.25e\n"  ,  Con_re[3]);
  printf ("E     = %30.25e\n"  ,  Con_re[4]);

  real err_D  = (Con_re[0] - Con[0]) / Con[0];
  real err_M1 = (Con_re[1] - Con[1]) / Con[1];
  real err_M2 = (Con_re[2] - Con[2]) / Con[2];
  real err_M3 = (Con_re[3] - Con[3]) / Con[3];
  real err_E  = (Con_re[4] - Con[4]) / Con[4];

  printf ("\nrelative error:\n");
  printf ("=========================================\n\n");

  printf ("err_D   = %+30.25e\n", err_D );

  if ( fabs (Con_re[1]) > TINY_NUMBER ) printf ("err_Mx  = %+30.25e\n", err_M1) ;
  if ( fabs (Con_re[2]) > TINY_NUMBER ) printf ("err_My  = %+30.25e\n", err_M2) ;
  if ( fabs (Con_re[3]) > TINY_NUMBER ) printf ("err_Mz  = %+30.25e\n", err_M3) ;

  printf ("err_E   = %+30.25e\n", err_E );


  return 0;
}
