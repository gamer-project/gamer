#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"
#include "SRHydroPrototypes.h"

#define GAMMA 1.33333333333
#define MIN_TEMP 0.0


real H_Tilde ( real T )
{
 real H_Tilde  = (real)2.5*T + (real)2.25*T*T;
      H_Tilde /= (real)1.0 + sqrt( (real)1.0 + (real)2.25*T*T );

 return H_Tilde;
}

int  main()
{
 
  real relative_error[8];
  real Con[5]           ;
  real Pri_re[5]        ;
  real Pri[5]           ;
  real HTilde;
  real HTilde_re;
  real H;
  real H_re;

//  FILE *fptr=fopen("ErrorMap1.dat", "w");
//  real U           = (real)1e6;
//  double Tmax      = 1e4;
//  double Tmin      = 1e-4;  // CMB temperature ~ 2.3e-8

  FILE *fptr=fopen("ErrorMap2.dat", "w");
  real U           = (real)1e-8;
  double Tmax      = 1;
  double Tmin      = 1e-8;  // CMB temperature ~ 2.3e-8

  real rho  = 1.0;
  real Ux   = U/sqrt(3.0);
  real Uy   = U/sqrt(3.0);
  real Uz   = U/sqrt(3.0);
  real LorentzFactor = SQRT((double)1.0 + SQR(Ux) + SQR(Uy) + SQR(Uz));


  int  Npoint    = 1000;

  double Temp = Tmin;

  fprintf(fptr,"#[0] error rho[1] error Ux  [2] error Uy  [3] error Uz  [ 4] error P  [ 5] error Temp [6] error HTilde [7] error H [8] dT\n");
  fprintf(fptr,"#[9] rho   [10] Ux       [11] Uy        [12] Uz        [13] P        [14] LorentzFactor  [15] T [16] HTilde [17] HTilde_re [18] H [19] H_re\n");
  fprintf(fptr,"####################################################\n");

  for ( int i=1 ;i<=Npoint; i++ )
  {
      double dT    = pow( Tmax/Tmin, 1.0/(double)(Npoint-1) );
      Temp *= dT;

      Pri[0] = rho;
      Pri[1] = Ux;
      Pri[2] = Uy; 
      Pri[3] = Uz; 
      Pri[4] = Pri[0]*Temp;
      HTilde = H_Tilde(Temp);
      H      = (real)1.0/((real)1.0 + HTilde);
    
      SRHydro_Pri2Con ( Pri, Con, GAMMA);
      SRHydro_Con2Pri ( Con, Pri_re, GAMMA, MIN_TEMP);

      real Temp_re = Pri_re[4] / Pri_re[0];	  
      HTilde_re = H_Tilde(Temp_re);
      H_re      = (real)1.0/((real)1.0 + HTilde_re);

      relative_error[0] = 1.0 - Pri_re[0] / Pri[0];
      relative_error[1] = 1.0 - Pri_re[1] / Pri[1];
      relative_error[2] = 1.0 - Pri_re[2] / Pri[2];
      relative_error[3] = 1.0 - Pri_re[3] / Pri[3];
      relative_error[4] = 1.0 - Pri_re[4] / Pri[4];
      relative_error[5] = 1.0 - Temp_re   / Temp;
      relative_error[6] = 1.0 - HTilde_re / HTilde;
      relative_error[7] = 1.0 - H_re      / H;
 
      fprintf(fptr,"%+20.16e  ",relative_error[0]);
      fprintf(fptr,"%+20.16e  ",relative_error[1]);
      fprintf(fptr,"%+20.16e  ",relative_error[2]);
      fprintf(fptr,"%+20.16e  ",relative_error[3]);
      fprintf(fptr,"%+20.16e  ",relative_error[4]);
      fprintf(fptr,"%+20.16e  ",relative_error[5]);
      fprintf(fptr,"%+20.16e  ",relative_error[6]);
      fprintf(fptr,"%+20.16e  ",relative_error[7]);
      fprintf(fptr,"%+20.16e  ",dT    );
      fprintf(fptr,"%+20.16e  ",Pri[0]);
      fprintf(fptr,"%+20.16e  ",Pri[1]);
      fprintf(fptr,"%+20.16e  ",Pri[2]);
      fprintf(fptr,"%+20.16e  ",Pri[3]);
      fprintf(fptr,"%+20.16e  ",Pri[4]);
      fprintf(fptr,"%+20.16e  ",LorentzFactor);
      fprintf(fptr,"%+20.16e  ",Temp);
      fprintf(fptr,"%+20.16e  ",HTilde);
      fprintf(fptr,"%+20.16e  ",HTilde_re);
      fprintf(fptr,"%+20.16e  ",H);
      fprintf(fptr,"%+20.16e  ",H_re);
      fprintf(fptr,"\n");
  }



return 0;
}
