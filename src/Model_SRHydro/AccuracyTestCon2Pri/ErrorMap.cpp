#include <stdio.h>
#include <math.h>
#include "GAMER.h"
#include "CUFLU.h"
#include "SRHydroPrototypes.h"

#define GAMMA 1.33333333333
#define MIN_TEMP 0.0

int  main()
{
  FILE *fptr=fopen("ErrorMap.dat", "w");
 
  real relative_error[6];
  real Con[5]           ;
  real Pri_re[5]        ;
  real Pri[5]           ;


  real rho  = 1e2;
  real Ux   = 1e6;
  real Uy   = 0.0;
  real Uz   = 0.0;
  real LorentzFactor = SQRT((double)1.0 + SQR(Ux) + SQR(Uy) + SQR(Uz));

  double Tmax      = 1e6  ;
  double Tmin      = 1e-18;

  int  Npoint    = 1000;

  double Temp = Tmin;

  fprintf(fptr,"#[0] error rho[1] error Ux  [2] error Uy  [3] error Uz  [ 4] error P  [ 5] error Temp [ 6] dT\n");
  fprintf(fptr,"#[7] rho   [8] Ux       [9] Uy        [10] Uz        [11] P        [12] LorentzFactor  [13] T\n");
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
    
      SRHydro_Pri2Con ( Pri, Con, GAMMA);
      SRHydro_Con2Pri ( Con, Pri_re, GAMMA, MIN_TEMP);

      real Temp_re = Pri_re[4] / Pri_re[0];	  

      relative_error[0] = 1.0 - Pri_re[0] / Pri[0];
      relative_error[1] = 1.0 - Pri_re[1] / Pri[1];
      relative_error[2] = 1.0 - Pri_re[2] / Pri[2];
      relative_error[3] = 1.0 - Pri_re[3] / Pri[3];
      relative_error[4] = 1.0 - Pri_re[4] / Pri[4];
      relative_error[5] = 1.0 - Temp_re   / Temp;
 
      fprintf(fptr,"%+20.16e  ",relative_error[0]);
      fprintf(fptr,"%+20.16e  ",relative_error[1]);
      fprintf(fptr,"%+20.16e  ",relative_error[2]);
      fprintf(fptr,"%+20.16e  ",relative_error[3]);
      fprintf(fptr,"%+20.16e  ",relative_error[4]);
      fprintf(fptr,"%+20.16e  ",relative_error[5]);
      fprintf(fptr,"%+20.16e  ",dT    );
      fprintf(fptr,"%+20.16e  ",Pri[0]);
      fprintf(fptr,"%+20.16e  ",Pri[1]);
      fprintf(fptr,"%+20.16e  ",Pri[2]);
      fprintf(fptr,"%+20.16e  ",Pri[3]);
      fprintf(fptr,"%+20.16e  ",Pri[4]);
      fprintf(fptr,"%+20.16e  ",LorentzFactor);
      fprintf(fptr,"%+20.16e  ",Temp);
      fprintf(fptr,"\n");
  }



return 0;
}
