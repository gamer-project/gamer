#include "Plummer_calculator.h"

Plummer_calculator::Plummer_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );
}

Plummer_calculator::~Plummer_calculator()
{

}

void Plummer_calculator::init(double newton_g,double rho,double r){
  NEWTON_G=newton_g;
  Plummer_Rho0=rho;
  Plummer_R0=r;
  
}
double Plummer_calculator::set_vel(double r){
  const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
  const double Vmax_Fac    = sqrt( 2.0*NEWTON_G*TotM_Inf );

  double  Vmax, RanV, RanProb, Prob;

  Vmax = Vmax_Fac*pow( SQR(Plummer_R0) + SQR(r*Plummer_R0), -0.25 );

  //       randomly determine the velocity amplitude (ref: Aarseth, S. et al. 1974, A&A, 37, 183: Eq. [A4,A5])
  do
  {
    RanV    = RNG->GetValue( 0, 0.0, 1.0 );         // (0.0, 1.0)
    RanProb = RNG->GetValue( 0, 0.0, 0.1 );         // (0.0, 0.1)
    Prob    = SQR(RanV)*pow( 1.0-SQR(RanV), 3.5 );  // < 0.1
  }
  while ( RanProb > Prob );

  //       randomly set the velocity vector with the given amplitude (RanV*Vmax)
  return RanV*Vmax;

}  
