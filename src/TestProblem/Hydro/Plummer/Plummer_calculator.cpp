#include "Plummer_calculator.h"

Plummer_calculator::Plummer_calculator()
{
<<<<<<< HEAD
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );
=======
  
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
}

Plummer_calculator::~Plummer_calculator()
{

}
<<<<<<< HEAD

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
=======
//#define DEBUG
/***Gas***/
//x=r/Plummer_R0

double Plummer_rho_gas(double x){
  return Plummer_Rho0*Plummer_GasMFrac/pow((1+x*x)*x,2.5);
}
double Plummer_potential_gas(double x){
   double s=-(4*M_PI*NEWTON_G*Plummer_Rho0*Plummer_GasMFrac*pow(Plummer_R0,2)/3)*pow(1+x*x,-0.5);
   return s;
}
double Plummer_g_gas(double x){
  double s=(4*M_PI*NEWTON_G*Plummer_Rho0*Plummer_GasMFrac/3)*pow(1+x*x,-1.5);
  return s;
}
double Plummer_calculator::pressure(double x){
  //cout<<"Plummer_Rho0:"<<Plummer_Rho0<<",  Plummer_Rho0*Plummer_GasMFrac:"<<Plummer_Rho0*Plummer_GasMFrac<<endl;
  return (2*M_PI*NEWTON_G*Plummer_Rho0*Plummer_Rho0*Plummer_GasMFrac*pow(Plummer_R0,2))*pow(1+x*x,-3)/9;
}

/***Particle***/
double Plummer_potential(double x){
   double s=-(4*M_PI*NEWTON_G*Plummer_Rho0*pow(Plummer_R0,2)/3)*pow(1+x*x,-0.5);
   return s;
}
double Plummer_psi(double x){
  return -Plummer_potential(x);
}
double Plummer_epsilon(double v,void *r){
  double r0=*(double *)r;
  return Plummer_psi(r0)-0.5*v*v;
}
double Plummer_calculator::psi(double x){
  return -Plummer_potential(x);
}
double Plummer_distri_f(double e){
  double result=pow(e,3.5);
  return result;
}
double Plummer_f_v_base(double v,void *r){
  double e=Plummer_epsilon(v,r);
  return Plummer_distri_f(e);
}

double Plummer_distribution(double v,void *r){
  return v*v*Plummer_f_v_base(v,r);
}
double Plummer_deriv_dis(double v,void *r){/****/
  gsl_function F;
      double result, abserr;
      
      F.function = &Plummer_distribution;
      F.params=r;
      gsl_deriv_central (&F, v, 1e-8, &result, &abserr);
      return result;     
}
double Plummer_best_v(double r){
  double delta=0.0001;

  double r0=r;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double v0= 0;
  double x_lo =delta, x_hi = pow(2*Plummer_psi(r),0.5)-delta;
  gsl_function F;
  
  F.function = &Plummer_deriv_dis;
  F.params=&r0;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do
  {
  iter++;
  status = gsl_root_fsolver_iterate (s);
  v0= gsl_root_fsolver_root (s);
  x_lo = gsl_root_fsolver_x_lower (s);
  x_hi = gsl_root_fsolver_x_upper (s);
  status = gsl_root_test_interval (x_lo, x_hi,
                                                      0, 0.001);

  }
  while (status == GSL_CONTINUE && iter < max_iter);        
  gsl_root_fsolver_free (s);
  return v0;
}

double Plummer_calculator::prob(double v,void* r){
  return Plummer_distribution(v,r);
}
double Plummer_calculator::max_prob(double r){
  double *r0=new double(r);
  return Plummer_distribution(Plummer_best_v(r),r0);
}
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
