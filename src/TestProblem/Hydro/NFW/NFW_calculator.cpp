#include "NFW_calculator.h"
double NFW_NEWTON_G;
double NFW_Rho0;
double NFW_R0;
double NFW_MaxR;
double NFW_mass;
NFW_calculator::NFW_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );
}

NFW_calculator::~NFW_calculator()
{

}

//Statistics
double NFW_calculator::ave(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=a[k];
  }
  return sum/(fin-start);
}
double NFW_calculator::var_n(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=(a[k])*(a[k]);
  }
  sum=sum-(fin-start)*pow(ave(a,start,fin),2);
  return sum;
}
double NFW_calculator::cor(double* x,double* y,int start,int fin){
  double up=0,down = pow(var_n(x,start,fin)*var_n(y,start,fin),0.5);
  double ave_x = ave(x,start,fin),ave_y = ave(y,start,fin);
  for(int k=start;k<fin;k++){
    up+=(x[k]-ave_x)*(y[k]-ave_y);
  }
  return up/down;
}
void NFW_calculator::mask(double* x,int start,int fin){
  double standard=3;
  for(int j=start;j<fin;j++){
    bool flag=0;
    for(int k=start;k<fin;k++){
      double test_fac;
      if(x[k]!=0){
        test_fac=fabs(x[j]/x[k]);
        if(test_fac>standard)flag=1;
      }
      if(flag)x[j]=0;
    }
  }
}
void NFW_calculator::add_num(double* x,int start,int fin){
  double sum=0;
  int num=0;
  for(int j=start;j<fin;j++){
    if(x[j]!=0)sum+=x[j];
    num++;
  }
  double ave_x;
  if(num!=0){
    ave_x=sum/(num+0.0);
    for(int j=start;j<fin;j++){
      if(x[j]==0)x[j]=ave_x;
    }
  }
}
void NFW_calculator::smooth_all(double* x,int start,int fin){
  int num=10;
  for(int k=start;k<fin-num+1;k++){
    mask(x,k,k+num);
  }
  for(int k=start;k<fin-num+1;k++){
    add_num(x,k,k+num);
  }
}
double NFW_calculator::slope(double* x,double* y,int start,int fin){
  double cor_ = cor(x,y,start,fin);
  double var_n_x =var_n(x,start,fin), var_n_y =var_n(y,start,fin);
  double s =cor_*pow(var_n_y,0.5)/pow(var_n_x,0.5);

  return s;
}

//Physical Properties
double potential_NFW(double x){
  /*double M = NFW_mass;//4*M_PI*NFW_Rho0*pow(NFW_R0,3)*(log(1+x)-x/(1+x));
  double x0=NFW_MaxR/NFW_R0;
  if(x>x0) return -NFW_NEWTON_G*M/(NFW_R0*x);
  return - (4*M_PI*NFW_NEWTON_G*NFW_Rho0*NFW_R0*NFW_R0) * (log(1+x) / x- log(1+x0) / x0)-NFW_NEWTON_G*M/(NFW_R0*x0);*/
  return- (4*M_PI*NFW_NEWTON_G*NFW_Rho0*NFW_R0*NFW_R0) * (log(1+x) / x);
}
double mass_base_NFW(double x,void* nothing){
  double x0 = 0.7*NFW_MaxR/NFW_R0;
  double xmax = NFW_MaxR/NFW_R0;
  
  if(x<x0) return 4*M_PI*pow(NFW_R0,3)*(NFW_Rho0*(x/((1+x)*(1+x))));
  else {
    double rho0 = NFW_Rho0*(1/(x0*(1+x0)*(1+x0)));
    double rho = rho0 *pow(x/x0,-10); //( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(NFW_R0*x,2)* NFW_R0 *rho;
  }
  //return 4*M_PI*pow(NFW_R0,3)*(NFW_Rho0*(x/((1+x)*(1+x))));
}
double NFW_calculator::set_mass(double x){
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);

  double  error;
  double result;

  gsl_function F;
  F.function = &mass_base_NFW;
  gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
  gsl_integration_workspace_free (w);
  return result;
}
double rho_dx_NFW(double x){
  return - NFW_Rho0*(1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
}
double de_rho_over_de_psi_NFW(double x) {
  double rho_dx_NFW = - NFW_Rho0 *(1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
  double psi_dx_NFW = 4*M_PI*NFW_NEWTON_G*NFW_R0*NFW_R0*NFW_Rho0*( -log(1+x)/(x*x) + 1/(x*(1+x)) );
      
  return rho_dx_NFW/psi_dx_NFW;
}

//GSL functions
double psi_potential_NFW(double r, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_NFW(r);
}
double inverse_psi_to_x_NFW (double psi) {
  double psi1=psi;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r0= 0;
  double x_lo =0.0001, x_hi = 100000000000.0;
  gsl_function F;
  
  F.function = &psi_potential_NFW;
  F.params=&psi1;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r0= gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  return r0;
}

//Probability Density
double integration_NFW(double eng){
  double min = 0;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    double x0 = inverse_psi_to_x_NFW(min+(i+0.5)*dx);
    if(i==num-1)result += 2* de_rho_over_de_psi_NFW(x0) * ( pow(eng-psi_l,0.5) );
    else result += 2* de_rho_over_de_psi_NFW(x0) * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
    
  }
  cout<<result<<endl;
  return result;
}
double integration_eng_base_NFW(double eng){
  //double x0 =inverse_psi_to_x_NFW(eng);
  return integration_NFW(eng);
}
void NFW_calculator::initialize_prob_dens(){
  double min =-potential_NFW(100),max =-potential_NFW(0.001);/***difference***/
  delta =(max-min)/size_NFW;
  double eng=min;

  
  for(int k =0;k<size_NFW;k++){
    psi[k] = eng;
    int_prob_dens[k] = integration_eng_base_NFW(eng);
    
    
    eng +=delta;
  
  }
  for(int k =0;k<size_NFW;k++){

    if(k==0)prob_dens[k]=slope(psi,int_prob_dens,k,k+5);
    else if(k==1)prob_dens[k]=slope(psi,int_prob_dens,k-1,k+4);

    else if(k==size_NFW-2)prob_dens[k]=slope(psi,int_prob_dens,k-3,k+2);
    else if(k==size_NFW-1)prob_dens[k]=slope(psi,int_prob_dens,k-4,k+1);

    else prob_dens[k]=slope(psi,int_prob_dens,k-2,k+3);

    if(prob_dens[k]<0)prob_dens[k]=0;
    
  }
  smooth_all(prob_dens,0,size_NFW);
}

//Principal functions
void NFW_calculator::init(double newton_g,double rho,double r0,double maxr){

  NFW_NEWTON_G=newton_g;
  NFW_Rho0=rho;
  NFW_R0=r0;
  NFW_MaxR=maxr;
  NFW_mass=set_mass(NFW_MaxR/NFW_R0);

  
  initialize_prob_dens();
  
}
double NFW_calculator::set_vel(double r){  
  double index,sum=0;
  double psi_per =-potential_NFW(r);//+potential_NFW(NFW_MaxR/NFW_R0);
  for(int k =0;k<size_NFW;k++){
    if(psi[k]>psi_per){
      index =k-1;
      break;
    }
    sum += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
  }
  double sum_rad,sum_mes=0,par,psi_ass;
  int index_ass;

  sum_rad = RNG->GetValue( 0, 0.0, 1.0 ); 
  sum_rad*=sum;

  for(int k =0;k<size_NFW;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(prob_dens[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
  }
  psi_ass = psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_NFW(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_NFW(r)),0.5);

  //Truncation
  /*double eng = pow(v,2)/2 + potential_NFW(r);
  double eng_trunc = potential_NFW(NFW_MaxR/NFW_R0);
  if(eng>eng_trunc)cout<<"Trunc!"<<endl;*/

  return v;
}  