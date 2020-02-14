#include "Jaffe_calculator.h"
double Jaffe_NEWTON_G;
double Jaffe_Rho0;
double Jaffe_R0;
double Jaffe_MaxR;


//Quantities
double x_Jaffe[nbin_Jaffe];
double masses_Jaffe[nbin_Jaffe];
double pot_Jaffe[nbin_Jaffe];
double g_Jaffe[nbin_Jaffe];
Jaffe_calculator::Jaffe_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );

  prob_dens=new double*[size_Jaffe];
  int_prob_dens=new double*[size_Jaffe];
  psi=new double*[size_Jaffe];
  for(int i=0;i<size_Jaffe;i++){
    prob_dens[i]=new double();
    int_prob_dens[i]=new double();
    psi[i]=new double();
  }
}

Jaffe_calculator::~Jaffe_calculator()
{

}

//Statistics
double Jaffe_calculator::ave(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=a[k];
  }
  return sum/(fin-start);
}
double Jaffe_calculator::var_n(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=(a[k])*(a[k]);
  }
  sum=sum-(fin-start)*pow(ave(a,start,fin),2);
  return sum;
}
double Jaffe_calculator::cor(double* x,double* y,int start,int fin){
  double up=0,down = pow(var_n(x,start,fin)*var_n(y,start,fin),0.5);
  double ave_x = ave(x,start,fin),ave_y = ave(y,start,fin);
  for(int k=start;k<fin;k++){
    up+=(x[k]-ave_x)*(y[k]-ave_y);
  }
  return up/down;
}
void Jaffe_calculator::mask(double* x,int start,int fin){
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
void Jaffe_calculator::add_num(double* x,int start,int fin){
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
void Jaffe_calculator::smooth_all(double* x,int start,int fin){
  int num=10;
  for(int k=start;k<fin-num+1;k++){
    mask(x,k,k+num);
  }
  for(int k=start;k<fin-num+1;k++){
    add_num(x,k,k+num);
  }
}
double Jaffe_calculator::slope(double* x,double* y,int start,int fin){
  double cor_ = cor(x,y,start,fin);
  double var_n_x =var_n(x,start,fin), var_n_y =var_n(y,start,fin);
  double s =cor_*pow(var_n_y,0.5)/pow(var_n_x,0.5);

  return s;
}

//Physical Properties
double rho_Jaffe(double x){
  return Jaffe_Rho0 *pow(x*(1+x),-2)/(4*M_PI);
}
double mass_base_Jaffe(double x,void *nothing){
  return Jaffe_Rho0 *pow((1+x),-2)*pow(Jaffe_R0,3);
}
double mass_Jaffe(double x){
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);

  double  error;
  double result;

  gsl_function F;
  F.function = &mass_base_Jaffe;
  gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
  gsl_integration_workspace_free (w);
  return result;
}
/*double potential_Jaffe(double x){
  if(x>double(Jaffe_MaxR/Jaffe_R0)){
    return pot_Jaffe[nbin_Jaffe-1]*(Jaffe_MaxR)/(x*Jaffe_R0);
  }

  else{
    double dr = Jaffe_MaxR / (nbin_Jaffe-1);
    double r = x*Jaffe_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return pot_Jaffe[ind+1] + (pot_Jaffe[ind+1] - pot_Jaffe[ind])*par/dr;       
  }
}
double rho_dx_Jaffe(double x){
  return - Jaffe_Rho0 * (1/(x*x)) * (1/((1+x)*(1+x))) * (1/x + 1/(1+x))/(2*M_PI);
}
double de_rho_over_de_psi_Jaffe(double x) {
  double *s;
  double rho_dx_Jaffe = Jaffe_Rho0 * (1/(x*x)) * (1/((1+x)*(1+x))) * (1/x + 1/(1+x))/(2*M_PI);
  double psi_dx_Jaffe;
  if(x>double(Jaffe_MaxR/Jaffe_R0)){
    psi_dx_Jaffe=g_Jaffe[nbin_Jaffe-1]*pow(Jaffe_MaxR,2)/pow(x*Jaffe_R0,2);
    psi_dx_Jaffe*=Jaffe_R0; 
  }

  else{
    double dr = Jaffe_MaxR / (nbin_Jaffe-1);
    double r = x*Jaffe_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    psi_dx_Jaffe=g_Jaffe[ind+1] + (g_Jaffe[ind+1] - g_Jaffe[ind])*par/dr;   
    psi_dx_Jaffe*=Jaffe_R0;    
  }
    
  return rho_dx_Jaffe/psi_dx_Jaffe;
}*/
//NFW
/*
double potential_Jaffe(double x){
  return - (Jaffe_NEWTON_G*Jaffe_Rho0*Jaffe_R0*Jaffe_R0) * log(1+x) / x;
}
double rho_dx_Jaffe(double x){
  return - Jaffe_Rho0*(1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
}
double de_rho_over_de_psi_Jaffe(double x) {
  double rho_dx_Jaffe = - Jaffe_Rho0 *(1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
  double psi_dx_Jaffe = Jaffe_NEWTON_G*Jaffe_R0*Jaffe_R0*Jaffe_Rho0*( -log(1+x)/(x*x) + 1/(x*(1+x)) );
      
  return rho_dx_Jaffe/psi_dx_Jaffe;
}*/
//Jaffe

double potential_Jaffe(double x){
   double s=(Jaffe_NEWTON_G*Jaffe_Rho0*pow(Jaffe_R0,2))*log(x/(1+x));
   return s;
}
double rho_dx_Jaffe(double x){
  return -2 *Jaffe_Rho0 *(1/x + 1/(1+x))/pow(x*(1+x),2);
}
double de_rho_over_de_psi_Jaffe(double x) {
  double rho_dx_Jaffe = Jaffe_Rho0 * (1/(x*x)) * (1/((1+x)*(1+x))) * (1/x + 1/(1+x))/(2*M_PI);
  double psi_dx_Jaffe = Jaffe_NEWTON_G*Jaffe_R0*Jaffe_R0*Jaffe_Rho0*( 1/(x*(1+x)) );
      
  return rho_dx_Jaffe/psi_dx_Jaffe;
}
//Analytical Prob_Dens
double Jaffe_d_minus(double x){
  return pow(M_PI,0.5)*exp(x*x)*gsl_sf_erf(x)/2;
}
double Jaffe_distri_f(double e){
  double result=Jaffe_d_minus(pow(2*e,0.5))-pow(2,0.5)*Jaffe_d_minus(pow(e,0.5))-pow(2,0.5)*gsl_sf_dawson(pow(e,0.5))+gsl_sf_dawson(pow(2*e,0.5));
  double m = Jaffe_Rho0 * pow(Jaffe_R0,3);
  result *= pow(2,0.5)/(M_PI*pow(Jaffe_NEWTON_G*Jaffe_R0,1.5)*pow(m,0.5));
  return result;
}

//GSL functions
double psi_potential_Jaffe(double r, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_Jaffe(r);
}
double inverse_psi_to_x_Jaffe (double psi) {
  double psi1=psi;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r0= 0;
  double x_lo =0.0001, x_hi = 100000000000.0;
  gsl_function F;
  
  F.function = &psi_potential_Jaffe;
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
double integration_Jaffe(double eng){
  //double epson=0.1;
  double min = 0;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    double x0 = inverse_psi_to_x_Jaffe(min+(i+0.5)*dx);
    if(i==num-1)result += 2* de_rho_over_de_psi_Jaffe(x0) * ( pow(eng-psi_l,0.5) );
    else result += 2* de_rho_over_de_psi_Jaffe(x0) * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
    
  }
  
  return result;
}
double integration_eng_base_Jaffe(double eng){
  //double x0 =inverse_psi_to_x_Jaffe(eng);
  return integration_Jaffe(eng);
}
void Jaffe_calculator::initialize_mass(){
  double dx = Jaffe_MaxR/(Jaffe_R0*(nbin_Jaffe-1));
  for(int i=0;i<nbin_Jaffe;i++){
    x_Jaffe[i] = dx *i;
    masses_Jaffe[i] = mass_Jaffe(x_Jaffe[i]);
    g_Jaffe[i] = masses_Jaffe[i]*Jaffe_NEWTON_G/pow(x_Jaffe[i]*Jaffe_R0,2);
    
  }
  g_Jaffe[0]=0;
}
void Jaffe_calculator::initialize_pot(){
  double dx = Jaffe_MaxR/(Jaffe_R0*(nbin_Jaffe-1));
  double pot_out = - masses_Jaffe[nbin_Jaffe-1]*Jaffe_NEWTON_G/Jaffe_MaxR;
  for(int i=nbin_Jaffe-2;i>=0;i--){
    pot_Jaffe[i]=pot_Jaffe[i+1] - (g_Jaffe[i]+g_Jaffe[i+1])*Jaffe_R0*dx/2;
   
  }
}
void Jaffe_calculator::initialize_prob_dens(){
  double min =0,max =-potential_Jaffe(1e-2);/***difference***/
  delta =(max-min)/size_Jaffe;
  
  double eng=min;

  
  for(int k =0;k<size_Jaffe;k++){
    *psi[k] = eng;
    //*int_prob_dens[k] = integration_eng_base_Jaffe(eng);
    
    
    eng +=delta;
  
  }
  for(int k =0;k<size_Jaffe;k++){

    /*if(k==0)*prob_dens[k]=slope(psi,int_prob_dens,k,k+5);
    else if(k==1)*prob_dens[k]=slope(psi,int_prob_dens,k-1,k+4);

    else if(k==size_Jaffe-2)*prob_dens[k]=slope(psi,int_prob_dens,k-3,k+2);
    else if(k==size_Jaffe-1)*prob_dens[k]=slope(psi,int_prob_dens,k-4,k+1);

    else *prob_dens[k]=slope(psi,int_prob_dens,k-2,k+3);

    if(*prob_dens[k]<0)*prob_dens[k]=0;*/
    double fac = Jaffe_NEWTON_G * Jaffe_Rho0 * pow(Jaffe_R0,2);
  
    *prob_dens[k] = Jaffe_distri_f(*psi[k]/fac);
    //cout<<*prob_dens[k]<<endl;
  }
  //smooth_all(prob_dens,0,size_Jaffe);
}

//Principal functions
void Jaffe_calculator::init(double newton_g,double rho,double r0,double maxr){
  
  Jaffe_NEWTON_G=newton_g;
  Jaffe_Rho0=rho;
  Jaffe_R0=r0;
  Jaffe_MaxR=maxr;
  

  initialize_mass();
  initialize_pot();
  initialize_prob_dens();
  
}
double Jaffe_calculator::set_vel(double r){  
  int index=0;
  double sum=0;
  double psi_per =-potential_Jaffe(r);
  for(int k =0;k<size_Jaffe;k++){
    if(*psi[k]>psi_per){
      index =k-1;
      break;
    }
    else if(k==size_Jaffe-1)index=size_Jaffe-1;
    sum += *prob_dens[k] *pow(psi_per-*psi[k],0.5) *delta;
  }
  double sum_rad,sum_mes=0,par,psi_ass;
  int index_ass=0;
  
  sum_rad = RNG->GetValue( 0, 0.0, 1.0 ); 
  sum_rad*=sum;

  for(int k =0;k<size_Jaffe;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(*prob_dens[index_ass] *pow(psi_per-*psi[index_ass],0.5) *delta);
      break;
      }
      else if(k==size_Jaffe-1)index_ass=size_Jaffe-1;
    sum_mes += *prob_dens[k] *pow(psi_per-*psi[k],0.5) *delta;
  }
  //cout<<"x:"<<r<<'\t'<<"index:"<<index<<'\t'<<"index_ass:"<<index_ass<<endl;
  if(index<index_ass)cout<<"x:"<<r<<'\t'<<"index:"<<index<<'\t'<<"index_ass:"<<index_ass<<endl;
  psi_ass = *psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_Jaffe(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_Jaffe(r)),0.5);
  return v;
}  
double Jaffe_calculator::set_mass(double x){
  return mass_Jaffe(x);
}