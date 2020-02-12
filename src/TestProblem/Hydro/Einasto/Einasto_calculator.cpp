#include "Einasto_calculator.h"
double Einasto_NEWTON_G;
double Einasto_Rho0;
double Einasto_R0;
double Einasto_MaxR;
double alpha;

//Quantities
double x_Einasto[nbin_Einasto];
double masses_Einasto[nbin_Einasto];
double pot_Einasto[nbin_Einasto];
double g_Einasto[nbin_Einasto];
Einasto_calculator::Einasto_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );
}

Einasto_calculator::~Einasto_calculator()
{

}

//Statistics
double Einasto_calculator::ave(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=a[k];
  }
  return sum/(fin-start);
}
double Einasto_calculator::var_n(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=(a[k])*(a[k]);
  }
  sum=sum-(fin-start)*pow(ave(a,start,fin),2);
  return sum;
}
double Einasto_calculator::cor(double* x,double* y,int start,int fin){
  double up=0,down = pow(var_n(x,start,fin)*var_n(y,start,fin),0.5);
  double ave_x = ave(x,start,fin),ave_y = ave(y,start,fin);
  for(int k=start;k<fin;k++){
    up+=(x[k]-ave_x)*(y[k]-ave_y);
  }
  return up/down;
}
void Einasto_calculator::mask(double* x,int start,int fin){
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
void Einasto_calculator::add_num(double* x,int start,int fin){
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
void Einasto_calculator::smooth_all(double* x,int start,int fin){
  int num=10;
  for(int k=start;k<fin-num+1;k++){
    mask(x,k,k+num);
  }
  for(int k=start;k<fin-num+1;k++){
    add_num(x,k,k+num);
  }
}
double Einasto_calculator::slope(double* x,double* y,int start,int fin){
  double cor_ = cor(x,y,start,fin);
  double var_n_x =var_n(x,start,fin), var_n_y =var_n(y,start,fin);
  double s =cor_*pow(var_n_y,0.5)/pow(var_n_x,0.5);

  return s;
}

//Physical Properties
double rho_Einasto(double x){
  return Einasto_Rho0 *exp(-pow(x,alpha));
}
double mass_base_Einasto(double x,void *nothing){
  return 4*M_PI*Einasto_Rho0*pow(Einasto_R0,3)*pow(x,2) *exp(-pow(x,alpha));
}
double mass_Einasto(double x){
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);

  double  error;
  double result;

  gsl_function F;
  F.function = &mass_base_Einasto;
  gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
  gsl_integration_workspace_free (w);
  return result;
}
double potential_Einasto(double x){
  if(x>double(Einasto_MaxR/Einasto_R0)){
    return pot_Einasto[nbin_Einasto-1]*(Einasto_MaxR)/(x*Einasto_R0);
  }

  else{
    double dr = Einasto_MaxR / (nbin_Einasto-1);
    double r = x*Einasto_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return pot_Einasto[ind+1] + (pot_Einasto[ind+1] - pot_Einasto[ind])*par/dr;       
  }
}
double rho_dx_Einasto(double x){
  return - alpha*pow(x,alpha-1)*Einasto_Rho0 *exp(-pow(x,alpha));
}
double de_rho_over_de_psi_Einasto(double x) {
  double *s;
  double rho_dx_Einasto =  alpha*pow(x,alpha-1)*Einasto_Rho0 *exp(-pow(x,alpha));
  double psi_dx_Einasto;
  if(x>double(Einasto_MaxR/Einasto_R0)){
    psi_dx_Einasto=g_Einasto[nbin_Einasto-1]*pow(Einasto_MaxR,2)/pow(x*Einasto_R0,2);
    psi_dx_Einasto*=Einasto_R0; 
  }

  else{
    double dr = Einasto_MaxR / (nbin_Einasto-1);
    double r = x*Einasto_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    psi_dx_Einasto=g_Einasto[ind+1] + (g_Einasto[ind+1] - g_Einasto[ind])*par/dr;   
    psi_dx_Einasto*=Einasto_R0;    
  }
    
  return rho_dx_Einasto/psi_dx_Einasto;
}

//GSL functions
double psi_potential_Einasto(double r, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_Einasto(r);
}
double inverse_psi_to_x_Einasto (double psi) {
  double psi1=psi;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r0= 0;
  double x_lo =0.0001, x_hi = 100000000000.0;
  gsl_function F;
  
  F.function = &psi_potential_Einasto;
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
double integration_Einasto(double eng){
  //double epson=0.1;
  double min = 0;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    double x0 = inverse_psi_to_x_Einasto(min+(i+0.5)*dx);
    if(i==num-1)result += 2* de_rho_over_de_psi_Einasto(x0) * ( pow(eng-psi_l,0.5) );
    else result += 2* de_rho_over_de_psi_Einasto(x0) * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
    
  }
  cout<<result<<endl;
  return result;
}
double integration_eng_base_Einasto(double eng){
  //double x0 =inverse_psi_to_x_Einasto(eng);
  return integration_Einasto(eng);
}
void Einasto_calculator::initialize_mass(){
  double dx = Einasto_MaxR/(Einasto_R0*(nbin_Einasto-1));
  for(int i=0;i<nbin_Einasto;i++){
    x_Einasto[i] = dx *i;
    masses_Einasto[i] = mass_Einasto(x_Einasto[i]);
    g_Einasto[i] = masses_Einasto[i]*Einasto_NEWTON_G/(x_Einasto[i]*Einasto_R0);
    
  }
  g_Einasto[0]=0;
}
void Einasto_calculator::initialize_pot(){
  double dx = Einasto_MaxR/(Einasto_R0*(nbin_Einasto-1));
  double pot_out = - masses_Einasto[nbin_Einasto-1]*Einasto_NEWTON_G/Einasto_MaxR;
  for(int i=nbin_Einasto-2;i>=0;i--){
    pot_Einasto[i]=pot_Einasto[i+1] - (g_Einasto[i]+g_Einasto[i+1])*Einasto_R0*dx/2;
   
  }
}
void Einasto_calculator::initialize_prob_dens(){
  double min =-potential_Einasto(100),max =-potential_Einasto(0.001);/***difference***/
  delta =(max-min)/size_Einasto;
  double eng=min;

  
  for(int k =0;k<size_Einasto;k++){
    psi[k] = eng;
    int_prob_dens[k] = integration_eng_base_Einasto(eng);
    
    
    eng +=delta;
  
  }
  for(int k =0;k<size_Einasto;k++){

    if(k==0)prob_dens[k]=slope(psi,int_prob_dens,k,k+5);
    else if(k==1)prob_dens[k]=slope(psi,int_prob_dens,k-1,k+4);

    else if(k==size_Einasto-2)prob_dens[k]=slope(psi,int_prob_dens,k-3,k+2);
    else if(k==size_Einasto-1)prob_dens[k]=slope(psi,int_prob_dens,k-4,k+1);

    else prob_dens[k]=slope(psi,int_prob_dens,k-2,k+3);

    if(prob_dens[k]<0)prob_dens[k]=0;
    
  }
  smooth_all(prob_dens,0,size_Einasto);
}

//Principal functions
void Einasto_calculator::init(double al,double newton_g,double rho,double r0,double maxr){

  Einasto_NEWTON_G=newton_g;
  Einasto_Rho0=rho;
  Einasto_R0=r0;
  Einasto_MaxR=maxr;
  alpha=al;

  initialize_mass();
  initialize_pot();
  initialize_prob_dens();
  
}
double Einasto_calculator::set_vel(double r){  
  double index,sum=0;
  double psi_per =-potential_Einasto(r);
  for(int k =0;k<size_Einasto;k++){
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

  for(int k =0;k<size_Einasto;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(prob_dens[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
  }
  psi_ass = psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_Einasto(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_Einasto(r)),0.5);
  return v;
}  
double Einasto_calculator::set_mass(double x){
  return mass_Einasto(x);
}