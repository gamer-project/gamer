#include "NFW_calculator.h"
<<<<<<< HEAD
double NFW_NEWTON_G;
double NFW_Rho0;
double NFW_R0;
NFW_calculator::NFW_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );
=======

NFW_calculator::NFW_calculator()
{

>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
}

NFW_calculator::~NFW_calculator()
{

}
<<<<<<< HEAD

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
  return - (4*M_PI*NFW_NEWTON_G*NFW_Rho0*NFW_R0*NFW_R0) * log(1+x) / x;
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
=======
double potential_NFW(double x){
  return - (4*M_PI*NEWTON_G*NFW_Rho0*NFW_R0*NFW_R0) * log(1+x) / x;
}
double rho_dx_NFW(double x){
  return -NFW_Rho0 * (1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
}
double de_rho_over_de_psi_NFW(double x) {
  double rho_dx_NFW = -NFW_Rho0 * (1/(x*x*(1+x)*(1+x)) + 2/(x*(1+x)*(1+x)*(1+x)));
  double psi_dx_NFW = 4*M_PI*NEWTON_G*NFW_R0*NFW_R0*NFW_Rho0*( -log(1+x)/(x*x) + 1/(x*(1+x)) );
      
  return rho_dx_NFW/psi_dx_NFW;
}
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
double psi_potential_NFW(double r, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_NFW(r);
}
<<<<<<< HEAD
double inverse_psi_to_x_NFW (double psi) {
=======
double inverse_psi_to_r_NFW (double psi) {
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
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
<<<<<<< HEAD

//Probability Density
double integration_NFW(double eng){
  //double epson=0.1;
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
void NFW_calculator::init(double newton_g,double rho,double r){

  NFW_NEWTON_G=newton_g;
  NFW_Rho0=rho;
  NFW_R0=r;

  
  initialize_prob_dens();
  
}
double NFW_calculator::set_vel(double r){  
  double index,sum=0;
  double psi_per =-potential_NFW(r);
  for(int k =0;k<size_NFW;k++){
=======
double basis_NFW (double r, void * energy_negative) {
  double e= *(double *)energy_negative;
  return abs(rho_dx_NFW(r))*pow(e+potential_NFW(r),-0.5);
}
double integration_NFW(double r0, void * energy_negative){
  double epson=0.1;
  double max=10000;
  /***part I***/
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);
  double e= *(double *)energy_negative;

  double  error;
  double result;
  gsl_function F;
  F.function = &basis_NFW;
  F.params = &e;
  
  if(e+potential_NFW(r0+epson)<0) return 0;
  gsl_integration_qags (&F, r0+epson, max, 0, 1e-7, 1000,
                        w, &result, &error); 
  gsl_integration_workspace_free (w);


  /***part II***/
  double mean = (de_rho_over_de_psi_NFW(r0) + de_rho_over_de_psi_NFW(r0+epson))/2;
  double much=2*mean*pow(e+potential_NFW(r0+epson),0.5);
  
  result=result+much;
  return result;
}
double integration_eng_base_NFW(double eng,void *s){
  double r0 =inverse_psi_to_r_NFW(eng);
  return integration_NFW(r0,&eng);
}
double prob_dens_NFW(double eng){
  gsl_function F;
  double result, abserr;
  F.function = &integration_eng_base_NFW;
  gsl_deriv_central (&F, eng, 1e-8, &result, &abserr);

  return result;
}
void NFW_calculator::init(){
  string line;

  double min =-potential_NFW(20),max =-potential_NFW(0.1);
  delta =(max-min)/size;
  double eng=min;
  for(int k =0;k<size;k++){

    psi[k] = eng;
    
    f[k] = prob_dens_NFW(eng);
    if(f[k]<0)f[k]=0;
    
    int sec_num=8;
    int s=k*sec_num/size;

    double d=0.2;
    eng +=delta*(1+d*(s+0.5-sec_num/2));
  
  }
  
}
double NFW_calculator::set_vel(double r){
    
  double index,sum=0;
  double psi_per =-potential_NFW(r);
  for(int k =0;k<size;k++){
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
    if(psi[k]>psi_per){
      index =k-1;
      break;
    }
<<<<<<< HEAD
    sum += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
=======
    sum += f[k] *pow(psi_per-psi[k],0.5) *delta;
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
  }
  double sum_rad,sum_mes=0,par,psi_ass;
  int index_ass;

<<<<<<< HEAD
  sum_rad = RNG->GetValue( 0, 0.0, 1.0 ); 
  sum_rad*=sum;

  for(int k =0;k<size_NFW;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(prob_dens[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
=======
  random_device rd;  
  mt19937 gen(rd()); 
  uniform_real_distribution<> pro_dis(0.0, sum);
  sum_rad = pro_dis(gen);

  for(int k =0;k<size;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(f[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += f[k] *pow(psi_per-psi[k],0.5) *delta;
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
  }
  psi_ass = psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_NFW(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_NFW(r)),0.5);
  return v;
}  