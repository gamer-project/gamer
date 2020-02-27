#include "UNKNOWN_calculator.h"
double UNKNOWN_NEWTON_G;
double UNKNOWN_Rho0;
double UNKNOWN_R0;
double UNKNOWN_MaxR;

double *Table_MassProf_r_UNKNOWN;
double *Table_MassProf_M_UNKNOWN;
double *Table_MassProf_rho_UNKNOWN;
double *Table_MassProf_rhodx_UNKNOWN;
double *Table_MassProf_derho_overdx_UNKNOWN;
double *Table_MassProf_g_UNKNOWN;
double *Table_MassProf_pot_UNKNOWN;

int    UNKNOWN_MassProfNBin;
UNKNOWN_calculator::UNKNOWN_calculator()
{
  RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, 123 );

}

UNKNOWN_calculator::~UNKNOWN_calculator()
{

}
//statistics
double UNKNOWN_calculator::ave(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=a[k];
  }
  return sum/(fin-start);
}
double UNKNOWN_calculator::var_n(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=(a[k])*(a[k]);
  }
  sum=sum-(fin-start)*pow(ave(a,start,fin),2);
  return sum;
}
double UNKNOWN_calculator::cor(double* x,double* y,int start,int fin){
  double up=0,down = pow(var_n(x,start,fin)*var_n(y,start,fin),0.5);
  double ave_x = ave(x,start,fin),ave_y = ave(y,start,fin);
  for(int k=start;k<fin;k++){
    up+=(x[k]-ave_x)*(y[k]-ave_y);
  }
  return up/down;
}
void UNKNOWN_calculator::mask(double* x,int start,int fin){
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
void UNKNOWN_calculator::add_num(double* x,int start,int fin){
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
void UNKNOWN_calculator::smooth_all(double* x,int start,int fin){
  int num=10;
  for(int k=start;k<fin-num+1;k++){
    mask(x,k,k+num);
  }
  for(int k=start;k<fin-num+1;k++){
    add_num(x,k,k+num);
  }
}

double UNKNOWN_calculator::slope(double* x,double* y,int start,int fin){
  double cor_ = cor(x,y,start,fin);
  double var_n_x =var_n(x,start,fin), var_n_y =var_n(y,start,fin);
  double s =cor_*pow(var_n_y,0.5)/pow(var_n_x,0.5);

  return s;
}

double potential_UNKNOWN(double x){

  if(x>double(UNKNOWN_MaxR/UNKNOWN_R0)){
    return Table_MassProf_pot_UNKNOWN[UNKNOWN_MassProfNBin-1]*(UNKNOWN_MaxR)/(x*UNKNOWN_R0);
  }

  else{
    double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);
    double r = x*UNKNOWN_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return Table_MassProf_pot_UNKNOWN[ind+1] + (Table_MassProf_pot_UNKNOWN[ind+1] - Table_MassProf_pot_UNKNOWN[ind])*par/dr;       
  }
  
}
double rho_UNKNOWN(double x){
  
  if(x>UNKNOWN_MaxR/UNKNOWN_R0){
    return 0;
  }

  else{
    double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);
    double r = x*UNKNOWN_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return Table_MassProf_rho_UNKNOWN[ind+1] + (Table_MassProf_rho_UNKNOWN[ind+1] - Table_MassProf_rho_UNKNOWN[ind])*par/dr;       
  }
}

double rho_dx_UNKNOWN(double x){
  double rd=10;
  if(x>(UNKNOWN_MaxR/UNKNOWN_R0+rd)){
    return Table_MassProf_rhodx_UNKNOWN[UNKNOWN_MassProfNBin-1]*exp(-x+UNKNOWN_MaxR/UNKNOWN_R0+rd);
  }
  else if(x>UNKNOWN_MaxR/UNKNOWN_R0 and x<(UNKNOWN_MaxR/UNKNOWN_R0+rd)){
    
    return  Table_MassProf_rhodx_UNKNOWN[UNKNOWN_MassProfNBin-1];
  }

  else{
    double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);
    double r = x*UNKNOWN_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return Table_MassProf_rhodx_UNKNOWN[ind+1] + (Table_MassProf_rhodx_UNKNOWN[ind+1] - Table_MassProf_rhodx_UNKNOWN[ind])*par/dr;       
  }
}
double de_rho_over_de_psi_UNKNOWN(double x){
  
  if(x>UNKNOWN_MaxR/UNKNOWN_R0){
      return Table_MassProf_derho_overdx_UNKNOWN[UNKNOWN_MassProfNBin-1]*pow(UNKNOWN_MaxR,2)/pow(UNKNOWN_R0*x,2);
  }

  else{
    double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);
    double r = x*UNKNOWN_R0;
    int ind = r/dr;
    double par = r - ind * dr;
    
    return Table_MassProf_derho_overdx_UNKNOWN[ind+1] + (Table_MassProf_derho_overdx_UNKNOWN[ind+1] - Table_MassProf_derho_overdx_UNKNOWN[ind])*par/dr;         
  }
}


double psi_potential_UNKNOWN(double x, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_UNKNOWN(x);
}
double inverse_psi_to_x_UNKNOWN (double psi) {
  double psi1=psi;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double x0= 0;
  double x_lo =0.0001, x_hi = 100000000000.0;
  gsl_function F;
  
  F.function = &psi_potential_UNKNOWN;
  F.params=&psi1;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    x0= gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,0, 0.001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  return x0;
}


double basis_UNKNOWN (double x, void * energy_negative) {
  double e= *(double *)energy_negative;
  return fabs(rho_dx_UNKNOWN(x))*pow(e+potential_UNKNOWN(x),-0.5);
}

double integration_UNKNOWN(double x0, void * energy_negative){
  double epson=0.1;
  double max=10000;
  /***part I***/
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);
  double e= *(double *)energy_negative;

  double  error;
  double result_all,result_high;
  gsl_function F;
  F.function = &basis_UNKNOWN;
  F.params = &e;
  
  if(e+potential_UNKNOWN(x0+epson)<0) return 0;
  gsl_integration_qags (&F, x0+epson, max, 0, 1e-7, 1000,
                        w, &result_all, &error); 
  gsl_integration_qags (&F, UNKNOWN_MaxR/UNKNOWN_R0, max, 0, 1e-7, 1000,
                        w, &result_high, &error); 
  gsl_integration_workspace_free (w);


  /***part II***/
  double mean = (de_rho_over_de_psi_UNKNOWN(x0) + de_rho_over_de_psi_UNKNOWN(x0+epson))/2;
  double much=2*mean*pow(e+potential_UNKNOWN(x0+epson),0.5);
  
  
  return result_all-result_high+much;
}
double integration_simpson(double eng){
  //double epson=0.1;
  double min = 0;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  double result_right=0,result_left=0,result_simpson=0;
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    double x0 = inverse_psi_to_x_UNKNOWN(min+(i+0.5)*dx);
    if(i==num-1)result += 2* de_rho_over_de_psi_UNKNOWN(x0) * ( pow(eng-psi_l,0.5) );
    else result += 2* de_rho_over_de_psi_UNKNOWN(x0) * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
    
  }
  return result;
}

double integration_eng_base_UNKNOWN(double eng){
  //double x0 =inverse_psi_to_x_UNKNOWN(eng);
  return integration_simpson(eng);
}

/*double UNKNOWN_calculator::MassProf_UNKNOWN( const double r )
{ 
  //main
  const double x = r / UNKNOWN_R0;
  return 4.0/3.0*M_PI*UNKNOWN_Rho0*CUBE(r)*pow( 1.0+x*x, -1.5 );
}*/
double mass_base_UNKNOWN(double x,void* nothing){
  double x0 = 0.7*UNKNOWN_MaxR/UNKNOWN_R0;
  double xmax = UNKNOWN_MaxR/UNKNOWN_R0;
  
  if(x<x0) return 4*M_PI*pow(UNKNOWN_R0,3)*(UNKNOWN_Rho0*(x/((1+x)*(1+x))));
  else {
    double rho0 = UNKNOWN_Rho0*(1/(x0*(1+x0)*(1+x0)));
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) ); //pow(x/x0,-10);
    return 4*M_PI*pow(UNKNOWN_R0*x,2)* UNKNOWN_R0 *rho;
  }
  //return 4*M_PI*pow(UNKNOWN_R0,3)*(UNKNOWN_Rho0*(x/((1+x)*(1+x))));
}
double UNKNOWN_calculator::set_mass(double x){
  gsl_integration_workspace * w 
  = gsl_integration_workspace_alloc (1000);

  double  error;
  double result;

  gsl_function F;
  F.function = &mass_base_UNKNOWN;
  gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
  gsl_integration_workspace_free (w);
  return result;
}
void UNKNOWN_calculator::initialize_mass(){
  
  
  Table_MassProf_r_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  Table_MassProf_M_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  Table_MassProf_rho_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  Table_MassProf_rhodx_UNKNOWN = new double [UNKNOWN_MassProfNBin];

  double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);

  for (int b=0; b<UNKNOWN_MassProfNBin; b++)
  {
    Table_MassProf_r_UNKNOWN[b] = dr*b;
    Table_MassProf_M_UNKNOWN[b] = set_mass( Table_MassProf_r_UNKNOWN[b]/UNKNOWN_R0);//main
    if(Table_MassProf_M_UNKNOWN[b]<0)cout<<Table_MassProf_M_UNKNOWN[b]<<endl;
  }

  //Rho
  for (int b=1; b<UNKNOWN_MassProfNBin; b++)
  {
    double x = dr*b/UNKNOWN_R0;
    double x0 = 0.7*UNKNOWN_MaxR/UNKNOWN_R0;
    double xmax = UNKNOWN_MaxR/UNKNOWN_R0;
  
    if(x<x0) Table_MassProf_rho_UNKNOWN[b] =UNKNOWN_Rho0*(1/(x*(1+x)*(1+x)));
    else {
      double rho0 = UNKNOWN_Rho0*(1/(x0*(1+x0)*(1+x0)));
      double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) ); //pow(x/x0,-10);
      Table_MassProf_rho_UNKNOWN[b] = rho;
    }
  }
  Table_MassProf_rho_UNKNOWN[0] = Table_MassProf_rho_UNKNOWN[1];
  //Table_MassProf_rho_UNKNOWN[UNKNOWN_MassProfNBin-1]=0;

  //Rhodx
  
  Table_MassProf_rhodx_UNKNOWN[0]=(Table_MassProf_rho_UNKNOWN[1]-Table_MassProf_rho_UNKNOWN[0])*UNKNOWN_R0/(dr);
  for (int b=1; b<UNKNOWN_MassProfNBin-1; b++)
  {
    //Table_MassProf_rhodx_UNKNOWN[b] =(Table_MassProf_rho_UNKNOWN[b+1]-Table_MassProf_rho_UNKNOWN[b-1])*UNKNOWN_R0/(2*dr);
    int num=3;
    if(b==0)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,0,num/2+1);
    else if(b==1)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,0,num/2+2);
    
    else if(b==UNKNOWN_MassProfNBin-2)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,UNKNOWN_MassProfNBin-num/2-1,UNKNOWN_MassProfNBin);
    //else if(b==UNKNOWN_MassProfNBin-1)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,UNKNOWN_MassProfNBin-num/2-1,UNKNOWN_MassProfNBin);

    else Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,b-num/2,b+num/2+1);
    Table_MassProf_rhodx_UNKNOWN[b] *= -UNKNOWN_R0;

    
  }Table_MassProf_rhodx_UNKNOWN[UNKNOWN_MassProfNBin-1]=Table_MassProf_rhodx_UNKNOWN[UNKNOWN_MassProfNBin-2];
  
}

void UNKNOWN_calculator::initialize_pot(){
  Table_MassProf_g_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  Table_MassProf_pot_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  Table_MassProf_derho_overdx_UNKNOWN = new double [UNKNOWN_MassProfNBin];
  double dr = UNKNOWN_MaxR / (UNKNOWN_MassProfNBin-1);

  Table_MassProf_g_UNKNOWN[0] =0;
  for (int b=1; b<UNKNOWN_MassProfNBin; b++)
  {
    Table_MassProf_g_UNKNOWN[b] = -UNKNOWN_NEWTON_G*Table_MassProf_M_UNKNOWN[b]/pow(Table_MassProf_r_UNKNOWN[b],2);
    
  }
  //Pot
  Table_MassProf_pot_UNKNOWN[UNKNOWN_MassProfNBin-1] = -UNKNOWN_NEWTON_G*Table_MassProf_M_UNKNOWN[UNKNOWN_MassProfNBin-1]/Table_MassProf_r_UNKNOWN[UNKNOWN_MassProfNBin-1];
  for (int b=UNKNOWN_MassProfNBin-2;b>0;b--)
  {
    Table_MassProf_pot_UNKNOWN[b] = Table_MassProf_pot_UNKNOWN[b+1] + Table_MassProf_g_UNKNOWN[b] * dr;
    
  }Table_MassProf_pot_UNKNOWN[0]=Table_MassProf_pot_UNKNOWN[1];
    
  //derho_overdx
  for (int b=0; b<UNKNOWN_MassProfNBin; b++)
  {
    Table_MassProf_derho_overdx_UNKNOWN[b] = -Table_MassProf_rhodx_UNKNOWN[b]/(Table_MassProf_g_UNKNOWN[b]*UNKNOWN_R0);
    
  }
  
}
void UNKNOWN_calculator::initialize_prob_dens(){
  double min =-potential_UNKNOWN(100),max =-potential_UNKNOWN(0.001);/***difference***/
  delta =(max-min)/size_UNKNOWN;
  double eng=min;

  
  for(int k =0;k<size_UNKNOWN;k++){
    psi[k] = eng;
    int_prob_dens[k] = integration_eng_base_UNKNOWN(eng);
    
    
    eng +=delta;
  
  }
  for(int k =0;k<size_UNKNOWN;k++){

    if(k==0)prob_dens[k]=slope(psi,int_prob_dens,k,k+5);
    else if(k==1)prob_dens[k]=slope(psi,int_prob_dens,k-1,k+4);

    else if(k==size_UNKNOWN-2)prob_dens[k]=slope(psi,int_prob_dens,k-3,k+2);
    else if(k==size_UNKNOWN-1)prob_dens[k]=slope(psi,int_prob_dens,k-4,k+1);

    else prob_dens[k]=slope(psi,int_prob_dens,k-2,k+3);

    if(prob_dens[k]<0)prob_dens[k]=0;
    
  }
  smooth_all(prob_dens,0,size_UNKNOWN);
}

void UNKNOWN_calculator::init(double newton_g,double rho,double r,int nbin,double rmax){
  Table_MassProf_r_UNKNOWN = NULL;
  Table_MassProf_M_UNKNOWN= NULL;
  Table_MassProf_g_UNKNOWN = NULL;
  Table_MassProf_pot_UNKNOWN= NULL;

  UNKNOWN_NEWTON_G=newton_g;
  UNKNOWN_Rho0=rho;
  UNKNOWN_R0=r;
  UNKNOWN_MassProfNBin=1000000;
  UNKNOWN_MaxR=rmax;

  initialize_mass();
  initialize_pot();
  initialize_prob_dens();
  
}
double UNKNOWN_calculator::set_vel(double r){  
  double index,sum=0;
  double psi_per =-potential_UNKNOWN(r);
  for(int k =0;k<size_UNKNOWN;k++){
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

  for(int k =0;k<size_UNKNOWN;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(prob_dens[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
  }
  psi_ass = psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_UNKNOWN(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_UNKNOWN(r)),0.5);
  return v;
}  
