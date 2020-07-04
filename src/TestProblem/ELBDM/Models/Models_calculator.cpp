#include "Models_calculator.h"
double Models_NEWTON_G;
double Models_rho;
double Models_r;
double Models_maxr;
double alpha;
int row_Models;
int  Models_massprofnbin;

double *Table_MassProf_r_Models;
double *Table_MassProf_M_Models;
double *Table_MassProf_rho_Models;
double *Table_MassProf_rhodx_Models;
double *Table_MassProf_derho_overdx_Models;
double *Table_MassProf_g_Models;
double *Table_MassProf_pot_Models;


Models_calculator::Models_calculator()
{

}

Models_calculator::~Models_calculator()
{

}
//statistics
double Models_calculator::ave(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=a[k];
  }
  return sum/(fin-start);
}
double Models_calculator::var_n(double* a,int start,int fin){
  double sum=0;
  for(int k=start;k<fin;k++){
    sum+=(a[k])*(a[k]);
  }
  sum=sum-(fin-start)*pow(ave(a,start,fin),2);
  return sum;
}
double Models_calculator::cor(double* x,double* y,int start,int fin){
  double up=0,down = pow(var_n(x,start,fin)*var_n(y,start,fin),0.5);
  double ave_x = ave(x,start,fin),ave_y = ave(y,start,fin);
  for(int k=start;k<fin;k++){
    up+=(x[k]-ave_x)*(y[k]-ave_y);
  }
  return up/down;
}
void Models_calculator::mask(double* x,int start,int fin){
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
void Models_calculator::add_num(double* x,int start,int fin){
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
void Models_calculator::smooth_all(double* x,int start,int fin){
  int num=10;
  for(int k=start;k<fin-num+1;k++){
    mask(x,k,k+num);
  }
  for(int k=start;k<fin-num+1;k++){
    add_num(x,k,k+num);
  }
}

double Models_calculator::slope(double* x,double* y,int start,int fin){
  double cor_ = cor(x,y,start,fin);
  double var_n_x =var_n(x,start,fin), var_n_y =var_n(y,start,fin);
  double s =cor_*pow(var_n_y,0.5)/pow(var_n_x,0.5);

  return s;
}

double potential_Models(double x){
  
  return Mis_InterpolateFromTable( row_Models, Table_MassProf_r_Models, Table_MassProf_pot_Models, x );  
}
double de_rho_over_de_psi_Models(double x){
  
  return Mis_InterpolateFromTable( row_Models, Table_MassProf_r_Models, Table_MassProf_derho_overdx_Models, x );;
}


double psi_potential_Models(double x, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_Models(x);
}
double inverse_psi_to_ind_Models (double psi) {
  int max=row_Models-1;
  int min =0;
  int mid =(max+min)/2;
  double mid_psi;

  while(true){
    if(max==min+1){
      return max;
    }
    mid_psi=- Table_MassProf_pot_Models[mid] ;
    if(psi>mid_psi){
      max=mid;
      mid =(max+min)/2;
    }
    else{
      min=mid;
      mid =(max+min)/2;
    }
  }

}
double inverse_psi_to_x_Models (double psi) {
  double psi1=psi;
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double x0= 0;
  double x_lo =0.0001, x_hi = 100000000000.0;
  gsl_function F;
  
  F.function = &psi_potential_Models;
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
double Models_calculator::integration_eng_base_Models(double eng){
  double min =  eng_min_Models;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  double result_right=0,result_left=0,result_simpson=0;
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    int ind0 = inverse_psi_to_ind_Models(min+(i+0.5)*dx);
    if(i==num-1)result += -2* Table_MassProf_derho_overdx_Models[ind0] * ( pow(eng-psi_l,0.5) );
    else result += -2* Table_MassProf_derho_overdx_Models[ind0] * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
  }
  return result;
}

//Calculate Mass
//Different Model Type
//Plummer
double mass_base_Plummer(double x,void* nothing){
    return 4*M_PI*pow(Models_r,3)*(Models_rho*pow(x,2)*pow(1+x,-2.5));
}
double mass_base_Plummer_trunc(double x,void* trunc_fac){
  double fac = *(double *) trunc_fac;
  double x0 = fac*Models_maxr/Models_r;
  double xmax = Models_maxr/Models_r;
  
  if(x<x0) return 4*M_PI*pow(Models_r,3)*(Models_rho*pow(x,2)*pow(1+x,-2.5));
  else {
    double rho0 = Models_rho*pow(x0,2)*pow(1+x0,-2.5);
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(Models_r*x,2)* Models_r *rho;
  }
}

//NFW
double mass_base_NFW(double x,void* nothing){
    return 4*M_PI*pow(Models_r,3)*(Models_rho*(x/((1+x)*(1+x))));
}
double mass_base_NFW_trunc(double x,void* trunc_fac){
  double fac = *(double *) trunc_fac;
  double x0 = fac*Models_maxr/Models_r;
  double xmax = Models_maxr/Models_r;
  
  if(x<x0) return 4*M_PI*pow(Models_r,3)*(Models_rho*(x/((1+x)*(1+x))));
  else {
    double rho0 = Models_rho*(1/(x0*(1+x0)*(1+x0)));
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(Models_r*x,2)* Models_r *rho;
  }
}

//Burkert
double mass_base_Burkert(double x,void* nothing){
    return 4*M_PI*pow(Models_r,3)*(Models_rho*(x*x*(1/((1+x)*(1+x*x)))));
}
double mass_base_Burkert_trunc(double x,void* trunc_fac){
  double fac = *(double *) trunc_fac;
  double x0 = fac*Models_maxr/Models_r;
  double xmax = Models_maxr/Models_r;
  
  if(x<x0) return 4*M_PI*pow(Models_r,3)*(Models_rho*(x*x*(1/((1+x)*(1+x*x)))));
  else {
    double rho0 = Models_rho*(1/((1+x0)*(1+x0*x0)));
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(Models_r*x,2)* Models_r *rho;
  }
}

//Jaffe
double mass_base_Jaffe(double x,void* nothing){
    return 4*M_PI*pow(Models_r,3)*(Models_rho*(1/(1+x)));
}
double mass_base_Jaffe_trunc(double x,void* trunc_fac){
  double fac = *(double *) trunc_fac;
  double x0 = fac*Models_maxr/Models_r;
  double xmax = Models_maxr/Models_r;
  
  if(x<x0) return 4*M_PI*pow(Models_r,3)*(Models_rho*(1/(1+x)));
  else {
    double rho0 = Models_rho*(1/(x0*x0*(1+x0)));
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(Models_r*x,2)* Models_r *rho;
  }
}

//Hernquist
double mass_base_Hernquist(double x,void* nothing){
    return 4*M_PI*pow(Models_r,3)*(Models_rho*(x/((1+x)*(1+x)*(1+x))));
}
double mass_base_Hernquist_trunc(double x,void* trunc_fac){
  double fac = *(double *) trunc_fac;
  double x0 = fac*Models_maxr/Models_r;
  double xmax = Models_maxr/Models_r;
  
  if(x<x0) return 4*M_PI*pow(Models_r,3)*(Models_rho*(x/((1+x)*(1+x)*(1+x))));
  else {
    double rho0 = Models_rho*(1/(x0*(1+x0)*(1+x0)*(1+x0)));
    double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) );
    return 4*M_PI*pow(Models_r*x,2)* Models_r *rho;
  }
}

//Einasto
double mass_base_Einasto(double x,void *nothing){
  return 4*M_PI*Models_rho*pow(Models_r,3)*pow(x,2) *exp(-pow(x,alpha));
}
double mass_base_Einasto_trunc(double x,void *nothing){
  return 4*M_PI*Models_rho*pow(Models_r,3)*pow(x,2) *exp(-pow(x,alpha));
}

double Models_calculator::set_mass(double x){
  if (model_type == "UNKNOWN"){
    if(x>=Table_MassProf_r_Models[row_Models-1])return Table_MassProf_M_Models[row_Models-1];
    return Mis_InterpolateFromTable( row_Models, Table_MassProf_r_Models, Table_MassProf_M_Models, x );
  }

  else{
    gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);

    double  error;
    double result;
    gsl_function F;

    if(Trunc_Flag){
      if(model_type=="Plummer")F.function = &mass_base_Plummer_trunc;
      else if(model_type=="NFW")F.function = &mass_base_NFW_trunc;
      else if(model_type=="Burkert")F.function = &mass_base_Burkert_trunc;
      else if(model_type=="Jaffe")F.function = &mass_base_Jaffe_trunc;
      else if(model_type=="Hernquist")F.function = &mass_base_Hernquist_trunc;
      else if(model_type=="Einasto")F.function = &mass_base_Einasto_trunc;
    
      F.params =&Trunc_Fac;
      gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
      gsl_integration_workspace_free (w);
      return result;
    }
    else {
      if(model_type=="Plummer")F.function = &mass_base_Plummer;
      else if(model_type=="NFW")F.function = &mass_base_NFW;
      else if(model_type=="Burkert")F.function = &mass_base_Burkert;
      else if(model_type=="Jaffe")F.function = &mass_base_Jaffe;
      else if(model_type=="Hernquist")F.function = &mass_base_Hernquist;
      else if(model_type=="Einasto")F.function = &mass_base_Einasto;
      gsl_integration_qag  (&F, 0, x, 0, 1e-7, 1000, 1, w, &result,  &error);
      gsl_integration_workspace_free (w);
      return result;
    }
  }
  
}
void Models_calculator::initialize_mass_UNKNOWN(int row_Models){
  
  //Mass
  Table_MassProf_M_Models[0]=0;
  double rho,dr,r;
  for (int b=1; b<row_Models; b++)
  {
    rho = (Table_MassProf_rho_Models[b] + Table_MassProf_rho_Models[b-1])/2;
    dr = Table_MassProf_r_Models[b] - Table_MassProf_r_Models[b-1];
    r = (Table_MassProf_r_Models[b] + Table_MassProf_r_Models[b-1])/2;
    Table_MassProf_M_Models[b] = Table_MassProf_M_Models[b-1] + 4*M_PI*pow(r,2) *rho * dr;
  }

  //Rhodx
  Table_MassProf_rhodx_Models[0]=(Table_MassProf_rho_Models[1]-Table_MassProf_rho_Models[0])/(Table_MassProf_r_Models[1]-Table_MassProf_r_Models[0]);
  for (int b=1; b<row_Models-1; b++)
  {
    int num=3;
    if(b==0)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,0,num/2+1);
    else if(b==1)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,0,num/2+2);
    
    else if(b==row_Models-2)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,row_Models-num/2-1,row_Models);
    else Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,b-num/2,b+num/2+1);
    

    
  }Table_MassProf_rhodx_Models[row_Models-1]=Table_MassProf_rhodx_Models[row_Models-2];
  
}

void Models_calculator::initialize_pot_UNKNOWN(int row_Models){
  
  

  Table_MassProf_g_Models[0] =0;
  for (int b=1; b<row_Models; b++)
  {
    Table_MassProf_g_Models[b] = -Models_NEWTON_G*Table_MassProf_M_Models[b]/pow(Table_MassProf_r_Models[b],2); 
  }
  //Pot
  Table_MassProf_pot_Models[row_Models-1] = -Models_NEWTON_G*Table_MassProf_M_Models[row_Models-1]/Table_MassProf_r_Models[row_Models-1];
  eng_min_Models = -Table_MassProf_pot_Models[row_Models-1];
  for (int b=row_Models-2;b>0;b--)
  {
    double dr = Table_MassProf_r_Models[b+1]-Table_MassProf_r_Models[b];
    
    Table_MassProf_pot_Models[b] = Table_MassProf_pot_Models[b+1] + Table_MassProf_g_Models[b] * dr;
    
  }Table_MassProf_pot_Models[0]=Table_MassProf_pot_Models[1];
    
  //derho_overdx
  for (int b=0; b<row_Models; b++)
  {
    Table_MassProf_derho_overdx_Models[b] = -Table_MassProf_rhodx_Models[b]/(Table_MassProf_g_Models[b]);
  }
  
}

void Models_calculator::initialize_mass_others(){
  
  double dr = Models_maxr / (Models_massprofnbin-1);

  for (int b=0; b<Models_massprofnbin; b++)
  {
    Table_MassProf_r_Models[b] = dr*b;
    Table_MassProf_M_Models[b] = set_mass( Table_MassProf_r_Models[b]/Models_r);
  }
  
  //Rho
  if(Trunc_Flag){
    for (int b=1; b<Models_massprofnbin; b++)
    {
      double x = dr*b/Models_r;
      double x0 = Trunc_Fac*Models_maxr/Models_r;
      double xmax = Models_maxr/Models_r;
    
      if(x<x0) Table_MassProf_rho_Models[b] =Models_rho*(1/(x*x*(1+x)));
      else {
        double rho0 = Models_rho*(1/(x0*x0*(1+x0)));
        double rho = rho0 *( 1 - pow(1-pow((x-xmax)/(xmax-x0) ,2) ,0.5 ) ); 
        Table_MassProf_rho_Models[b] = rho;
      }
    }
    Table_MassProf_rho_Models[0] = Table_MassProf_rho_Models[1];
  }
  else{
    for (int b=1; b<Models_massprofnbin; b++)
    {
      double x = dr*b/Models_r;
      Table_MassProf_rho_Models[b] =Models_rho*(1/(x*x*(1+x)));
    }
    Table_MassProf_rho_Models[0] = Table_MassProf_rho_Models[1];
  }

  //Rhodx
  Table_MassProf_rhodx_Models[0]=(Table_MassProf_rho_Models[1]-Table_MassProf_rho_Models[0])*Models_r/(dr);
  for (int b=1; b<Models_massprofnbin-1; b++)
  {
    int num=3;
    if(b==0)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,0,num/2+1);
    else if(b==1)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,0,num/2+2);
    
    else if(b==Models_massprofnbin-2)Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,Models_massprofnbin-num/2-1,Models_massprofnbin);
    else Table_MassProf_rhodx_Models[b] = slope(Table_MassProf_r_Models,Table_MassProf_rho_Models,b-num/2,b+num/2+1);
    Table_MassProf_rhodx_Models[b] *= -Models_r;

    
  }Table_MassProf_rhodx_Models[Models_massprofnbin-1]=Table_MassProf_rhodx_Models[Models_massprofnbin-2];
  
}

void Models_calculator::initialize_pot_others(){

  double dr = Models_maxr / (Models_massprofnbin-1);

  Table_MassProf_g_Models[0] =0;
  for (int b=1; b<Models_massprofnbin; b++)
  {
    Table_MassProf_g_Models[b] = -Models_NEWTON_G*Table_MassProf_M_Models[b]/pow(Table_MassProf_r_Models[b],2);
    
  }
  //Pot
  Table_MassProf_pot_Models[Models_massprofnbin-1] = -Models_NEWTON_G*Table_MassProf_M_Models[Models_massprofnbin-1]/Table_MassProf_r_Models[Models_massprofnbin-1];
  eng_min_Models = -Table_MassProf_pot_Models[Models_massprofnbin-1];
  for (int b=Models_massprofnbin-2;b>0;b--)
  {
    Table_MassProf_pot_Models[b] = Table_MassProf_pot_Models[b+1] + Table_MassProf_g_Models[b] * dr;
  }Table_MassProf_pot_Models[0]=Table_MassProf_pot_Models[1];
    
  //derho_overdx
  for (int b=0; b<Models_massprofnbin; b++)
  {
    Table_MassProf_derho_overdx_Models[b] = -Table_MassProf_rhodx_Models[b]/(Table_MassProf_g_Models[b]*Models_r);
  }
  
}
void Models_calculator::initialize_prob_dens(){
  double min,max;
  if(model_type=="UNKNOWN")min=-Table_MassProf_pot_Models[row_Models-1];
  else min = -Table_MassProf_pot_Models[Models_massprofnbin-1];
  max =-Table_MassProf_pot_Models[1];
  delta =(max-min)/size_Models;
  double eng=min;

  

  for(int k =0;k<size_Models;k++){
    psi[k] = eng;
    int_prob_dens[k] = integration_eng_base_Models(eng);
    
    eng +=delta;
  }
  for(int k =0;k<size_Models;k++){

    if(k==0)prob_dens[k]=slope(psi,int_prob_dens,k,k+5);
    else if(k==1)prob_dens[k]=slope(psi,int_prob_dens,k-1,k+4);

    else if(k==size_Models-2)prob_dens[k]=slope(psi,int_prob_dens,k-3,k+2);
    else if(k==size_Models-1)prob_dens[k]=slope(psi,int_prob_dens,k-4,k+1);

    else prob_dens[k]=slope(psi,int_prob_dens,k-2,k+3);

    if(prob_dens[k]<0)prob_dens[k]=0;
    
  }
  smooth_all(prob_dens,0,size_Models);
}

void Models_calculator::init(string type,double al,double newton_g,double rho,double r,int nbin,double rmax,int rseed,bool trunc_flag,double trunc_fac,int r_col,int rho_col,const char* Filename){
  Table_MassProf_r_Models = NULL;
  Table_MassProf_M_Models= NULL;
  Table_MassProf_g_Models = NULL;
  Table_MassProf_pot_Models= NULL;

  model_type = type;
  alpha = al;
  Models_NEWTON_G=newton_g;
  Models_rho=rho;
  Models_r=r;
  Models_massprofnbin=nbin;
  Models_maxr=rmax;

  Trunc_Flag=trunc_flag;
  Trunc_Fac=trunc_fac;

  //Set random seeds
  if(RNG==NULL)RNG = new RandomNumber_t( 1 );
  RNG->SetSeed( 0, rseed );

  //Initialize densities with Table
  if(model_type=="UNKNOWN"){
    int Tcol_r[1]={r_col};
    int Tcol_rho[1]={rho_col};
    int r_row_Models;
    if(sizeof(Table_MassProf_r_Models)==0)r_row_Models= Aux_LoadTable( Table_MassProf_r_Models, Filename, 1, Tcol_r,true,true );
    else r_row_Models= Aux_LoadTable( Table_MassProf_r_Models, Filename, 1, Tcol_r,true,false );
    
    int density_row_Models;
    if(sizeof(Table_MassProf_rho_Models)==0)density_row_Models= Aux_LoadTable( Table_MassProf_rho_Models, Filename, 1, Tcol_rho,true,true );
    else density_row_Models= Aux_LoadTable( Table_MassProf_rho_Models, Filename, 1, Tcol_rho,true,false );
    
    row_Models=r_row_Models;

    Table_MassProf_M_Models = new double [row_Models];
    Table_MassProf_rhodx_Models = new double [row_Models];
    Table_MassProf_g_Models = new double [row_Models];
    Table_MassProf_pot_Models = new double [row_Models];
    Table_MassProf_derho_overdx_Models = new double [row_Models];

    initialize_mass_UNKNOWN(row_Models);
    initialize_pot_UNKNOWN(row_Models);
  }

  else{
    Table_MassProf_r_Models = new double [Models_massprofnbin];
    Table_MassProf_M_Models = new double [Models_massprofnbin];
    Table_MassProf_rho_Models = new double [Models_massprofnbin];
    Table_MassProf_rhodx_Models = new double [Models_massprofnbin];
    Table_MassProf_g_Models = new double [Models_massprofnbin];
    Table_MassProf_pot_Models = new double [Models_massprofnbin];
    Table_MassProf_derho_overdx_Models = new double [Models_massprofnbin];

    initialize_mass_others();
    initialize_pot_others();
  }


  initialize_prob_dens();
}
double Models_calculator::set_vel(double r){  
  double index,sum=0;
  double psi_per =-potential_Models(r);
  for(int k =0;k<size_Models;k++){
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

  for(int k =0;k<size_Models;k++){
    if(sum_mes>sum_rad){
      index_ass =k-1;
      par = (sum_mes-sum_rad)/(prob_dens[index_ass] *pow(psi_per-psi[index_ass],0.5) *delta);
      break;
      }
    sum_mes += prob_dens[k] *pow(psi_per-psi[k],0.5) *delta;
  }
  psi_ass = psi[index_ass] +delta *par;
  if(-2*(psi_ass+potential_Models(r))<0){
    return 0;
  }
  double v =pow(-2*(psi_ass+potential_Models(r)),0.5);
  return v;
}  

double Models_calculator::set_vel_test(double r){  
  const double TotM_Inf    = 4.0/3.0*M_PI*CUBE(Models_r)*Models_rho;
  const double Vmax_Fac    = sqrt( 2.0*Models_NEWTON_G*TotM_Inf );

  double  Vmax, RanV, RanProb, Prob;

  Vmax = Vmax_Fac*pow( SQR(Models_r) + SQR(r*Models_r), -0.25 );

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
