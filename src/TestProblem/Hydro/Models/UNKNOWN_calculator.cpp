#include "UNKNOWN_calculator.h"
double UNKNOWN_NEWTON_G;
int row_UNKNOWN;

double *Table_MassProf_r_UNKNOWN;
double *Table_MassProf_M_UNKNOWN;
double *Table_MassProf_rho_UNKNOWN;
double *Table_MassProf_rhodx_UNKNOWN;
double *Table_MassProf_derho_overdx_UNKNOWN;
double *Table_MassProf_g_UNKNOWN;
double *Table_MassProf_pot_UNKNOWN;

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
  
  return Mis_InterpolateFromTable( row_UNKNOWN, Table_MassProf_r_UNKNOWN, Table_MassProf_pot_UNKNOWN, x );  
}
double de_rho_over_de_psi_UNKNOWN(double x){
  
  return Mis_InterpolateFromTable( row_UNKNOWN, Table_MassProf_r_UNKNOWN, Table_MassProf_derho_overdx_UNKNOWN, x );;
}


double psi_potential_UNKNOWN(double x, void * parameters){
  double psi= *(double *)parameters;
  return psi+potential_UNKNOWN(x);
}
double inverse_psi_to_ind_UNKNOWN (double psi) {
  int max=row_UNKNOWN-1;
  int min =0;
  int mid =(max+min)/2;
  double mid_psi;

  while(true){
    if(max==min+1){
      return max;
    }
    mid_psi=- Table_MassProf_pot_UNKNOWN[mid] ;
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
double UNKNOWN_calculator::integration_eng_base_UNKNOWN(double eng){
  double min =  eng_min_UNKNOWN;
  double max = eng;
  int num=1000;

  double dx=(max-min)/num;
  double result_right=0,result_left=0,result_simpson=0;
  double result = 0;
  for(int i=0;i<num;i++){
    double psi_l = min+i*dx,psi_r = min+(i+1)*dx;
    int ind0 = inverse_psi_to_ind_UNKNOWN(min+(i+0.5)*dx);
    if(i==num-1)result += -2* Table_MassProf_derho_overdx_UNKNOWN[ind0] * ( pow(eng-psi_l,0.5) );
    else result += -2* Table_MassProf_derho_overdx_UNKNOWN[ind0] * ( pow(eng-psi_l,0.5) - pow(eng-psi_r,0.5) );
  }
  return result;
}
double UNKNOWN_calculator::set_mass(double x){
  if(x>=Table_MassProf_r_UNKNOWN[row_UNKNOWN-1])return Table_MassProf_M_UNKNOWN[row_UNKNOWN-1];
  return Mis_InterpolateFromTable( row_UNKNOWN, Table_MassProf_r_UNKNOWN, Table_MassProf_M_UNKNOWN, x );
  
}
void UNKNOWN_calculator::initialize_mass(int row_UNKNOWN){
  Table_MassProf_M_UNKNOWN = new double [row_UNKNOWN];
  Table_MassProf_rhodx_UNKNOWN = new double [row_UNKNOWN];
  //Mass
  Table_MassProf_M_UNKNOWN[0]=0;
  double rho,dr,r;
  for (int b=1; b<row_UNKNOWN; b++)
  {
    rho = (Table_MassProf_rho_UNKNOWN[b] + Table_MassProf_rho_UNKNOWN[b-1])/2;
    dr = Table_MassProf_r_UNKNOWN[b] - Table_MassProf_r_UNKNOWN[b-1];
    r = (Table_MassProf_r_UNKNOWN[b] + Table_MassProf_r_UNKNOWN[b-1])/2;
    Table_MassProf_M_UNKNOWN[b] = Table_MassProf_M_UNKNOWN[b-1] + 4*M_PI*pow(r,2) *rho * dr;
  }

  //Rhodx
  Table_MassProf_rhodx_UNKNOWN[0]=(Table_MassProf_rho_UNKNOWN[1]-Table_MassProf_rho_UNKNOWN[0])/(Table_MassProf_r_UNKNOWN[1]-Table_MassProf_r_UNKNOWN[0]);
  for (int b=1; b<row_UNKNOWN-1; b++)
  {
    int num=3;
    if(b==0)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,0,num/2+1);
    else if(b==1)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,0,num/2+2);
    
    else if(b==row_UNKNOWN-2)Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,row_UNKNOWN-num/2-1,row_UNKNOWN);
    else Table_MassProf_rhodx_UNKNOWN[b] = slope(Table_MassProf_r_UNKNOWN,Table_MassProf_rho_UNKNOWN,b-num/2,b+num/2+1);
    

    
  }Table_MassProf_rhodx_UNKNOWN[row_UNKNOWN-1]=Table_MassProf_rhodx_UNKNOWN[row_UNKNOWN-2];
  
}

void UNKNOWN_calculator::initialize_pot(int row_UNKNOWN){
  Table_MassProf_g_UNKNOWN = new double [row_UNKNOWN];
  Table_MassProf_pot_UNKNOWN = new double [row_UNKNOWN];
  Table_MassProf_derho_overdx_UNKNOWN = new double [row_UNKNOWN];
  

  Table_MassProf_g_UNKNOWN[0] =0;
  for (int b=1; b<row_UNKNOWN; b++)
  {
    Table_MassProf_g_UNKNOWN[b] = -UNKNOWN_NEWTON_G*Table_MassProf_M_UNKNOWN[b]/pow(Table_MassProf_r_UNKNOWN[b],2); 
  }
  //Pot
  Table_MassProf_pot_UNKNOWN[row_UNKNOWN-1] = -UNKNOWN_NEWTON_G*Table_MassProf_M_UNKNOWN[row_UNKNOWN-1]/Table_MassProf_r_UNKNOWN[row_UNKNOWN-1];
  eng_min_UNKNOWN = -Table_MassProf_pot_UNKNOWN[row_UNKNOWN-1];
  for (int b=row_UNKNOWN-2;b>0;b--)
  {
    double dr = Table_MassProf_r_UNKNOWN[b+1]-Table_MassProf_r_UNKNOWN[b];
    
    Table_MassProf_pot_UNKNOWN[b] = Table_MassProf_pot_UNKNOWN[b+1] + Table_MassProf_g_UNKNOWN[b] * dr;
    
  }Table_MassProf_pot_UNKNOWN[0]=Table_MassProf_pot_UNKNOWN[1];
    
  //derho_overdx
  for (int b=0; b<row_UNKNOWN; b++)
  {
    Table_MassProf_derho_overdx_UNKNOWN[b] = -Table_MassProf_rhodx_UNKNOWN[b]/(Table_MassProf_g_UNKNOWN[b]);
  }
  
}
void UNKNOWN_calculator::initialize_prob_dens(){
  double min =-Table_MassProf_pot_UNKNOWN[row_UNKNOWN-1],max =-Table_MassProf_pot_UNKNOWN[1];
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

void UNKNOWN_calculator::init(double newton_g,int r_col,int rho_col,const char* Filename){
  Table_MassProf_r_UNKNOWN = NULL;
  Table_MassProf_M_UNKNOWN= NULL;
  Table_MassProf_g_UNKNOWN = NULL;
  Table_MassProf_pot_UNKNOWN= NULL;

  UNKNOWN_NEWTON_G=newton_g;

  int Tcol_r[1]={0};
  int Tcol_rho[1]={0};
  int r_row_UNKNOWN= Aux_LoadTable( Table_MassProf_r_UNKNOWN, "r.txt", 1, Tcol_r,true,true );
  int density_row_UNKNOWN= Aux_LoadTable( Table_MassProf_rho_UNKNOWN, "rho.txt", 1, Tcol_rho,true,true );
  
  row_UNKNOWN=r_row_UNKNOWN;
  
  initialize_mass(row_UNKNOWN);
  initialize_pot(row_UNKNOWN);
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
