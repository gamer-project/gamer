#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

#include "Par_EquilibriumIC.h"



Par_EquilibriumIC::Par_EquilibriumIC()
{
}

Par_EquilibriumIC::~Par_EquilibriumIC()
{
}



//-------------------------------------------------------------------------------------------------------
// Function    :  Read_Filenames
// Description :  Read in file names of physical parameters
//
// Note        :
//
// Parameter   :  filename_para : pointer of the file name
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Read_Filenames( const char *filename_para )
{

   vector <string> EMPTY;

   filenames.Cloud_Num = GetParams( filename_para, "ParEqmIC_Cloud_Num", 1, "int", EMPTY );

   GetParams( filename_para, "ParEqmIC_Params_Filenames", filenames.Cloud_Num, "string", filenames.Params_Filenames );


   Aux_Message( stdout, "Checking Params_Filenames\n" );

   fstream file;

   for (int k=0; k<filenames.Cloud_Num; k++)
   {

      const char * c = filenames.Params_Filenames[k].c_str();

      file.open( c, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "Test Problem parameter file %s cannot be found !!\n", filenames.Params_Filenames[k].c_str() );
         Aux_Error( ERROR_INFO, "Error in the input of Params_Filenames !!\n" );
      }

      file.close();
   }

} // FUNCTION : Read_Filenames



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Physical_Params
// Description :  Load the physical parameters from a file
//
// Note        :
//
// Parameter   :  filename_para : the file name loader's FP (file parameters) structure
//                cloud_idx     : index of the cloud
//                NPar_AllRank  : particle total number (including all clouds)
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Load_Physical_Params( const FP filename_para, const int cloud_idx, const long NPar_AllRank )
{

   params.Cloud_Center   = new double[3]; // central coordinates
   params.Cloud_BulkVel  = new double[3]; // bulk velocity

   Aux_Message( stdout, "Reading physical parameters input file: %s\n", filename_para.Params_Filenames[cloud_idx].c_str() );

   // (1) load the problem-specific runtime parameters
   params.Cloud_Num        = filename_para.Cloud_Num;
   params.Params_Filenames = filename_para.Params_Filenames[cloud_idx];

   const char* FileName  = filename_para.Params_Filenames[cloud_idx].c_str();
   ReadPara_t *ReadPara  = new ReadPara_t;
   double ratio;

   // (1-1) add parameters in the following format:
   // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
   // --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,                          DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "Cloud_RSeed",                &params.Cloud_RSeed,                123,           0,                NoMax_int         );
   ReadPara->Add( "Cloud_Rho0",                 &params.Cloud_Rho0,                 1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_R0",                   &params.Cloud_R0,                   0.1,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_MaxR",                 &params.Cloud_MaxR,                 0.375,         Eps_double,       NoMax_double      );
   ReadPara->Add( "Cloud_CenterX",              &params.Cloud_Center[0],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_CenterY",              &params.Cloud_Center[1],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_CenterZ",              &params.Cloud_Center[2],            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelX",             &params.Cloud_BulkVel[0],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelY",             &params.Cloud_BulkVel[1],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_BulkVelZ",             &params.Cloud_BulkVel[2],           0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Cloud_MassProfNBin",         &params.Cloud_MassProfNBin,         1000,          2,                NoMax_int         );
   ReadPara->Add( "Cloud_Par_Num_Ratio",        &ratio,                             0.,            0.,               1.0               );
   ReadPara->Add( "Cloud_Type",                  params.Cloud_Type,                 NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "Density_Table_Name",          params.Density_Table_Name,         NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "AddExtPot",                  &params.AddExtPot,                  0,             0,                1                 );
   ReadPara->Add( "ExtPot_Table_Name",           params.ExtPot_Table_Name,          NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "Cloud_Einasto_Power_Factor", &params.Cloud_Einasto_Power_Factor, 1.0,           0.1,              10.0              );

   ReadPara->Read( FileName );
   delete ReadPara;

   // Convert Cloud_Par_Num_Ratio to Cloud_Par_Num
   params.Cloud_Par_Num = long(ratio*NPar_AllRank);

   // Check whether user forgot to fill in Cloud_Par_Num_Ratio
   if ( params.Cloud_Par_Num == 0 )
   {
      Aux_Error( ERROR_INFO, "Cloud_Par_Num_Ratio is 0! There is no particle in this cloud!!" );
   }

   // (1-2) set the default values
   for (int d=0; d<3; d++)
      if ( params.Cloud_Center[d] == NoDef_double )  params.Cloud_Center[d] = 0.5*amr->BoxSize[d];

   // (2) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  random seed for setting particle position = %d\n",     params.Cloud_RSeed );
      Aux_Message( stdout, "  peak density                              = %13.7e\n", params.Cloud_Rho0 );
      Aux_Message( stdout, "  scale radius                              = %13.7e\n", params.Cloud_R0 );
      Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", params.Cloud_MaxR );

      for (int d=0; d<3; d++)
      {
      Aux_Message( stdout, "  central coordinate [%d]                    = %13.7e\n", d, params.Cloud_Center[d] );
      Aux_Message( stdout, "  bulk velocity [%d]                         = %13.7e\n", d, params.Cloud_BulkVel[d] );
      }

      if ( strcmp( params.Cloud_Type, "Table" ) != 0 )
      Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     params.Cloud_MassProfNBin );

      Aux_Message( stdout, "  Cloud_Type                                = %s\n",     params.Cloud_Type );
      Aux_Message( stdout, "=============================================================================\n" );
   }//if ( MPI_Rank == 0 )

   // (3) Warn against small R0
   if ( params.Cloud_R0 < amr->dh[MAX_LEVEL] )
      Aux_Message( stdout, "WARNING : Characteristic length R0:%f is smaller than spatial resolution %f!\n", params.Cloud_R0, amr->dh[MAX_LEVEL] );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

   // (4) Check Cloud_Type and table filenames
   // Checking Cloud_Type
   Aux_Message( stdout, "Checking Cloud_Type\n" );

   int flag = 0;

   if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 )   flag = 1;
   else if ( strcmp( params.Cloud_Type, "Table"     ) == 0 )   flag = 1;

   if ( flag == 0 )
   {
      Aux_Message( stdout, "%s is not a Model Type\n", params.Cloud_Type );
      Aux_Error( ERROR_INFO, "Error in the input of Cloud_Type !!\n" );
   }

   // Checking Density_Table_Name
   Aux_Message( stdout, "Checking Density_Table_Name\n" );

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      fstream file;
      file.open( params.Density_Table_Name, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "Density profile %s cannot be found !!\n", params.Density_Table_Name );
         Aux_Error( ERROR_INFO, "Error in the input of Density_Table_Name !!\n" );
      }

      file.close();
   }

   // Checking ExtPot_Table_Name
   Aux_Message( stdout, "Checking ExtPot_Table_Name\n" );

   if ( params.AddExtPot )
   {

      fstream file;
      file.open( params.ExtPot_Table_Name, ios::in );

      if ( !file )
      {
         Aux_Message( stdout, "External potential profile %s cannot be found !!\n", params.ExtPot_Table_Name );
         Aux_Error( ERROR_INFO, "Error in the input of ExtPot_Table_Name!!\n" );
      }

      file.close();
   }

} // FUNCTION : Load_Physical_Params



//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialize all necessary tables of physical parameters, including radius, mass, density, gravitational potential
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init()
{

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "Must enable SUPPORT_GSL for Par_EquilibriumIC !!\n" );
#  endif

   Table_r                 = NULL;
   Table_Enclosed_Mass     = NULL;
   Table_Density           = NULL;
   Table_dRho_dr           = NULL;
   Table_dRho_dx           = NULL;
   Table_Gravity_Field     = NULL;
   Table_Gravity_Potential = NULL;

   prob_dens               = NULL;
   int_prob_dens           = NULL;
   psi                     = NULL;

   //Set random seeds
   Random_Num_Gen = new RandomNumber_t( 1 );
   Random_Num_Gen->SetSeed( 0, params.Cloud_RSeed );

   //Initialize densities with Table
   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {

      int Tcol_r[1]   =  {0};
      int Tcol_rho[1] =  {1};
      int Row_r_Table;

      Aux_Message( stdout, "Loading Density Profile Table: %s\n", params.Density_Table_Name );

      Row_r_Table = Aux_LoadTable( Table_r, params.Density_Table_Name, 1, Tcol_r,true,true );

      int Row_Density_Table;
      Row_Density_Table = Aux_LoadTable( Table_Density, params.Density_Table_Name , 1, Tcol_rho,true,true );

      if ( Row_r_Table != Row_Density_Table )
         Aux_Error( ERROR_INFO, "Density row number is not equal to radius row number in the profile file !! Please check this file.\n" );

      params.Cloud_MassProfNBin = Row_r_Table;

      // Radii in the density table must be no less than Cloud_MaxR
      if ( Table_r[params.Cloud_MassProfNBin-1] < params.Cloud_MaxR )
      {
         Aux_Error( ERROR_INFO, "Maximum radius in your density table is smaller then Cloud_MaxR! Please check!\n" );
      }

      Table_Enclosed_Mass     = new double [params.Cloud_MassProfNBin];
      Table_dRho_dr           = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Field     = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Potential = new double [params.Cloud_MassProfNBin];
      Table_dRho_dx           = new double [params.Cloud_MassProfNBin];

      prob_dens               = new double [params.Cloud_MassProfNBin];
      int_prob_dens           = new double [params.Cloud_MassProfNBin];
      psi                     = new double [params.Cloud_MassProfNBin];

      Init_Mass_Table();
      Init_Pot_Table ();
      Add_Ext_Pot    ();
      Init_Prob_Dens ();

   }

   else
   {
      Table_r                 = new double [params.Cloud_MassProfNBin];
      Table_Enclosed_Mass     = new double [params.Cloud_MassProfNBin];
      Table_Density           = new double [params.Cloud_MassProfNBin];
      Table_dRho_dr           = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Field     = new double [params.Cloud_MassProfNBin];
      Table_Gravity_Potential = new double [params.Cloud_MassProfNBin];
      Table_dRho_dx           = new double [params.Cloud_MassProfNBin];

      prob_dens               = new double [params.Cloud_MassProfNBin];
      int_prob_dens           = new double [params.Cloud_MassProfNBin];
      psi                     = new double [params.Cloud_MassProfNBin];

      Init_Mass     ();
      Init_Pot      ();
      Add_Ext_Pot   ();
      Init_Prob_Dens();
   }

} // FUNCTION : Init



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SetEquilibriumIC
// Description :  Set particle's initial conditions (IC) for a cloud that is in equilibrium state
//
// Note        :
//
// Parameter   :  Mass_AllRank : An array of all particles' masses
//                Pos_AllRank  : An array of all particles' position vectors
//                Vel_AllRank  : An array of all particles' velocity vectors
//                Par_Idx0     : Starting index of particles in this cloud
//
// Return      :  Mass_AllRank
//                Pos_AllRank
//                Vel_AllRank
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Par_SetEquilibriumIC( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx0 )
{

   double *Table_MassProf_r = NULL;
   double *Table_MassProf_M = NULL;
   double  TotM, ParM, dr, RanM, RanR, EstM, ErrM, ErrM_Max=-1.0, RanVec[3];
   double  Vmax, RanV, RanProb, Prob;


   // determine the total enclosed mass within the maximum radius
   TotM = Set_Mass( params.Cloud_MaxR );
   ParM = TotM / (params.Cloud_Par_Num);

   // construct the mass profile table
   Table_MassProf_r = new double [params.Cloud_MassProfNBin];
   Table_MassProf_M = new double [params.Cloud_MassProfNBin];

   dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = Set_Mass( Table_MassProf_r[b] );
   }

   // set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+params.Cloud_Par_Num; p++)
   {
      // mass
      Mass_AllRank[p] = ParM;

      // position
      // --> sample from the cumulative mass profile with linear interpolation
      RanM = Random_Num_Gen->GetValue( 0, 0.0, 1.0 )*TotM;
      RanR = Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_MassProf_M, Table_MassProf_r, RanM );

      // record the maximum error
      EstM     = Set_Mass( RanR );
      ErrM     = fabs( (EstM-RanM)/RanM );
      ErrM_Max = fmax( ErrM, ErrM_Max );

      // randomly set the position vector with a given radius
      RanVec_FixRadius( RanR, RanVec );
      for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + params.Cloud_Center[d];

      // check periodicity
      for (int d=0; d<3; d++)
      {
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
            Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
      }

      // velocity
      double a3 = RanR/params.Cloud_R0;

      RanV = Set_Velocity(a3);

      // randomly set the velocity vector with the given amplitude (RanV*Vmax)
      RanVec_FixRadius( RanV, RanVec );
      for (int d=0; d<3; d++)    Vel_AllRank[d][p] = RanVec[d] + params.Cloud_BulkVel[d];

   } // for (long p=0; p<NPar_AllRank; p++)

   Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotM );
   Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
   Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max );

   // free memory
   delete [] Table_MassProf_r;
   delete [] Table_MassProf_M;

} // FUNCTION : Par_SetEquilibriumIC


// Parameters for the intergration of mass profile
struct mass_integrand_params
{
   double Cloud_R0;
   double Cloud_Rho0;
};


struct mass_integrand_params_Einasto
{
   double Cloud_R0;
   double Cloud_Rho0;
   double Cloud_Einasto_Power_Factor;
};


//Different Model Type
//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_Plummer
// Description :  Analytical density profile of the Plummer model
//
// Note        :  1. \rho(r) = \rho_0 ( 1 + (\frac{r}{a})^2 )^{-5/2}
//                2. Reference: Plummer H. C., MNRAS, 1991, doi:10.1093/mnras/71.5.460
//
// Parameter   :  r    : input radius
//                R0   : Plummer scale radius, a
//                Rho0 : Plummer scale density, \rho_0
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_Plummer( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0*pow( 1+x*x, -2.5 );

} // FUNCTION : AnalyticalDensProf_Plummer


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_Plummer
// Description :  Integrand for the enclosed mass profile of the Plummer model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_Plummer( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params *p = (struct mass_integrand_params *) parameters;
   double R0   = p->Cloud_R0;
   double Rho0 = p->Cloud_Rho0;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_Plummer( r, R0, Rho0 );

} // FUNCTION : MassIntegrand_Plummer


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Plummer
// Description :  Analytical enclosed mass profile of the Plummer model
//
// Note        :  1. M(r) = M_0 \frac{ r^3 }{ (r^2 + a^2)^{3/2} }
//                        = \frac{4\pi}{3} a^3 \rho_0 \frac{ r^3 }{ (r^2 + a^2)^{3/2} }
//                2. Reference: Plummer H. C., MNRAS, 1991, doi:10.1093/mnras/71.5.460
//
// Parameter   :  r    : input radius
//                R0   : Plummer scale radius, a
//                Rho0 : Plummer scale density, \rho_0
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalMassProf_Plummer( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return (4.0/3.0)*M_PI*Rho0*CUBE(r)*pow( 1+x*x, -1.5 );

} // FUNCTION : AnalyticalMassProf_Plummer


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_NFW
// Description :  Analytical density profile of the NFW model
//
// Note        :  1. \rho(r) = \rho_s \frac{1}{ ( \frac{r}{r_s} ) ( 1 + \frac{r}{r_s} )^2 }
//                2. Reference: Navarro J.~F., Frenk C.~S., White S.~D.~M., 1996, ApJ, doi:10.1086/177173
//                              Li P. et al., 2020, ApJS, doi:10.3847/1538-4365/ab700e
//
// Parameter   :  r    : input radius
//                R0   : NFW scale radius, r_s
//                Rho0 : NFW scale density, \rho_s
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_NFW( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0/( x*SQR( 1+x ) );

} // FUNCTION : AnalyticalDensProf_NFW


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_NFW
// Description :  Integrand for the enclosed mass profile of the NFW model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_NFW( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params *p = (struct mass_integrand_params *) parameters;
   double R0   = p->Cloud_R0;
   double Rho0 = p->Cloud_Rho0;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_NFW( r, R0, Rho0 );

} // FUNCTION : MassIntegrand_NFW


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_NFW
// Description :  Analytical enclosed mass profile of the NFW model
//
// Note        :  1. M(r) = 4\pi \rho_s r_s^3 [ \ln( \frac{ r_s + r }{ r_s } ) - \frac{ r }{ r_s + r } ]
//                2. Reference: Navarro J.~F., Frenk C.~S., White S.~D.~M., 1996, ApJ, doi:10.1086/177173
//                              Li P. et al., 2020, ApJS, doi:10.3847/1538-4365/ab700e
//
// Parameter   :  r    : input radius
//                R0   : NFW scale radius, r_s
//                Rho0 : NFW scale density, \rho_s
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalMassProf_NFW( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return 4*M_PI*Rho0*CUBE(R0)*( log( 1+x ) - x/(1+x) );

} // FUNCTION : AnalyticalMassProf_NFW


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_Burkert
// Description :  Analytical density profile of the Burkert model
//
// Note        :  1. \rho(r) = \rho_s \frac{1}{ ( 1 + \frac{r}{r_s} ) ( 1 + (\frac{r}{r_s})^2 ) }
//                2. Reference: Burkert A., 1995, ApJL, doi:10.1086/309560
//                              Li P. et al., 2020, ApJS, doi:10.3847/1538-4365/ab700e
//
// Parameter   :  r    : input radius
//                R0   : Plummer scale radius, a
//                Rho0 : Plummer scale density, \rho_0
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_Burkert( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0/( (1+x)*(1+x*x) );

} // FUNCTION : AnalyticalDensProf_Burkert


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_Burkert
// Description :  Integrand for the enclosed mass profile of the Burkert model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_Burkert( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params *p = (struct mass_integrand_params *) parameters;
   double R0   = p->Cloud_R0;
   double Rho0 = p->Cloud_Rho0;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_Burkert( r, R0, Rho0 );

} // FUNCTION : MassIntegrand_Burkert


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Burkert
// Description :  Analytical enclosed mass profile of the Burkert model
//
// Note        :  1. M(r) = 2\pi \rho_s r_s^3 [ \frac{1}{2}\ln( 1 + (\frac{r}{r_s})^2 ) + \ln( 1 + \frac{r}{r_s} ) - \arctan( \frac{r}{r_s} ) ]
//                2. Reference: Burkert A., 1995, ApJL, doi:10.1086/309560
//                              Li P. et al., 2020, ApJS, doi:10.3847/1538-4365/ab700e
//
// Parameter   :  r    : input radius
//                R0   : Burkert scale radius, r_s
//                Rho0 : Burkert scale density, \rho_s
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalMassProf_Burkert( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return 2*M_PI*Rho0*CUBE(R0)*( 0.5*log( 1+x*x ) + log( 1+x ) - atan( x ) );

} // FUNCTION : AnalyticalMassProf_Burkert


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_Jaffe
// Description :  Analytical density profile of the Jaffe model
//
// Note        :  1. \rho(r) = \frac{ \rho_0 }{ 4\pi } \frac{ r_J^4 }{ r^2 (r_J + r)^2}
//                           = \frac{ \rho_0 }{ 4\pi } \frac{1}{ ( \frac{r}{r_J} )^2 ( 1 + \frac{r}{r_J} )^2 }
//                2. Reference: Jaffe W., 1983, MNRAS, doi:10.1093/mnras/202.4.995
//                              Ciotti L. and Ziaee Lorzad A., 2018, MNRAS, doi:10.1093/mnras/stx2771
//
// Parameter   :  r    : input radius
//                R0   : Jaffe scale radius, r_J
//                Rho0 : Jaffe scale density, \rho_0
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_Jaffe( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0/( 4*M_PI*SQR(x)*SQR(1+x) ); // return Rho0/(x*(1+x)); //previous one

} // FUNCTION : AnalyticalDensProf_Jaffe


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_Jaffe
// Description :  Integrand for the enclosed mass profile of the Jaffe model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_Jaffe( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params *p = (struct mass_integrand_params *) parameters;
   double R0   = p->Cloud_R0;
   double Rho0 = p->Cloud_Rho0;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_Jaffe( r, R0, Rho0 );

} // FUNCTION : MassIntegrand_Jaffe


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Jaffe
// Description :  Analytical enclosed mass profile of the Jaffe model
//
// Note        :  1. M(r) = M_J \frac{ r }{ r_J + r }
//                   ,where M_J = \rho_0*r_J^3 is the total mass
//                2. Reference: Jaffe W., 1983, MNRAS, doi:10.1093/mnras/202.4.995
//                              Ciotti L. and Ziaee Lorzad A., 2018, MNRAS, doi:10.1093/mnras/stx2771
//
// Parameter   :  r    : input radius
//                R0   : Jaffe scale radius, r_J
//                Rho0 : Jaffe scale density, \rho_0
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalMassProf_Jaffe( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0*CUBE(R0)*x/(1+x);

} // FUNCTION : AnalyticalMassProf_Jaffe


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_Hernquist
// Description :  Analytical density profile of the Hernquist model
//
// Note        :  1. \rho(r) = \rho_0 \frac{ a^4 }{ r (r+a)^3 }
//                           = \rho_0 \frac{1}{ \frac{r}{a} ( 1+\frac{r}{a} )^3 }
//                2. Reference: Hernquist L., 1990, ApJ, doi:10.1086/168845
//
// Parameter   :  r    : input radius
//                R0   : Hernquist scale radius, a
//                Rho0 : Hernquist scale density, \rho_0
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_Hernquist( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return Rho0/( x*CUBE( 1+x ) );

} // FUNCTION : AnalyticalDensProf_Hernquist


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_Hernquist
// Description :  Integrand for the enclosed mass profile of the Hernquist model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_Hernquist( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params *p = (struct mass_integrand_params *) parameters;
   double R0   = p->Cloud_R0;
   double Rho0 = p->Cloud_Rho0;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_Hernquist( r, R0, Rho0 );

} // FUNCTION : MassIntegrand_Hernquist


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Hernquist
// Description :  Analytical enclosed mass profile of the Hernquist model
//
// Note        :  1. M(r) = M_0 \frac{ r^2 }{ (r+a)^2 }
//                   ,where M_0 = 2\pi*\rho_0*a^3 is the total mass
//                2. Reference: Hernquist L., 1990, ApJ, doi:10.1086/168845
//
// Parameter   :  r    : input radius
//                R0   : Hernquist scale radius, a
//                Rho0 : Hernquist scale density, \rho_0
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalMassProf_Hernquist( const double r, const double R0, const double Rho0 )
{
   const double x = r/R0;

   return 2*M_PI*Rho0*CUBE(R0)*SQR(x)/SQR(1+x);

} // FUNCTION : AnalyticalMassProf_Hernquist


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalDensProf_Einasto
// Description :  Analytical density profile of the Einasto model
//
// Note        :  1. \rho(r) = \rho_s    \exp{ -d_n [ (\frac{r}{r_s})^{1/n}    - 1 ] }, where r_s is the radius at which contains half of the total mass
//                           = \rho_{-2} \exp{ -2n  [ (\frac{r}{r_{-2}})^{1/n} - 1 ] }, where r_{-2} is the radius at which \rho(r) \propto r^{-2}
//                           = \rho_0    \exp{ -[ (\frac{r}{h})^{1/n} ] }
//                2. Reference: Einasto J., 1965, TrAlm
//                              Retana-Montenegro E. et al., 2012, A&A, doi:10.1051/0004-6361/201118543
//
// Parameter   :  r                    : input radius
//                R0                   : Einasto scale radius, h = \frac{r_s}{d_n^n} = \frac{r_{-2}}{(2n)^n}
//                Rho0                 : Einasto central density, \rho_0 = \rho_s\exp{d_n} = \rho_{-2}\exp{2n}
//                Einasto_Power_Factor : Einasto power factor, 1/n
//
// Return      :  density at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalDensProf_Einasto( const double r, const double R0, const double Rho0, const double Einasto_Power_Factor )
{
   const double x = r/R0;

   return Rho0*exp( -pow( x, Einasto_Power_Factor ) );

} // FUNCTION : AnalyticalDensProf_Einasto


//-------------------------------------------------------------------------------------------------------
// Function    :  MassIntegrand_Einasto
// Description :  Integrand for the enclosed mass profile of the Einasto model
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
//-------------------------------------------------------------------------------------------------------
double MassIntegrand_Einasto( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct mass_integrand_params_Einasto *p = (struct mass_integrand_params_Einasto *) parameters;
   double R0                           = p->Cloud_R0;
   double Rho0                         = p->Cloud_Rho0;
   double Einasto_Power_Factor         = p->Cloud_Einasto_Power_Factor;

   return 4*M_PI*SQR(r)*AnalyticalDensProf_Einasto( r, R0, Rho0, Einasto_Power_Factor );

} // FUNCTION : MassIntegrand_Einasto


//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Einasto
// Description :  Analytical enclosed mass profile of the Einasto model
//
// Note        :  1. M(r) = M_0 [1 - \frac{ \Gamma( 3n, (r/h)^{1/n} ) }{ \Gamma(3n) }]
//                   ,where M_0 = 4\pi \rho_0 h^3 n \Gamma(3n) is the total mass
//                2. Reference: Einasto J., 1965, TrAlm
//                              Retana-Montenegro E. et al., 2012, A&A, doi:10.1051/0004-6361/201118543
//
// Parameter   :  r                    : input radius
//                R0                   : Einasto scale radius, h
//                Rho0                 : Einasto central density, \rho_0
//                Einasto_Power_Factor : Einasto power factor, 1/n
//
// Return      :  enclosed mass at the given radius
//-------------------------------------------------------------------------------------------------------
// double AnalyticalMassProf_Einasto( const double r, const double R0, const double Rho0, const double Einasto_Power_Factor )
// {
//    const double x = r/R0;
//
//    //Gamma function:                  Gamma( s    ) = \int_{0}^{\infty} t^{s-1}*e^{-t} dt
//    //Upper incomplete Gamma function: Gamma( s, x ) = \int_{x}^{\infty} t^{s-1}*e^{-t} dt
//
//    return 4*M_PI*Rho0*CUBE(R0)/Einasto_Power_Factor*( Gamma(3/Einasto_Power_Factor) - UpperIncompleteGamma( 3/Einasto_Power_Factor, pow( x, Einasto_Power_Factor ) ) );
//
// } // FUNCTION : AnalyticalMassProf_Einasto


//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Mass
// Description :  Calculate the enclosed mass of this cloud within radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Enclosed mass of this cloud within radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Mass( const double r )
{

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {
      if      ( r >= Table_r[params.Cloud_MassProfNBin-1] )   return Table_Enclosed_Mass[params.Cloud_MassProfNBin-1];
      else if ( r <= Table_r[0] )                             return Table_Enclosed_Mass[0];
      else                                                    return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Enclosed_Mass, r );
   }
   else
   {
      double enclosed_mass = NULL_REAL;
      double abs_error;

#     ifdef SUPPORT_GSL
      // arguments for the gsl integration
      const double lower_bound  = 0.0;
      const double upper_bound  = r;
      const double abs_err_lim  = 0;
      const double rel_err_lim  = 1e-7;
      const int    limit_size   = 1000;
      const int    integ_rule   = 1;    // 1 = GSL_INTEG_GAUSS15, the 15 point Gauss-Kronrod rule

      gsl_integration_workspace * w = gsl_integration_workspace_alloc( limit_size );

      gsl_function F;

      // integrand for the integration
      if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 ) F.function = &MassIntegrand_Plummer;
      else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 ) F.function = &MassIntegrand_NFW;
      else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 ) F.function = &MassIntegrand_Burkert;
      else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 ) F.function = &MassIntegrand_Jaffe;
      else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 ) F.function = &MassIntegrand_Hernquist;
      else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 ) F.function = &MassIntegrand_Einasto;

      // parameters for the integrand
      struct mass_integrand_params         integrand_params         = { params.Cloud_R0, params.Cloud_Rho0 };
      struct mass_integrand_params_Einasto integrand_params_Einasto = { params.Cloud_R0, params.Cloud_Rho0, params.Cloud_Einasto_Power_Factor };

      if      ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 ) F.params   = &integrand_params_Einasto;
      else                                                      F.params   = &integrand_params;

      // integration
      gsl_integration_qag( &F, lower_bound, upper_bound, abs_err_lim, rel_err_lim, limit_size, integ_rule, w, &enclosed_mass, &abs_error );

      gsl_integration_workspace_free( w );
#     endif // #ifdef SUPPORT_GSL

      return enclosed_mass;
   }

} // FUNCTION : Set_Mass



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Density
// Description :  Calculate the density of this cloud at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Density of this cloud at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Density( const double r )
{

   if ( strcmp( params.Cloud_Type, "Table" ) == 0 )
   {
      if ( r >= Table_r[params.Cloud_MassProfNBin-1] )
         return Table_Density[params.Cloud_MassProfNBin-1];
      else
         return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Density, r );
   }
   else
   {
      double rho;

      if      ( strcmp( params.Cloud_Type, "Plummer"   ) == 0 ) rho = AnalyticalDensProf_Plummer  ( r, params.Cloud_R0, params.Cloud_Rho0 );
      else if ( strcmp( params.Cloud_Type, "NFW"       ) == 0 ) rho = AnalyticalDensProf_NFW      ( r, params.Cloud_R0, params.Cloud_Rho0 );
      else if ( strcmp( params.Cloud_Type, "Burkert"   ) == 0 ) rho = AnalyticalDensProf_Burkert  ( r, params.Cloud_R0, params.Cloud_Rho0 );
      else if ( strcmp( params.Cloud_Type, "Jaffe"     ) == 0 ) rho = AnalyticalDensProf_Jaffe    ( r, params.Cloud_R0, params.Cloud_Rho0 );
      else if ( strcmp( params.Cloud_Type, "Hernquist" ) == 0 ) rho = AnalyticalDensProf_Hernquist( r, params.Cloud_R0, params.Cloud_Rho0 );
      else if ( strcmp( params.Cloud_Type, "Einasto"   ) == 0 ) rho = AnalyticalDensProf_Einasto  ( r, params.Cloud_R0, params.Cloud_Rho0, params.Cloud_Einasto_Power_Factor );

      return rho;

   }

} // FUNCTION : Set_Density



//-------------------------------------------------------------------------------------------------------
// Function    :  Set_Velocity
// Description :  Set the velocity of a particle at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Particle velocity
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Set_Velocity( const double x )
{

   double index, sum = 0;
   double psi_per = -potential(x);

   for (int k=0; k<params.Cloud_MassProfNBin; k++)
   {

      if ( psi[k] > psi_per )
      {
         index = k-1;
         break;
      }

      sum += prob_dens[k] *pow( psi_per-psi[k], 0.5 ) *delta;
   }

   double sum_rad, sum_mes=0, par, psi_ass;
   int index_ass = 0;

   sum_rad = Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

   sum_rad *= sum;

   for (int k=0; k<params.Cloud_MassProfNBin; k++)
   {

      if ( sum_mes > sum_rad )
      {
         index_ass = k-1;
         par = (sum_mes-sum_rad)/( prob_dens[index_ass] *pow( psi_per-psi[index_ass], 0.5 ) *delta );
         break;
      }

      sum_mes += prob_dens[k] *pow( psi_per-psi[k], 0.5 ) *delta;

      if ( k == params.Cloud_MassProfNBin-1 )
         index_ass = params.Cloud_MassProfNBin-1;
   }

   psi_ass = psi[index_ass] + delta*par;

   double kim = -2*(psi_ass-psi_per);

   if ( kim < 0.0 )
   {
      return 0;
   }

   double v = pow( kim, 0.5 );

   return v;

} // FUNCTION : Set_Velocity



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
// Solve Eddington's equation
double Par_EquilibriumIC::potential( const double x )
{
   const double r = x*params.Cloud_R0;

   if ( r >= Table_r[params.Cloud_MassProfNBin-1] )
   {
      return Table_Gravity_Potential[params.Cloud_MassProfNBin-1]*Table_r[params.Cloud_MassProfNBin-1]/r;
   }

   if ( r <= Table_r[0] )
   {
      return Table_Gravity_Potential[0];
   }

   return Mis_InterpolateFromTable( params.Cloud_MassProfNBin, Table_r, Table_Gravity_Potential, r );

} // FUNCTION : potential



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Integration_Eng_base( const double Eng, const int N_points )
{
   const double dEng = (Eng-Eng_min)/N_points;

   double Integration_result = 0;

   for (int i=0; i<N_points; i++)
   {
      const double Psi_L = Eng_min +        i*dEng;
      const double Psi_M = Eng_min +  (i+0.5)*dEng;
      const double Psi_R = Eng_min +    (i+1)*dEng;

      const int    index_Psi  = Mis_BinarySearch_Real( Table_Gravity_Potential, 0, params.Cloud_MassProfNBin-1,
                                                       -Psi_M ) + 1;

      if ( i == N_points-1 )   Integration_result += -2*Table_dRho_dx[index_Psi]*( sqrt( Eng-Psi_L ) );
      else                     Integration_result += -2*Table_dRho_dx[index_Psi]*( sqrt( Eng-Psi_L ) - sqrt( Eng-Psi_R ) );
   }

   return Integration_result;

} // FUNCTION : Integration_Eng_base



// Initialize physical parameter tables
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Mass
// Description :  Calculate the table of enclosed masses (vs. radius) of the cloud
//                by giving a known analytical density function of the cloud
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Mass()
{

   double dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   //Radius & Mass
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_r[b] = dr*b;
      Table_Enclosed_Mass[b] = Set_Mass( Table_r[b] );
   }

   //Rho
   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      Table_Density[b] = Set_Density( Table_r[b] );
   }

   Table_Density[0] = Table_Density[1];

   //Rhodr
   Table_dRho_dr[0] = ( Table_Density[1]-Table_Density[0] )/dr;
   for (int b=1; b<params.Cloud_MassProfNBin-1; b++)
   {
      int num=3;

      if      ( b == 0 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+1 );

      else if ( b == 1 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+2 );

      else if ( b == params.Cloud_MassProfNBin-2 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, params.Cloud_MassProfNBin-num/2-1, num/2+1 );

      else
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, b-num/2, num+1 );

   }

   Table_dRho_dr[params.Cloud_MassProfNBin-1] = Table_dRho_dr[params.Cloud_MassProfNBin-2];

} // FUNCTION : Init_Mass



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Pot
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by giving a known analytical density function of the cloud
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Pot()
{

   double dr = params.Cloud_MaxR / (params.Cloud_MassProfNBin-1);

   Table_Gravity_Field[0] =0;

   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      Table_Gravity_Field[b] = -NEWTON_G*Table_Enclosed_Mass[b]/pow( Table_r[b], 2 );
   }

   //Pot
   Table_Gravity_Potential[params.Cloud_MassProfNBin-1] = -NEWTON_G*Table_Enclosed_Mass[params.Cloud_MassProfNBin-1]/Table_r[params.Cloud_MassProfNBin-1];
   Eng_min = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];

   for (int b=params.Cloud_MassProfNBin-2; b>0; b--)
   {
      Table_Gravity_Potential[b] = Table_Gravity_Potential[b+1] + Table_Gravity_Field[b]*dr;
   }

   Table_Gravity_Potential[0] = Table_Gravity_Potential[1];

   //derho_overdx
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_dRho_dx[b] = -Table_dRho_dr[b]/Table_Gravity_Field[b];
   }

} // FUNCTION : Init_Pot



// Initialization through loading a file of table
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Mass_Table
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by loading a file of table
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Mass_Table()
{

   //Mass
   Table_Enclosed_Mass[0] = 0;
   double rho,dr,r;

   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      rho = (Table_Density[b] + Table_Density[b-1])/2;
      dr  = Table_r[b] - Table_r[b-1];
      r   = (Table_r[b] + Table_r[b-1])/2;

      Table_Enclosed_Mass[b] = Table_Enclosed_Mass[b-1] + 4*M_PI*pow(r,2) *rho * dr;
   }

   //Rhodr
   Table_dRho_dr[0] = (Table_Density[1]-Table_Density[0])/(Table_r[1]-Table_r[0]);

   for (int b=1; b<params.Cloud_MassProfNBin-1; b++)
   {
      int num = 3;
      if      ( b == 0 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+1 );

      else if ( b == 1 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, 0, num/2+2 );

      else if ( b == params.Cloud_MassProfNBin-2 )
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, params.Cloud_MassProfNBin-num/2-1, num/2+1 );
      else
         Table_dRho_dr[b] = Slope_LinearRegression( Table_r, Table_Density, b-num/2, num+1 );

   }

   Table_dRho_dr[params.Cloud_MassProfNBin-1] = Table_dRho_dr[params.Cloud_MassProfNBin-2];

} // FUNCTION : Init_Mass_Table



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Pot_Table
// Description :  Calculate the table of potential (vs. radius) of the cloud
//                by loading a file of table
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Pot_Table()
{

   Table_Gravity_Field[0] = 0;
   for (int b=1; b<params.Cloud_MassProfNBin; b++)
   {
      Table_Gravity_Field[b] = -NEWTON_G*Table_Enclosed_Mass[b]/pow( Table_r[b], 2 );
   }

   //Pot
   Table_Gravity_Potential[params.Cloud_MassProfNBin-1] = -NEWTON_G*Table_Enclosed_Mass[params.Cloud_MassProfNBin-1]/Table_r[params.Cloud_MassProfNBin-1];

   Eng_min = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];

   for (int b=params.Cloud_MassProfNBin-2; b>0; b--)
   {
      double dr = Table_r[b+1] - Table_r[b];

      Table_Gravity_Potential[b] = Table_Gravity_Potential[b+1] + Table_Gravity_Field[b]*dr;
   }

   Table_Gravity_Potential[0] = Table_Gravity_Potential[1];

   //derho_overdx
   for (int b=0; b<params.Cloud_MassProfNBin; b++)
   {
      Table_dRho_dx[b] = -Table_dRho_dr[b]/Table_Gravity_Field[b];
   }

} // FUNCTION : Init_Pot_Table



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Prob_Dens
// Description :  Calculate the probability density function of particles' velocities
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_Prob_Dens()
{

   double min, max;
   min   = -Table_Gravity_Potential[params.Cloud_MassProfNBin-1];
   max   = -Table_Gravity_Potential[1];
   delta = (max-min)/params.Cloud_MassProfNBin;

   double eng = min;

   for (int k =0; k<params.Cloud_MassProfNBin; k++)
   {

      psi[k] = eng;

      int_prob_dens[k] = Integration_Eng_base( eng, 1000 );

      eng += delta;

   }

   for (int k =0; k<params.Cloud_MassProfNBin; k++)
   {

      if      ( k == 0 )                           prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k,   5 );
      else if ( k == 1 )                           prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-1, 5 );
      else if ( k == params.Cloud_MassProfNBin-2 ) prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-3, 5 );
      else if ( k == params.Cloud_MassProfNBin-1 ) prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-4, 5 );
      else                                         prob_dens[k] = Slope_LinearRegression( psi, int_prob_dens, k-2, 5 );

      if ( prob_dens[k] < 0 )                      prob_dens[k] = 0;

   }

   SmoothArray( prob_dens, 0, params.Cloud_MassProfNBin );

} // FUNCTION : Init_Prob_Dens



//-------------------------------------------------------------------------------------------------------
// Function    :  Add_Ext_Pot
// Description :  Add exteranl potential through loading a file of table
//
// Note        :  Only will be activated if AddExtPot is turned on
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Add_Ext_Pot()
{

   if ( ! bool(params.AddExtPot) )  return;

   Aux_Message( stdout, "Loading External Potential Table...\n");

   int Tcol_r[1]={0};
   int Tcol_Pot[1]={1};

   double* Radius=NULL;
   double* Ext_Pot=NULL;

   int Row_r_Table;
   Row_r_Table = Aux_LoadTable( Radius, params.ExtPot_Table_Name, 1, Tcol_r, true, true );

   int Row_Ext_Pot_Table;
   Row_Ext_Pot_Table = Aux_LoadTable( Ext_Pot, params.ExtPot_Table_Name, 1, Tcol_Pot, true, true );

   Aux_Message( stdout, "Loading Ext_Pot Profile Table: %s\n", params.ExtPot_Table_Name );

   if ( Row_r_Table != Row_Ext_Pot_Table )
      Aux_Error( ERROR_INFO, "Ext_Pot row number is not equal to radius row number in the profile file !! Please check this file.\n" );

   if ( params.Cloud_MassProfNBin != Row_r_Table )
      Aux_Error( ERROR_INFO, "Cloud_MassProfNBin is not equal to the row number in profile file !!\n" );

   for (int i=0; i<params.Cloud_MassProfNBin; i++)
   {
      Table_Gravity_Potential[i] += Ext_Pot[i];
   }

} // FUNCTION : Add_Ext_Pot



// Auxiliary functions
//-------------------------------------------------------------------------------------------------------
// Function    :  GetParams
// Description :  Get the parameters from a file
//
// Note        :
//
// Parameter   :  filename   : file name of the input file
//                keyword    : the name of the parameter
//                para_num   : number of input parameters
//                containaer : container of return values
//
// Return      :  containaer
//
//-------------------------------------------------------------------------------------------------------
int Par_EquilibriumIC::GetParams( const char *filename, const char *keyword, const int para_num,
                                  const char *para_type, vector <string> &container )
{

   fstream file;
   file.open( filename, ios::in );

   string line;
   string para;

   if ( !file )
   {
      Aux_Error( ERROR_INFO, "Failed to open file : %s", filename );
   }
   else
   {
      do
      {
         getline( file, line );

         istringstream templine( line );

         string first_word;

         getline( templine, first_word, ' ' );

         if ( strcmp( first_word.c_str(), keyword ) == 0 )
         {
            for (int i=0; i<para_num; i++)
            {
               do
               {
                  getline( templine, para, ' ' );

               } while( para.length() == 0 );

               if ( strcmp( para_type, "string" ) == 0 )
                  container.push_back(para);
               else
                  return atoi( para.c_str() );
            }

            break;
         }

      } while( !file.eof() );
   }

   file.close();

   return 0;

} // FUNCTION : GetParams



//-------------------------------------------------------------------------------------------------------
// Function    :  RanVec_FixRadius
// Description :  Compute a random 3D vector with a fixed radius
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  r      : Input radius
//                RanVec : Array to store the random 3D vector
//
// Return      :  RanVec
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::RanVec_FixRadius( const double r, double RanVec[] )
{

   double Norm, RanR2;

   do
   {
      RanR2 = 0.0;

      for (int d=0; d<3; d++)
      {
         RanVec[d]  = Random_Num_Gen->GetValue( 0, -1.0, +1.0 );
         RanR2     += SQR( RanVec[d] );
      }

   } while ( RanR2 > 1.0 );

   Norm = r / sqrt( RanR2 );

   for (int d=0; d<3; d++)    RanVec[d] *= Norm;

} // FUNCTION : RanVec_FixRadius



//-------------------------------------------------------------------------------------------------------
// Function    : SmoothArray
// Description :
//
// Note        :
//
// Parameter   : array_x     : array of x data
//               index_start : the first index in the array for the smoothing
//               index_end   : the last index in the array for the smoothing
//
// Return      : array_x
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::SmoothArray( double* array_x, int index_start, int index_end )
{

   int    smoothing_n_elements = 10; // smoothing every "smoothing_n_elements" elements
   double smoothing_criterion  =  3; // smoothing when the ratio is larger than this criterion

   // set the elements as zero if its ratio to other elements is larger than smoothing_criterion
   for (int i=index_start; i<index_end-smoothing_n_elements+1; i++)
   {
      for (int j=i; j<i+smoothing_n_elements; j++)
      for (int k=i; k<i+smoothing_n_elements; k++)
      {
         if ( array_x[k] != 0  &&  fabs( array_x[j]/array_x[k] ) > smoothing_criterion )   array_x[j] = 0;
      }
   }

   // set those zero elements as the average of non-zero elements
   for (int i=index_start; i<index_end-smoothing_n_elements+1; i++)
   {
      double sum_of_nonzero = 0;
      int    num_of_nonzero = 0;
      double ave_of_nonzero = 0;

      // sum the non-zero elements
      for (int j=i; j<i+smoothing_n_elements; j++)
      {
         if ( array_x[j] != 0 )
         {
            sum_of_nonzero += array_x[j];
            num_of_nonzero ++;
         }
      }

      // average of non-zero elements
      if ( num_of_nonzero != 0 )   ave_of_nonzero = sum_of_nonzero/num_of_nonzero;

      // assign the average of non-zero elements to the zero element
      for (int j=i; j<i+smoothing_n_elements; j++)
      {
         if ( array_x[j] == 0 )
         {
            array_x[j] = ave_of_nonzero;
         }
      }
   }

} // FUNCTION : SmoothArray



//-------------------------------------------------------------------------------------------------------
// Function    : ArrayCovariance
// Description : Get the covariance between two arrays
//
// Note        : if x==y, then the covariance is the variance of x
//
// Parameter   : array_x     : array of x data
//               array_y     : array of y data
//               index_start : the first index in the array for the linear regression
//               n_elements  : number of elements for the linear regression
//
// Return      : covariance_xy
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::ArrayCovariance( const double* array_x, const double* array_y,
                                           const int index_start, const int n_elements )
{
   const double normalized_factor = 1.0/n_elements;

   // average
   double average_x = 0.0;
   double average_y = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)
   {
      average_x += array_x[i];
      average_y += array_y[i];
   }

   average_x *= normalized_factor;
   average_y *= normalized_factor;


   // covariance
   double covariance_xy = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)
   {
      covariance_xy += (array_x[i]-average_x)*(array_y[i]-average_y);
   }

   covariance_xy *= normalized_factor;


   return covariance_xy;

} // FUNCTION : ArrayCovariance



//-------------------------------------------------------------------------------------------------------
// Function    : Slope_LinearRegression
// Description : Get the slope of y-x using linear regression
//
// Note        :
//
// Parameter   : array_x     : array of x data
//               array_y     : array of y data
//               index_start : the first index in the array for the linear regression
//               n_elements  : number of elements for the linear regression
//
// Return      : slope_y
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::Slope_LinearRegression( const double* array_x, const double* array_y,
                                                  const int index_start, const int n_elements )
{

   const double variance_x    = ArrayCovariance( array_x, array_x, index_start, n_elements );
   const double covariance_xy = ArrayCovariance( array_x, array_y, index_start, n_elements );

   const double slope_y       = covariance_xy/variance_x;

   return slope_y;

} // FUNCTION : Slope_LinearRegression



#endif // #ifdef PARTICLE
