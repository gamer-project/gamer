#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

#include "Par_EquilibriumIC.h"



Par_EquilibriumIC::Par_EquilibriumIC( const char* Type )
{
   strcpy( Cloud_Type, Type );

   if      ( strcmp( Cloud_Type, "Table"     ) == 0 ) Cloud_Model = CLOUD_MODEL_TABLE;
   else if ( strcmp( Cloud_Type, "Plummer"   ) == 0 ) Cloud_Model = CLOUD_MODEL_PLUMMER;
   else if ( strcmp( Cloud_Type, "NFW"       ) == 0 ) Cloud_Model = CLOUD_MODEL_NFW;
   else if ( strcmp( Cloud_Type, "Burkert"   ) == 0 ) Cloud_Model = CLOUD_MODEL_BURKERT;
   else if ( strcmp( Cloud_Type, "Jaffe"     ) == 0 ) Cloud_Model = CLOUD_MODEL_JAFFE;
   else if ( strcmp( Cloud_Type, "Hernquist" ) == 0 ) Cloud_Model = CLOUD_MODEL_HERNQUIST;
   else if ( strcmp( Cloud_Type, "Einasto"   ) == 0 ) Cloud_Model = CLOUD_MODEL_EINASTO;
   else
      Aux_Error( ERROR_INFO, "Unsupported Cloud_Type \"%s\" for Par_EquilibriumIC !!\n", Cloud_Type );
}


Par_EquilibriumIC::~Par_EquilibriumIC()
{
   if ( Cloud_Model == CLOUD_MODEL_TABLE )
   {
      delete [] InputTable_DensProf_radius;
      delete [] InputTable_DensProf_density;
      delete [] InputTable_DensProf_enclosedmass;
   }

   delete Random_Num_Gen;

   delete [] RadArray_dr;
   delete [] RadArray_Radius;
   delete [] RadArray_Density;
   delete [] RadArray_EnclosedMass;
   delete [] RadArray_DensitySlope;
   delete [] RadArray_GraviField;
   delete [] RadArray_GraviPotential;
   delete [] RadArray_dRho_dPsi;

   delete [] EngArray_DistriFunc;
   delete [] EngArray_IntegDistriFunc;
   delete [] EngArray_BindingEnergy;

   if ( AddExtPot )
   {
      delete [] InputTable_ExtPot_radius;
      delete [] InputTable_ExtPot_potential;
      delete [] RadArray_ExternalPotential;
   }
}


void Par_EquilibriumIC::setCenter( const double Center_X, const double Center_Y, const double Center_Z )
{
   Cloud_Center[0] = Center_X;
   Cloud_Center[1] = Center_Y;
   Cloud_Center[2] = Center_Z;
}


void Par_EquilibriumIC::setBulkVel( const double BulkVel_X, const double BulkVel_Y, const double BulkVel_Z )
{
   Cloud_BulkVel[0] = BulkVel_X;
   Cloud_BulkVel[1] = BulkVel_Y;
   Cloud_BulkVel[2] = BulkVel_Z;
}


void Par_EquilibriumIC::setModelParameters( const double Rho0, const double R0 )
{
   Cloud_Rho0 = Rho0;
   Cloud_R0   = R0;
}


void Par_EquilibriumIC::setEinastoPowerFactor( const double EinastoPowerFactor )
{
   Cloud_Einasto_Power_Factor = EinastoPowerFactor;
}


void Par_EquilibriumIC::setDensProfTableFilename( const char* DensProfTableFilename )
{
   strcpy( DensProf_Table_Name, DensProfTableFilename );
}


void Par_EquilibriumIC::setParticleParameters( const long ParNum, const double MaxR, const int ArrayNBin, const int RSeed )
{
   Cloud_Par_Num      = ParNum;
   Cloud_MaxR         = MaxR;
   Cloud_Rad_NBin     = ArrayNBin;
   Cloud_RSeed        = RSeed;
}


long Par_EquilibriumIC::getParticleNumber( )
{
   return Cloud_Par_Num;
}


void Par_EquilibriumIC::setExternalPotential( const int AddingExternalPotential, const char* ExtPotTableFilename )
{
   AddExtPot = AddingExternalPotential;

   if ( AddExtPot )
      strcpy( ExtPot_Table_Name, ExtPotTableFilename );
}


void Par_EquilibriumIC::loadInputDensProfTable()
{
   if ( MPI_Rank == 0 ) Aux_Message( stdout, "Loading Density Profile Table: %s\n", DensProf_Table_Name );

   int NRowR;
   int NRowD;
   const int  Col_R[1] = {0};   // target column: radius
   const int  Col_D[1] = {1};   // target column: density

   NRowR = Aux_LoadTable( InputTable_DensProf_radius,  DensProf_Table_Name, 1, Col_R, true, true );
   NRowD = Aux_LoadTable( InputTable_DensProf_density, DensProf_Table_Name, 1, Col_D, true, true );

   // Check the number of rows are consistent
   if ( NRowR != NRowD )
      Aux_Error( ERROR_INFO, "The number of rows of density (%d) is not equal to the number of rows of radii (%d) in density table %s !!\n",
                             NRowD, NRowR, DensProf_Table_Name );
   else
      InputTable_DensProf_nbin = NRowR;

   // Check maximum radius in the density table must be larger than Cloud_MaxR
   if ( InputTable_DensProf_radius[InputTable_DensProf_nbin-1] < Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Maximum radius (%14.7e) in density table %s is smaller then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_DensProf_radius[InputTable_DensProf_nbin-1], DensProf_Table_Name, Cloud_MaxR );

   // Check minimum radius in the density table must be smaller than Cloud_MaxR
   if ( InputTable_DensProf_radius[0] > Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Minimum radius (%14.7e) in density table %s is larger then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_DensProf_radius[0], DensProf_Table_Name, Cloud_MaxR );

   // Set the default number of bins the same as the input table
   if ( Cloud_Rad_NBin < 0 )
      Cloud_Rad_NBin = Mis_BinarySearch_Real( InputTable_DensProf_radius, 0, InputTable_DensProf_nbin-1, Cloud_MaxR ) + 2;

   // Set the enclosed mass profile // TODO: this is not good, maybe remove it
   InputTable_DensProf_enclosedmass = new double [InputTable_DensProf_nbin];

   InputTable_DensProf_enclosedmass[0] = 0;
   for (int b=1; b<InputTable_DensProf_nbin; b++)
   {
      double dr      =      InputTable_DensProf_radius[b]  - InputTable_DensProf_radius[b-1];
      double r_mid   = 0.5*(InputTable_DensProf_radius[b]  + InputTable_DensProf_radius[b-1]);
      double rho_mid = 0.5*(InputTable_DensProf_density[b] + InputTable_DensProf_density[b-1]);

      InputTable_DensProf_enclosedmass[b] = InputTable_DensProf_enclosedmass[b-1] + 4*M_PI*SQR(r_mid)*rho_mid*dr;
   }
}


void Par_EquilibriumIC::loadInputExtPotTable()
{
   if ( MPI_Rank == 0 ) Aux_Message( stdout, "Loading ExtPot Profile Table: %s\n", ExtPot_Table_Name );

   int NRowR;
   int NRowP;
   const int  Col_R[1] = {0};   // target column: radius
   const int  Col_P[1] = {1};   // target column: potential

   NRowR = Aux_LoadTable( InputTable_ExtPot_radius,     ExtPot_Table_Name, 1, Col_R, true, true );
   NRowP = Aux_LoadTable( InputTable_ExtPot_potential,  ExtPot_Table_Name, 1, Col_P, true, true );

   // Check the number of rows are consistent
   if ( NRowR != NRowP )
      Aux_Error( ERROR_INFO, "The number of rows of potential (%d) is not equal to the number of rows of radii (%d) in ExtPot table %s !!\n",
                             NRowP, NRowR, ExtPot_Table_Name );
   else
      InputTable_ExtPot_nbin = NRowR;

   if ( Cloud_Rad_NBin != NRowR )
      Aux_Error( ERROR_INFO, "Cloud_Rad_NBin is not equal to the row number in profile file !!\n" );
}


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

   // Set random seeds
   Random_Num_Gen = new RandomNumber_t( 1 );
   Random_Num_Gen->SetSeed( 0, Cloud_RSeed );

   // Load the input density table
   if ( Cloud_Model == CLOUD_MODEL_TABLE )   loadInputDensProfTable();

   if ( Cloud_Rad_NBin < 2 )   Aux_Error( ERROR_INFO, "Cloud_Rad_NBin = %d is less than 2 !!\n", Cloud_Rad_NBin );

   // allocate memory
   RadArray_dr              = new double [Cloud_Rad_NBin];
   RadArray_Radius          = new double [Cloud_Rad_NBin];
   RadArray_Density         = new double [Cloud_Rad_NBin];
   RadArray_EnclosedMass    = new double [Cloud_Rad_NBin];
   RadArray_DensitySlope    = new double [Cloud_Rad_NBin];
   RadArray_GraviField      = new double [Cloud_Rad_NBin];
   RadArray_GraviPotential  = new double [Cloud_Rad_NBin];
   RadArray_dRho_dPsi       = new double [Cloud_Rad_NBin];

   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      Cloud_Eng_NBin = Cloud_Rad_NBin-1;
   else
      Cloud_Eng_NBin = Cloud_Rad_NBin;

   EngArray_DistriFunc      = new double [Cloud_Eng_NBin];
   EngArray_IntegDistriFunc = new double [Cloud_Eng_NBin];
   EngArray_BindingEnergy   = new double [Cloud_Eng_NBin];

   Rad_LastIdx = Cloud_Rad_NBin-1;
   Eng_LastIdx = Cloud_Eng_NBin-1;

   setRadArray_Radius();
   setRadArray_Density();
   setRadArray_EnclosedMass();
   setRadArray_DensitySlope();
   setRadArray_GraviField();
   setRadArray_GraviPotential();
   setRadArray_dRho_dPsi();

   if ( AddExtPot )
   {
      loadInputExtPotTable();
      RadArray_ExternalPotential = new double [Cloud_Rad_NBin];
      setRadArray_ExternalPotential();
      addExternalPotential();
   }

   printf( "Radius: [" );
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_Radius[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_Radius[b] );
   printf( "]\n");

   printf( "Dens:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_Density[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_Density[b] );
   printf( "]\n");

   printf( "Mass:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_EnclosedMass[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_EnclosedMass[b] );
   printf( "]\n");

   printf( "dDdr:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_DensitySlope[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_DensitySlope[b] );
   printf( "]\n");

   printf( "Pote:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_GraviPotential[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_GraviPotential[b] );
   printf( "]\n");

   printf( "GreF:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_GraviField[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_GraviField[b] );
   printf( "]\n");

   printf( "dDdx:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", RadArray_dRho_dPsi[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin-10; b<Cloud_Rad_NBin; b++)  printf(" %21.14e,", RadArray_dRho_dPsi[b] );
   printf( "]\n");

   Init_EngArray();

   printf( "BEng:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", EngArray_BindingEnergy[b] );
   printf( " ... ");
   for (int b=Cloud_Eng_NBin-10; b<Cloud_Eng_NBin; b++)  printf(" %21.14e,", EngArray_BindingEnergy[b] );
   printf( "]\n");

   printf( "InDF:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", EngArray_IntegDistriFunc[b] );
   printf( " ... ");
   for (int b=Cloud_Eng_NBin-10; b<Cloud_Eng_NBin; b++)  printf(" %21.14e,", EngArray_IntegDistriFunc[b] );
   printf( "]\n");

   printf( "DisF:   [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", EngArray_DistriFunc[b] );
   printf( " ... ");
   for (int b=Cloud_Eng_NBin-10; b<Cloud_Eng_NBin; b++)  printf(" %21.14e,", EngArray_DistriFunc[b] );
   printf( "]\n");

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

   int Cloud_Rad_NBin_2;

   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      Cloud_Rad_NBin_2 = Cloud_Rad_NBin-1;
   else
      Cloud_Rad_NBin_2 = Cloud_Rad_NBin;

   // determine the total enclosed mass within the maximum radius
   const double TotCloudMass = getEnclosedMass( Cloud_MaxR );
   const double ParticleMass = TotCloudMass/Cloud_Par_Num;

   // construct the mass profile table
   Table_MassProf_r = new double [Cloud_Rad_NBin_2];
   Table_MassProf_M = new double [Cloud_Rad_NBin_2];

   const double dr = Cloud_MaxR / (Cloud_Rad_NBin_2-1);

   for (int b=0; b<Cloud_Rad_NBin; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = getEnclosedMass( Table_MassProf_r[b] );
   }

   printf( "T_r:    [" );
   for (int b=0; b<10; b++)  printf(" %21.14e,", Table_MassProf_r[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin_2-10; b<Cloud_Rad_NBin_2; b++)  printf(" %21.14e,", Table_MassProf_r[b] );
   printf( "]\n");

   printf( "T_M:    [");
   for (int b=0; b<10; b++)  printf(" %21.14e,", Table_MassProf_M[b] );
   printf( " ... ");
   for (int b=Cloud_Rad_NBin_2-10; b<Cloud_Rad_NBin_2; b++)  printf(" %21.14e,", Table_MassProf_M[b] );
   printf( "]\n");


   double  RandomSampleM;
   double  RandomSampleR;
   double  RandomSampleV;
   double  ErrM_Max =-1.0;
   double  RandomVectorR[3];
   double  RandomVectorV[3];

   // set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)
   {
      // mass
      Mass_AllRank[p] = ParticleMass;

      // position
      // --> sample from the cumulative mass profile with linear interpolation
      RandomSampleM = TotCloudMass*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

      if ( Cloud_Model == CLOUD_MODEL_TABLE )
         RandomSampleR = Mis_InterpolateFromTable( Cloud_Rad_NBin_2, Table_MassProf_M, Table_MassProf_r, RandomSampleM );
      else
         RandomSampleR = Mis_InterpolateFromTable( Cloud_Rad_NBin, RadArray_EnclosedMass, RadArray_Radius, RandomSampleM );

      // record the maximum error
      ErrM_Max = fmax( fabs( ( getEnclosedMass( RandomSampleR ) - RandomSampleM )/RandomSampleM ), ErrM_Max );

      // randomly set the position vector with a given radius
      RanVec_FixRadius( RandomSampleR, RandomVectorR );
      for (int d=0; d<3; d++)   Pos_AllRank[d][p] = Cloud_Center[d]  + RandomVectorR[d];

      // check periodicity
      for (int d=0; d<3; d++)
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )   Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );

      // velocity
      RandomSampleV = getRandomSampleVelocity( RandomSampleR );

      // randomly set the velocity vector with the given amplitude
      RanVec_FixRadius( RandomSampleV, RandomVectorV );
      for (int d=0; d<3; d++)   Vel_AllRank[d][p] = Cloud_BulkVel[d] + RandomVectorV[d];

   } // for (long p=0; p<NPar_AllRank; p++)

   Aux_Message( stdout, "   Total enclosed mass within MaxR  = %13.7e\n",  TotCloudMass );
   Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParticleMass );
   Aux_Message( stdout, "   Maximum mass interpolation error = %13.7e\n",  ErrM_Max     );

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


struct mass_integrand_params_Table
{
   int     NBin;
   double* Table_R;
   double* Table_D;
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
// Function    :
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
double ExtendedInterpolatedTable( const double x, const int N, const double Table_x[], const double Table_y[] )
{
   if      ( x <= Table_x[0]   )   return Table_y[0];
   else if ( x >= Table_x[N-1] )   return Table_y[N-1];
   else                            return Mis_InterpolateFromTable( N, Table_x, Table_y, x );
}


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
double MassIntegrand_Table( const double r, void* parameters )
{
   if ( r == 0.0 ) return 0.0;

   struct  mass_integrand_params_Table *p = (struct mass_integrand_params_Table *) parameters;
   int     NBin                       = p->NBin;
   double* Table_R                    = p->Table_R;
   double* Table_D                    = p->Table_D;

   return 4*M_PI*SQR(r)*ExtendedInterpolatedTable( r, NBin, Table_R, Table_D );

} // FUNCTION : MassIntegrand_Table


//-------------------------------------------------------------------------------------------------------
// Function    :  getDensity
// Description :  Calculate the density of this cloud at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Density of this cloud at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getDensity( const double r )
{
   double rho;

   if      ( Cloud_Model == CLOUD_MODEL_TABLE     ) rho = ExtendedInterpolatedTable   ( r, InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_density );
   else if ( Cloud_Model == CLOUD_MODEL_PLUMMER   ) rho = AnalyticalDensProf_Plummer  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_NFW       ) rho = AnalyticalDensProf_NFW      ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_BURKERT   ) rho = AnalyticalDensProf_Burkert  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_JAFFE     ) rho = AnalyticalDensProf_Jaffe    ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST ) rho = AnalyticalDensProf_Hernquist( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_EINASTO   ) rho = AnalyticalDensProf_Einasto  ( r, Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor );
   else
      Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );

   return rho;

} // FUNCTION : getDensity


//-------------------------------------------------------------------------------------------------------
// Function    :  getEnclosedMass
// Description :  Calculate the enclosed mass of this cloud within radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Enclosed mass of this cloud within radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getEnclosedMass( const double r )
{

   double enclosed_mass = NULL_REAL;

   if ( Cloud_Model == CLOUD_MODEL_TABLE )
   {
      enclosed_mass = ExtendedInterpolatedTable( r, InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_enclosedmass );
   }
   else
   {
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
      if      ( Cloud_Model == CLOUD_MODEL_PLUMMER   ) F.function = &MassIntegrand_Plummer;
      else if ( Cloud_Model == CLOUD_MODEL_NFW       ) F.function = &MassIntegrand_NFW;
      else if ( Cloud_Model == CLOUD_MODEL_BURKERT   ) F.function = &MassIntegrand_Burkert;
      else if ( Cloud_Model == CLOUD_MODEL_JAFFE     ) F.function = &MassIntegrand_Jaffe;
      else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST ) F.function = &MassIntegrand_Hernquist;
      else if ( Cloud_Model == CLOUD_MODEL_EINASTO   ) F.function = &MassIntegrand_Einasto;
      else if ( Cloud_Model == CLOUD_MODEL_TABLE     ) F.function = &MassIntegrand_Table;
      else
         Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );

      // parameters for the integrand
      struct mass_integrand_params         integrand_params         = { Cloud_R0, Cloud_Rho0 };
      struct mass_integrand_params_Einasto integrand_params_Einasto = { Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor };
      struct mass_integrand_params_Table   integrand_params_Table   = { InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_density };

      if      ( Cloud_Model == CLOUD_MODEL_EINASTO   ) F.params     = &integrand_params_Einasto;
      else if ( Cloud_Model == CLOUD_MODEL_TABLE     ) F.params     = &integrand_params_Table;
      else                                             F.params     = &integrand_params;

      // integration
      gsl_integration_qag( &F, lower_bound, upper_bound, abs_err_lim, rel_err_lim, limit_size, integ_rule, w, &enclosed_mass, &abs_error );

      gsl_integration_workspace_free( w );
#     endif // #ifdef SUPPORT_GSL
   }

   return enclosed_mass;

} // FUNCTION : getEnclosedMass


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
double Par_EquilibriumIC::getGraviPotential( const double r )
{
   // Note this direction is differt: from table
   if      ( r >= RadArray_Radius[Rad_LastIdx] )   return RadArray_GraviPotential[Rad_LastIdx]*RadArray_Radius[Rad_LastIdx]/r;
   else if ( r <= RadArray_Radius[0] )             return RadArray_GraviPotential[0];
   else                                            return Mis_InterpolateFromTable( Cloud_Rad_NBin, RadArray_Radius, RadArray_GraviPotential, r );

} // FUNCTION : getGraviPotential


//-------------------------------------------------------------------------------------------------------
// Function    :  getExternalPotential
// Description :  Calculate the external potential at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  external potential at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getExternalPotential( const double r )
{
   double extpot;

   extpot = ExtendedInterpolatedTable   ( r, InputTable_ExtPot_nbin, InputTable_ExtPot_radius, InputTable_ExtPot_potential );

   return extpot;

} // FUNCTION : getExternalPotential



//-------------------------------------------------------------------------------------------------------
// Function    :  getRandomSampleVelocity
// Description :  Set the velocity of a particle at radius r
//
// Note        :
//
// Parameter   :  r : radius
//
// Return      :  Particle velocity
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getRandomSampleVelocity( const double r )
{

   // TODO: use a table of probability and use interpolation
   double TotalProbability = 0;
   const double Psi = -getGraviPotential(r);

   for (int k=0; k<Cloud_Eng_NBin; k++)
   {
      if ( EngArray_BindingEnergy[k] > Psi )   break;

      TotalProbability += EngArray_DistriFunc[k]*sqrt( Psi-EngArray_BindingEnergy[k] )*EngArray_dBindingEnergy; //TODO: factor 2 in sqrt
   }

   double sum_mes        =  0;
   double Fraction_dEng  =  0;
   int RandomSampleIndex = -1;

   const double RandomSampleProbability = TotalProbability*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

   for (int k=0; k<Cloud_Eng_NBin; k++)
   {

      if ( sum_mes > RandomSampleProbability )
      {
         RandomSampleIndex = k-1;
         Fraction_dEng = (sum_mes-RandomSampleProbability)/( EngArray_DistriFunc[RandomSampleIndex]*sqrt( Psi-EngArray_BindingEnergy[RandomSampleIndex] )*EngArray_dBindingEnergy ); // TODO: (1 - ...)
         break;
      }

      sum_mes += EngArray_DistriFunc[k]*sqrt( Psi-EngArray_BindingEnergy[k] )*EngArray_dBindingEnergy;
   }

   if ( RandomSampleIndex < 0 )   RandomSampleIndex = Eng_LastIdx;

   const double RandomSampleBindingEnergy = EngArray_BindingEnergy[RandomSampleIndex] + EngArray_dBindingEnergy*Fraction_dEng;

   const double RandomSampleKineticEnergy = 2*(Psi-RandomSampleBindingEnergy);

   const double RandomSampleVelocity      = ( RandomSampleKineticEnergy < 0.0 ) ? 0 : sqrt( RandomSampleKineticEnergy );

   return RandomSampleVelocity;

} // FUNCTION : getRandomSampleVelocity



//-------------------------------------------------------------------------------------------------------
// Function    :
// Description :
//
// Note        :   \int_{0}^{\mathcal{E}} \frac{ 1 }{ \sqrt{ \mathcal{E} -\Psi } } \frac{ d\rho }{ d\Psi } d\Psi
//               = \int_{0}^{\mathcal{E}} -2 \frac{ d\rho }{ d\Psi } d\sqrt{ \mathcal{E} -\Psi } }
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getIntegratedDistributionFunction( const double Psi_Min, const double Psi_Max, const int N_points )
{
   const double dPsi = ( Psi_Max - Psi_Min )/N_points;

   double integral = 0;

   for (int i=0; i<N_points; i++)
   {
      const double Psi_L = Psi_Min +       i*dPsi; // TODO: use an array psi
      const double Psi_M = Psi_Min + (i+0.5)*dPsi;
      const double Psi_R = Psi_Min +   (i+1)*dPsi;

      const int    index_Psi  = Mis_BinarySearch_Real( RadArray_GraviPotential, 0, Rad_LastIdx, -Psi_M ) + 1;

      // TODO: seperate dsqrt(Psi_Max-Psi)
      if ( i == N_points-1 )   integral += -2*RadArray_dRho_dPsi[index_Psi]*(                       - sqrt( Psi_Max-Psi_L ) );
      else                     integral += -2*RadArray_dRho_dPsi[index_Psi]*( sqrt( Psi_Max-Psi_R ) - sqrt( Psi_Max-Psi_L ) );
   }

   return integral;

} // FUNCTION : getIntegratedDistributionFunction



// Initialize physical parameter tables
//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_Radius
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_Radius()
{
   RadArray_Radius[0] = 0;

   const double delta_r_linear = Cloud_MaxR/(Cloud_Rad_NBin-1);
   for (int b=1; b<Cloud_Rad_NBin; b++)   RadArray_Radius[b] = delta_r_linear*b;

   //const double ratio_r_log    = POW( Cloud_MaxR, 1.0/(Cloud_Rad_NBin-1) );
   //for (int b=1; b<Cloud_Rad_NBin; b++)   RadArray_Radius[b] = POW( ratio_r_log, b );

   for (int b=0; b<Rad_LastIdx; b++)   RadArray_dr[b] = RadArray_Radius[b+1] - RadArray_Radius[b];

   RadArray_dr[Rad_LastIdx] = 0;

} // FUNCTION : setRadArray_Radius


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_Density
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_Density()
{
   for (int b=1; b<Cloud_Rad_NBin; b++)   RadArray_Density[b] = getDensity( RadArray_Radius[b] );

   // when r=0
   RadArray_Density[0] = RadArray_Density[1];

} // FUNCTION : setRadArray_Density


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_EnclosedMass
// Description :  Calculate the table of enclosed masses (vs. radius) of the cloud
//                by giving a known analytical density function of the cloud
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_EnclosedMass()
{

   for (int b=1; b<Cloud_Rad_NBin; b++)   RadArray_EnclosedMass[b] = getEnclosedMass( RadArray_Radius[b] );

   // when r=0
   RadArray_EnclosedMass[0] = 0;

} // FUNCTION : setRadArray_EnclosedMass


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_DensitySlope
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_DensitySlope()
{
   const int Npoints = 3; // number of points for the linear regression

   RadArray_DensitySlope[0] =                 (RadArray_Density[1] - RadArray_Density[0])/RadArray_dr[0];


   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      RadArray_DensitySlope[1] =              (RadArray_Density[2] - RadArray_Density[1])/RadArray_dr[1];
   else
      RadArray_DensitySlope[1] =              Slope_LinearRegression( RadArray_Radius, RadArray_Density,                           0, Npoints/2+2 );

   for (int b=2; b<Cloud_Rad_NBin-2; b++)
      RadArray_DensitySlope[b] =              Slope_LinearRegression( RadArray_Radius, RadArray_Density,                 b-Npoints/2,     Npoints );

   RadArray_DensitySlope[Cloud_Rad_NBin-2]  = Slope_LinearRegression( RadArray_Radius, RadArray_Density, Cloud_Rad_NBin-Npoints/2-1, Npoints/2+1 );

   RadArray_DensitySlope[Rad_LastIdx]       = RadArray_DensitySlope[Cloud_Rad_NBin-2];

} // FUNCTION : setRadArray_DensitySlope


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_GraviField
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_GraviField()
{
   RadArray_GraviField[0] = 0;

   for (int b=1; b<Cloud_Rad_NBin; b++)
      RadArray_GraviField[b] = -NEWTON_G*RadArray_EnclosedMass[b]/SQR( RadArray_Radius[b] );

} // FUNCTION : setRadArray_GraviField


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_GraviPotential
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_GraviPotential()
{
   RadArray_GraviPotential[Rad_LastIdx] = -NEWTON_G*RadArray_EnclosedMass[Rad_LastIdx]/RadArray_Radius[Rad_LastIdx];

   for (int b=Cloud_Rad_NBin-2; b>1; b--)
      RadArray_GraviPotential[b] = RadArray_GraviPotential[b+1] + RadArray_GraviField[b]*RadArray_dr[b];

   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      RadArray_GraviPotential[1] = RadArray_GraviPotential[2];
   else
      RadArray_GraviPotential[1] = RadArray_GraviPotential[2] + RadArray_GraviField[1]*RadArray_dr[1];

   RadArray_GraviPotential[0] = RadArray_GraviPotential[1];

} // FUNCTION : setRadArray_GraviPotential


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_ExternalPotential
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_ExternalPotential()
{
   for (int b=0; b<Cloud_Rad_NBin; b++)   RadArray_ExternalPotential[b] = getExternalPotential( RadArray_Radius[b] );

} // FUNCTION : setRadArray_GraviPotential


//-------------------------------------------------------------------------------------------------------
// Function    :  setRadArray_dRho_dPsi
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setRadArray_dRho_dPsi()
{
   for (int b=0; b<Cloud_Rad_NBin; b++)   RadArray_dRho_dPsi[b] = RadArray_DensitySlope[b]/RadArray_GraviField[b];

} // FUNCTION : setRadArray_dRho_dPsi


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_EngArray
// Description :  Calculate the probability density function of particles' velocities
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::Init_EngArray()
{
   EngArray_MinBindingEnergy = -RadArray_GraviPotential[Rad_LastIdx]; // TODO: why not zero
   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      EngArray_MaxBindingEnergy = -RadArray_GraviPotential[2];
   else
      EngArray_MaxBindingEnergy = -RadArray_GraviPotential[1];
   EngArray_dBindingEnergy   = (EngArray_MaxBindingEnergy-EngArray_MinBindingEnergy)/Cloud_Eng_NBin;

   // Set the binding energy
   EngArray_BindingEnergy[0] = EngArray_MinBindingEnergy;
   for (int k=1; k<Cloud_Eng_NBin; k++)   EngArray_BindingEnergy[k] = EngArray_BindingEnergy[k-1] + EngArray_dBindingEnergy;

   // Set the integrated distribution function
   for (int k=0; k<Cloud_Eng_NBin; k++)   EngArray_IntegDistriFunc[k] = getIntegratedDistributionFunction( EngArray_MinBindingEnergy, EngArray_BindingEnergy[k], 1000 ); // TODO: why not from zero

   // Set the distribution function
   for (int k =0; k<Cloud_Eng_NBin; k++)
   {
      if      ( k <= 1                )   EngArray_DistriFunc[k] = Slope_LinearRegression( EngArray_BindingEnergy, EngArray_IntegDistriFunc,                0, 5 );
      else if ( k >= Cloud_Eng_NBin-2 )   EngArray_DistriFunc[k] = Slope_LinearRegression( EngArray_BindingEnergy, EngArray_IntegDistriFunc, Cloud_Eng_NBin-5, 5 );
      else                                EngArray_DistriFunc[k] = Slope_LinearRegression( EngArray_BindingEnergy, EngArray_IntegDistriFunc,              k-2, 5 );

      // check negative distribution function
      if ( EngArray_DistriFunc[k] < 0 )   EngArray_DistriFunc[k] = 0;
   }

   // Smooth the distribution function
   SmoothArray( EngArray_DistriFunc, 0, Cloud_Eng_NBin );

} // FUNCTION : Init_EngArray


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
void Par_EquilibriumIC::addExternalPotential()
{
   if ( AddExtPot == 0 )  return;

   for (int b=0; b<Cloud_Rad_NBin; b++)   RadArray_GraviPotential[b] += RadArray_ExternalPotential[b];

} // FUNCTION : addExternalPotential



// Auxiliary functions
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
