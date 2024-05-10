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

   if ( AddExtPot_Table )
   {
      delete [] InputTable_ExtPot_radius;
      delete [] InputTable_ExtPot_potential;
   }

   delete Random_Num_Gen;

   delete [] RArray_R;
   delete [] RArray_Rho;
   delete [] RArray_M_Enc;
   delete [] RArray_dRho_dR;
   delete [] RArray_G;
   delete [] RArray_Phi;
   delete [] RArray_dRho_dPsi;

   delete [] EArray_DFunc;
   delete [] EArray_IntDFunc;
   delete [] EArray_E;
}


void Par_EquilibriumIC::setCenterAndBulkVel( const double Center_X, const double Center_Y, const double Center_Z,
                                             const double BulkVel_X, const double BulkVel_Y, const double BulkVel_Z )
{
   Cloud_Center[0]  = Center_X;
   Cloud_Center[1]  = Center_Y;
   Cloud_Center[2]  = Center_Z;

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


void Par_EquilibriumIC::setParticleParameters( const long ParNum, const double MaxR, const int Radial_NBin, const int Energy_NBin, const int RSeed )
{
   Cloud_Par_Num  = ParNum;
   Cloud_MaxR     = MaxR;
   RNBin          = Radial_NBin;
   ENBin          = Energy_NBin;
   Cloud_RSeed    = RSeed;
}


double Par_EquilibriumIC::getTotCloudMass()
{
   return TotCloudMass;
}


double Par_EquilibriumIC::getParticleMass()
{
   return ParticleMass;
}


double Par_EquilibriumIC::getMaxMassError()
{
   return MaxMassError;
}


void Par_EquilibriumIC::setExternalPotential( const int AddingExternalPotential_Analytical, const int AddingExternalPotential_Table, const char* ExtPotTableFilename )
{
   AddExtPot_Analytical = AddingExternalPotential_Analytical;
   AddExtPot_Table      = AddingExternalPotential_Table;

   if ( AddExtPot_Analytical  &&  AddExtPot_Table )
      Aux_Error( ERROR_INFO, "AddExtPot_Analytical and AddExtPot_Table in Par_EquilibriumIC cannot both be turned on !!\n" );

   if ( AddExtPot_Table )   strcpy( ExtPot_Table_Name, ExtPotTableFilename );
}


void Par_EquilibriumIC::loadInputDensProfTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ...\n", DensProf_Table_Name );

   const int  Col_R[1] = {0};   // target column: radius
   const int  Col_D[1] = {1};   // target column: density

   const int NRowR = Aux_LoadTable( InputTable_DensProf_radius,  DensProf_Table_Name, 1, Col_R, true, true );
   const int NRowD = Aux_LoadTable( InputTable_DensProf_density, DensProf_Table_Name, 1, Col_D, true, true );

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
   if ( RNBin < 0 )
      RNBin = Mis_BinarySearch_Real( InputTable_DensProf_radius, 0, InputTable_DensProf_nbin-1, Cloud_MaxR ) + 2;

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

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ... done\n", DensProf_Table_Name );
}


void Par_EquilibriumIC::loadInputExtPotTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading ExtPot Profile Table: %s ...\n", ExtPot_Table_Name );

   const int  Col_R[1] = {0};   // target column: radius
   const int  Col_P[1] = {1};   // target column: potential

   const int NRowR = Aux_LoadTable( InputTable_ExtPot_radius,     ExtPot_Table_Name, 1, Col_R, true, true );
   const int NRowP = Aux_LoadTable( InputTable_ExtPot_potential,  ExtPot_Table_Name, 1, Col_P, true, true );

   // Check the number of rows are consistent
   if ( NRowR != NRowP )
      Aux_Error( ERROR_INFO, "The number of rows of potential (%d) is not equal to the number of rows of radii (%d) in ExtPot table %s !!\n",
                             NRowP, NRowR, ExtPot_Table_Name );
   else
      InputTable_ExtPot_nbin = NRowR;

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading ExtPot Profile Table: \"%s\" ... done\n", ExtPot_Table_Name );
}


//-------------------------------------------------------------------------------------------------------
// Function    :  initialize
// Description :  Initialize all necessary tables of physical parameters, including radius, mass, density, gravitational potential
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::initialize()
{

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "Must enable SUPPORT_GSL for Par_EquilibriumIC !!\n" );
#  endif

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Initializing Par_EquilibriumIC ...\n" );

   // Set random seeds
   Random_Num_Gen = new RandomNumber_t( 1 );
   Random_Num_Gen->SetSeed( 0, Cloud_RSeed );

   // Load the input density table
   if ( Cloud_Model == CLOUD_MODEL_TABLE )   loadInputDensProfTable();
   if ( AddExtPot_Table )                    loadInputExtPotTable();

   if ( RNBin < 2 )   Aux_Error( ERROR_INFO, "RNBin = %d is less than 2 !!\n", RNBin );

   // allocate memory
   RArray_R         = new double [RNBin];
   RArray_Rho       = new double [RNBin];
   RArray_M_Enc     = new double [RNBin];
   RArray_dRho_dR   = new double [RNBin];
   RArray_G         = new double [RNBin];
   RArray_Phi       = new double [RNBin];
   RArray_dRho_dPsi = new double [RNBin];

   RLastIdx = RNBin-1;
   constructRadialArray();

   EArray_DFunc     = new double [ENBin];
   EArray_IntDFunc  = new double [ENBin];
   EArray_E         = new double [ENBin];

   ELastIdx = ENBin-1;
   constructEnergyArray();

   //----------------------------------------------------------------------------------------------
   printf( "Radius: [" );
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_R[b] );
   printf( "]\n");

   printf( "Dens:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_Rho[b] );
   printf( "]\n");

   printf( "Mass:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_M_Enc[b] );
   printf( "]\n");

   printf( "dDdr:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_dRho_dR[b] );
   printf( "]\n");

   printf( "Pote:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_Phi[b] );
   printf( "]\n");

   printf( "GreF:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_G[b] );
   printf( "]\n");

   printf( "dDdx:   [");
   for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_dRho_dPsi[b] );
   printf( "]\n");

   printf( "BEng:   [");
   for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_E[b] );
   printf( "]\n");

   printf( "InDF:   [");
   for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_IntDFunc[b] );
   printf( "]\n");

   printf( "DisF:   [");
   for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_DFunc[b] );
   printf( "]\n");
   //----------------------------------------------------------------------------------------------

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Initializing Par_EquilibriumIC ... done\n" );

} // FUNCTION : initialize


//-------------------------------------------------------------------------------------------------------
// Function    :  constructParticles
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
void Par_EquilibriumIC::constructParticles( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx0 )
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Constructing Par_EquilibriumIC ...\n" );

   // TODO: remove them
   //------------------------------------------------------------------------------------------
   int RNBin_2;
   if ( Cloud_Model == CLOUD_MODEL_TABLE )   RNBin_2 = RNBin-1;
   else                                      RNBin_2 = RNBin;

   // construct the mass profile table
   double *Table_MassProf_r = new double [RNBin_2];
   double *Table_MassProf_M = new double [RNBin_2];

   const double dr = Cloud_MaxR / (RNBin_2-1);

   for (int b=0; b<RNBin; b++)
   {
      Table_MassProf_r[b] = dr*b;
      Table_MassProf_M[b] = getEnclosedMass( Table_MassProf_r[b] );
   }

   printf( "T_r:    [" );
   for (int b=0; b<RNBin_2; b++)  printf(" %21.14e,", Table_MassProf_r[b] );
   printf( "]\n");

   printf( "T_M:    [");
   for (int b=0; b<RNBin_2; b++)  printf(" %21.14e,", Table_MassProf_M[b] );
   printf( "]\n");
   //------------------------------------------------------------------------------------------


   // determine the total enclosed mass within the maximum radius
   TotCloudMass = getEnclosedMass( Cloud_MaxR );
   ParticleMass = TotCloudMass/Cloud_Par_Num;

   double  RandomSampleM;
   double  RandomSampleR;
   double  RandomSampleV;
   double  RandomVectorR[3];
   double  RandomVectorV[3];

   // set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)
   {
      // position, sample from the cumulative mass profile with linear interpolation
      RandomSampleM = TotCloudMass*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

      // TODO: combine them
      if ( Cloud_Model == CLOUD_MODEL_TABLE )
         RandomSampleR = Mis_InterpolateFromTable( RNBin_2, Table_MassProf_M, Table_MassProf_r, RandomSampleM );
      else
         RandomSampleR = Mis_InterpolateFromTable( RNBin, RArray_M_Enc, RArray_R, RandomSampleM );

      // randomly set the position vector with a given radius
      RandomVector_GivenLength( RandomSampleR, RandomVectorR );

      // velocity
      RandomSampleV = getRandomSampleVelocity( RandomSampleR );

      // randomly set the velocity vector with the given amplitude
      RandomVector_GivenLength( RandomSampleV, RandomVectorV );

      Mass_AllRank[p] = ParticleMass;
      for (int d=0; d<3; d++)   Pos_AllRank[d][p] = Cloud_Center[d]  + RandomVectorR[d];
      for (int d=0; d<3; d++)   Vel_AllRank[d][p] = Cloud_BulkVel[d] + RandomVectorV[d];

      // check periodicity
      for (int d=0; d<3; d++)
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )   Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );

      // record the maximum error
      MaxMassError = fmax( fabs( ( getEnclosedMass( RandomSampleR ) - RandomSampleM )/RandomSampleM ), MaxMassError );

   } // for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)

   // free memory
   delete [] Table_MassProf_r;
   delete [] Table_MassProf_M;

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing Par_EquilibriumIC ... done\n" );

} // FUNCTION : constructParticle


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
// Function    :  AnalyticalExternalPotential
// Description :  Analytical external potential
//
// Note        :
//
// Parameter   :  r                    : input radius
//
// Return      :  external potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalExternalPotential( const double r )
{
   return 0.0;

} // FUNCTION : AnalyticalExternalPotential



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
   if      ( Cloud_Model == CLOUD_MODEL_TABLE     ) return ExtendedInterpolatedTable   ( r, InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_density );
   else if ( Cloud_Model == CLOUD_MODEL_PLUMMER   ) return AnalyticalDensProf_Plummer  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_NFW       ) return AnalyticalDensProf_NFW      ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_BURKERT   ) return AnalyticalDensProf_Burkert  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_JAFFE     ) return AnalyticalDensProf_Jaffe    ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST ) return AnalyticalDensProf_Hernquist( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_EINASTO   ) return AnalyticalDensProf_Einasto  ( r, Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor );
   else
   {
      Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );
      return 0.0;
   }

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
      const double lower_limit  = 0.0;
      const double upper_limit  = r;
      const double abs_err_lim  = 0;
      const double rel_err_lim  = 1e-7;
      const int    integ_size   = 1000; // TODO: more efficient
      const int    integ_rule   = 1;    // 1 = GSL_INTEG_GAUSS15, the 15 point Gauss-Kronrod rule

      gsl_integration_workspace * w = gsl_integration_workspace_alloc( integ_size );

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
      gsl_integration_qag( &F, lower_limit, upper_limit, abs_err_lim, rel_err_lim, integ_size, integ_rule, w, &enclosed_mass, &abs_error );

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
   if      ( r >= RArray_R[RLastIdx] )   return RArray_Phi[RLastIdx]*RArray_R[RLastIdx]/r;
   else if ( r <= RArray_R[0] )             return RArray_Phi[0];
   else                                     return Mis_InterpolateFromTable( RNBin, RArray_R, RArray_Phi, r );

} // FUNCTION : getGraviPotential


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
   const double Psi = -getGraviPotential(r);

   double *CumulativeProbability = new double [ENBin];

   CumulativeProbability[0] = EArray_DFunc[0]*sqrt( Psi-EArray_E[0] )*EArray_dE;
   for (int k=1; k<ENBin; k++)
   {
      const double Probability = ( EArray_E[k] > Psi ) ? 0 : EArray_DFunc[k]*sqrt( Psi-EArray_E[k] )*EArray_dE; //TODO: factor 2 in sqrt
      CumulativeProbability[k] = CumulativeProbability[k-1] + Probability;
   }

   const double TotalProbability          = CumulativeProbability[ELastIdx];
   const double RandomSampleProbability   = TotalProbability*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );
   //const double RandomSampleE = Mis_InterpolateFromTable( ENBin, CumulativeProbability, EArray_E, RandomSampleProbability );

   //----------------------------
   double SumProbability =  0;
   double Fraction_dEng  =  0;
   int RandomSampleIndex = -1;

   for (int k=0; k<ENBin; k++)
   {

      if ( SumProbability > RandomSampleProbability )
      {
         RandomSampleIndex = k-1;
         Fraction_dEng = (SumProbability-RandomSampleProbability)/( EArray_DFunc[RandomSampleIndex]*sqrt( Psi-EArray_E[RandomSampleIndex] )*EArray_dE ); // TODO: (1 - ...)
         break;
      }

      SumProbability += EArray_DFunc[k]*sqrt( Psi-EArray_E[k] )*EArray_dE;
   }

   if ( RandomSampleIndex < 0 )   RandomSampleIndex = ELastIdx;

   const double RandomSampleE = EArray_E[RandomSampleIndex] + EArray_dE*Fraction_dEng;
   //----------------------------

   const double RandomSampleKineticEnergy = 2*(Psi-RandomSampleE);

   const double RandomSampleVelocity      = ( RandomSampleKineticEnergy < 0.0 ) ? 0 : sqrt( RandomSampleKineticEnergy );

   delete [] CumulativeProbability;

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
      const double Psi             = Psi_Min + i*dPsi;
      const int    index_Psi       = Mis_BinarySearch_Real( RArray_Phi, 0, RLastIdx, -(Psi+0.5*dPsi) ) + 1;
      const double dsqrt_EminusPsi = ( i == N_points-1 ) ? ( - sqrt( Psi_Max-Psi ) ) : ( sqrt( Psi_Max-(Psi+dPsi) ) - sqrt( Psi_Max-Psi ) ) ;

      integral += -2*RArray_dRho_dPsi[index_Psi]*dsqrt_EminusPsi;
   }

   return integral;

} // FUNCTION : getIntegratedDistributionFunction


//-------------------------------------------------------------------------------------------------------
// Function    :  constructRadialArray
// Description :  Calculate the probability density function of particles' velocities
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructRadialArray()
{
   // Radius
   RArray_dR = Cloud_MaxR/(RNBin-1);
   for (int b=0; b<RNBin; b++)   RArray_R[b] = RArray_dR*b;

   // Density
   for (int b=1; b<RNBin; b++)   RArray_Rho[b] = getDensity( RArray_R[b] );
   RArray_Rho[0] = RArray_Rho[1];   // when r=0

   // EnclosedMass
   for (int b=1; b<RNBin; b++)   RArray_M_Enc[b] = getEnclosedMass( RArray_R[b] );
   RArray_M_Enc[0] = 0;   // when r=0

   // DensitySlope
   RArray_dRho_dR[0] =               (RArray_Rho[1] - RArray_Rho[0])/RArray_dR;

   if ( Cloud_Model == CLOUD_MODEL_TABLE )
      RArray_dRho_dR[1] =            (RArray_Rho[2] - RArray_Rho[1])/RArray_dR;
   else
      RArray_dRho_dR[1] =             Slope_LinearRegression( RArray_R, RArray_Rho,                0, 3 );

   for (int b=2; b<RNBin-2; b++)
      RArray_dRho_dR[b] =             Slope_LinearRegression( RArray_R, RArray_Rho,              b-1, 3 );

   RArray_dRho_dR[RNBin-2] = Slope_LinearRegression( RArray_R, RArray_Rho, RNBin-2, 2 );

   RArray_dRho_dR[RLastIdx]      = RArray_dRho_dR[RNBin-2];

   // GraviField
   RArray_G[0] = 0;
   for (int b=1; b<RNBin; b++)   RArray_G[b] = -NEWTON_G*RArray_M_Enc[b]/SQR( RArray_R[b] );

   // GraviPotential
   RArray_Phi[RLastIdx] = -NEWTON_G*RArray_M_Enc[RLastIdx]/RArray_R[RLastIdx];
   for (int b=RNBin-2; b>1; b--)   RArray_Phi[b] = RArray_Phi[b+1] + RArray_G[b]*RArray_dR;

   if ( Cloud_Model == CLOUD_MODEL_TABLE )   RArray_Phi[1] = RArray_Phi[2];
   else                                      RArray_Phi[1] = RArray_Phi[2] + RArray_G[1]*RArray_dR;

   RArray_Phi[0] = RArray_Phi[1];


   if ( AddExtPot_Table )
      for (int b=0; b<RNBin; b++)   RArray_Phi[b] += ExtendedInterpolatedTable( RArray_R[b], InputTable_ExtPot_nbin, InputTable_ExtPot_radius, InputTable_ExtPot_potential );
   else if ( AddExtPot_Analytical )
      for (int b=0; b<RNBin; b++)   RArray_Phi[b] += AnalyticalExternalPotential( RArray_R[b] );

   // dRho_dPsi
   for (int b=0; b<RNBin; b++)   RArray_dRho_dPsi[b] = RArray_dRho_dR[b]/RArray_G[b];

} // FUNCTION : constructRadialArray


//-------------------------------------------------------------------------------------------------------
// Function    :  constructEnergyArray
// Description :  Calculate the probability density function of particles' velocities
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructEnergyArray()
{
   EArray_MinE = -RArray_Phi[RLastIdx]; // TODO: why not zero
   if ( Cloud_Model == CLOUD_MODEL_TABLE )   EArray_MaxE = -RArray_Phi[2];
   else                                      EArray_MaxE = -RArray_Phi[1];
   EArray_dE   = (EArray_MaxE-EArray_MinE)/ENBin;

   // Set the binding energy
   EArray_E[0] = EArray_MinE;
   for (int k=1; k<ENBin; k++)   EArray_E[k] = EArray_E[k-1] + EArray_dE;

   // Set the integrated distribution function
   for (int k=0; k<ENBin; k++)   EArray_IntDFunc[k] = getIntegratedDistributionFunction( EArray_MinE, EArray_E[k], 1000 ); // TODO: why not from zero

   // Set the distribution function
   for (int k =0; k<ENBin; k++)
   {
      if      ( k <= 1         )   EArray_DFunc[k] = Slope_LinearRegression( EArray_E, EArray_IntDFunc,       0, 5 );
      else if ( k >= ENBin-2   )   EArray_DFunc[k] = Slope_LinearRegression( EArray_E, EArray_IntDFunc, ENBin-5, 5 );
      else                         EArray_DFunc[k] = Slope_LinearRegression( EArray_E, EArray_IntDFunc,     k-2, 5 );

      // check negative distribution function
      if ( EArray_DFunc[k] < 0 )   EArray_DFunc[k] = 0;
   }

   // Smooth the distribution function
   SmoothArray( EArray_DFunc, 0, ENBin );

} // FUNCTION : constructEnergyArray


//-------------------------------------------------------------------------------------------------------
// Function    :  RandomVector_GivenLength
// Description :  Compute a random 3D vector with a given length
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  Length       : Input vector length
//                RandomVector : Array to store the random 3D vector
//
// Return      :  RandomVector
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::RandomVector_GivenLength( const double Length, double RandomVector[3] )
{
   do
   {
      for (int d=0; d<3; d++)   RandomVector[d] = Random_Num_Gen->GetValue( 0, -1.0, +1.0 );
   }
   while ( SQR(RandomVector[0]) + SQR(RandomVector[1]) + SQR(RandomVector[2]) > 1.0 );

   const double Normalization = Length / sqrt( SQR(RandomVector[0]) + SQR(RandomVector[1]) + SQR(RandomVector[2]) );

   for (int d=0; d<3; d++)   RandomVector[d] *= Normalization;

} // FUNCTION : RandomVector_GivenLength



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
