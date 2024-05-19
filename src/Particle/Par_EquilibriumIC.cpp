#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

#include "Par_EquilibriumIC.h"



static void   SmoothArray                 ( double* array_x, const int index_start, const int index_end );
static double ArrayCovariance             ( const double* array_x, const double* array_y, const int index_start, const int n_elements );
static double Slope_LinearRegression      ( const double* array_x, const double* array_y, const int index_start, const int n_elements );
static double ExtendedInterpolatedTable   ( const double x, const int N, const double Table_x[], const double Table_y[] );
static double LinearDensityShellMass      ( const double r0, const double r1, const double rho0, const double rho1 );
static double AnalyticalExternalPotential ( const double r );

static double AnalyticalDensProf_Plummer  ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Plummer  ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_NFW      ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_NFW      ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Burkert  ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Burkert  ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Jaffe    ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Jaffe    ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Hernquist( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Hernquist( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Einasto  ( const double r, const double R0, const double Rho0, const double Einasto_Power_Factor );

static double MassIntegrand_Einasto       ( const double r, void* parameters );
static double MassIntegrand_Table         ( const double r, void* parameters );

// Parameters for the intergration of mass profile
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



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_EquilibriumIC
// Description :  Constructor of the class Par_Equilibrium
//
// Note        :  1. Cloud_Type is determined during the construction
//
// Parameter   :  Type : Type of this particle cloud
//
// Return      :
//-------------------------------------------------------------------------------------------------------
Par_EquilibriumIC::Par_EquilibriumIC( const char* Cloud_Type )
{

   if      ( strcmp( Cloud_Type, "Table"     ) == 0 )   Cloud_Model = CLOUD_MODEL_TABLE;
   else if ( strcmp( Cloud_Type, "Plummer"   ) == 0 )   Cloud_Model = CLOUD_MODEL_PLUMMER;
   else if ( strcmp( Cloud_Type, "NFW"       ) == 0 )   Cloud_Model = CLOUD_MODEL_NFW;
   else if ( strcmp( Cloud_Type, "Burkert"   ) == 0 )   Cloud_Model = CLOUD_MODEL_BURKERT;
   else if ( strcmp( Cloud_Type, "Jaffe"     ) == 0 )   Cloud_Model = CLOUD_MODEL_JAFFE;
   else if ( strcmp( Cloud_Type, "Hernquist" ) == 0 )   Cloud_Model = CLOUD_MODEL_HERNQUIST;
   else if ( strcmp( Cloud_Type, "Einasto"   ) == 0 )   Cloud_Model = CLOUD_MODEL_EINASTO;
   else
      Aux_Error( ERROR_INFO, "Unsupported Cloud_Type \"%s\" for Par_EquilibriumIC !!\n", Cloud_Type );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  ~Par_EquilibriumIC
// Description :  Destructor of the class Par_Equilibrium
//
// Note        :  1. Memory is free here
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
Par_EquilibriumIC::~Par_EquilibriumIC()
{
   if ( Cloud_Model == CLOUD_MODEL_TABLE )
   {
      delete [] InputTable_DensProf_radius;
      delete [] InputTable_DensProf_density;
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
   delete [] RArray_Phi;
   delete [] RArray_dRho_dPsi;

   delete [] EArray_E;
   delete [] EArray_DFunc;
   delete [] EArray_IntDFunc;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  setCenterAndBulkVel
// Description :  Set the cloud center and bulk velocity from input parameters outside
//
// Note        :  1.
//
// Parameter   :  Center_X  : x coordinate of the center
//                Center_Y  : y coordinate of the center
//                Center_Z  : z coordinate of the center
//                BulkVel_X : x component of the bulk velocity
//                BulkVel_Y : y component of the bulk velocity
//                BulkVel_Z : z component of the bulk velocity
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
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



//-------------------------------------------------------------------------------------------------------
// Function    :  setModelParameters
// Description :  Set the parameters for the analytical models from input parameters outside
//
// Note        :  1. The scale density and scale radius are general but may have different names in
//                   different models. Check the definition in AnalyticalDensProf_* for details
//
// Parameter   :  Rho0 : scale density in the density profile
//                R0   : scale radius in the density profile
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setModelParameters( const double Rho0, const double R0 )
{
   Cloud_Rho0 = Rho0;
   Cloud_R0   = R0;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  setEinastoPowerFactor
// Description :  Set the power factor in the Einasto model from input parameters outside
//
// Note        :  1. See AnalyticalDensProf_Einasto for details
//
// Parameter   :  EinastoPowerFactor  : the power factor in the Einasto density profile
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setEinastoPowerFactor( const double EinastoPowerFactor )
{
   Cloud_Einasto_Power_Factor = EinastoPowerFactor;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  setDensProfTableFilename
// Description :  Set the filename for the density profile table from input parameters outside
//
// Note        :  1. The file has two columns, the first is the radius and the second is the density
//
// Parameter   :  DensProfTableFilename : filename for the density profile table
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setDensProfTableFilename( const char* DensProfTableFilename )
{
   strcpy( DensProf_Table_Name, DensProfTableFilename );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  setParticleParameters
// Description :  Set the parameters related to the construction of the particle cloud from input parameters outside
//
// Note        :  1.
//
// Parameter   :  ParNum      : Number of particles of the particle cloud
//                MaxR        : Maximum radius for the particles in this cloud
//                Radial_NBin : Number of bins in the radial direction for the profiles of density, enclosed mass, potential, etc
//                Energy_NBin : Number of bins in the energy space for the distribution function
//                RSeed       : Random seed for setting the particle position and velocity
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setParticleParameters( const long ParNum, const double MaxR, const int Radial_NBin, const int Energy_NBin, const int RSeed )
{
   Cloud_Par_Num  = ParNum;
   Cloud_MaxR     = MaxR;
   RNBin          = Radial_NBin;
   ENBin          = Energy_NBin;
   Cloud_RSeed    = RSeed;

   if ( RNBin < 2 )   Aux_Error( ERROR_INFO, "RNBin = %d is less than 2 !!\n", RNBin );

   RLastIdx = RNBin-1;
   ELastIdx = ENBin-1;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  setExternalPotential
// Description :  Set the parametes related to the external potential from input parameters outside
//
// Note        :  1. It supports adding either anlytical external potential or external potential from table
//                2. The analytical external potential can be set at AnalyticalExternalPotential
//                3. The external potential table has two columns, the first is the radius and the second is the potential
//                4. AddExtPot_Analytical and AddExtPot_Table in Par_EquilibriumIC cannot both be turned on
//
// Parameter   :  AddingExternalPotential_Analytical  : Whether adding the analytical external potential
//                AddingExternalPotential_Table       : Whether adding the external potential from table
//                ExtPotTableFilename                 : Filename for the external potential table
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setExternalPotential( const int AddingExternalPotential_Analytical, const int AddingExternalPotential_Table, const char* ExtPotTableFilename )
{
   AddExtPot_Analytical = AddingExternalPotential_Analytical;
   AddExtPot_Table      = AddingExternalPotential_Table;

   if ( AddExtPot_Analytical  &&  AddExtPot_Table )
      Aux_Error( ERROR_INFO, "AddExtPot_Analytical and AddExtPot_Table in Par_EquilibriumIC cannot both be turned on !!\n" );

   if ( AddExtPot_Table )   strcpy( ExtPot_Table_Name, ExtPotTableFilename );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  getTotCloudMass
// Description :  Get the total enclosed mass with the radius = MaxR for this cloud
//
// Note        :  1. The enclosed mass is computed from the sampled bins of the density profile
//
// Parameter   :  None
//
// Return      :  TotCloudMass
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getTotCloudMass()
{
   return TotCloudMass;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  getParticleMass
// Description :  Get the mass of each particle in this cloud
//
// Note        :  1. ParticleMass = TotCloudMass/ParNum
//
// Parameter   :  None
//
// Return      :  ParticleMass
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getParticleMass()
{
   return ParticleMass;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  getTotCloudMassError
// Description :  Get the relative error for the total enclosed mass
//
// Note        :  1. The enclosed mass is interpolated from RArray_M_Enc
//                2. The error is calculated by comparing the enclosed mass
//                   to the analytical models or the input table interpolation.
//
// Parameter   :  None
//
// Return      :  Total Cloud Mass Error
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getTotCloudMassError()
{
   const double TotCloudMass_Analytical = getAnalEnclosedMass( Cloud_MaxR );

   return ( TotCloudMass - TotCloudMass_Analytical )/TotCloudMass_Analytical;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  loadInputDensProfTable
// Description :  Load the density profile from the input table
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  InputTable_DensProf_radius, InputTable_DensProf_density
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::loadInputDensProfTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ...\n", DensProf_Table_Name );

   const int Col_R[1] = {0};   // target column: radius
   const int Col_D[1] = {1};   // target column: density

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

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ... done\n", DensProf_Table_Name );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  loadInputExtPotTable
// Description :  Load the external potential profile from the input table
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  InputTable_ExtPot_radius, InputTable_ExtPot_potential
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::loadInputExtPotTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading ExtPot Profile Table: %s ...\n", ExtPot_Table_Name );

   const int Col_R[1] = {0};   // target column: radius
   const int Col_P[1] = {1};   // target column: potential

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
// Description :  Initialization after reading input parameter and before constructing the cloud
//
// Note        :  1. Set the random number generator
//                2. Initialize the radial arrays of physical quantities, including radius, mass, density, gravitational potential
//                3. Initialize the distribution function in the energy space
//
// Parameter   :  None
//
// Return      :  None
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

// Allocate memory
   RArray_R         = new double [RNBin];
   RArray_Rho       = new double [RNBin];
   RArray_M_Enc     = new double [RNBin];
   RArray_Phi       = new double [RNBin];
   RArray_dRho_dPsi = new double [RNBin];

   constructRadialArray();

   EArray_DFunc     = new double [ENBin];
   EArray_IntDFunc  = new double [ENBin];
   EArray_E         = new double [ENBin];

   constructEnergyArray();

//----------------------------------------------------------------------------------------------
   printf( "Radius: [" ); for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_R[b]         ); printf( "]\n");
   printf( "Dens:   [" ); for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_Rho[b]       ); printf( "]\n");
   printf( "Mass:   [" ); for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_M_Enc[b]     ); printf( "]\n");
   printf( "Pote:   [" ); for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_Phi[b]       ); printf( "]\n");
   printf( "dDdx:   [" ); for (int b=0; b<RNBin; b++)  printf(" %21.14e,", RArray_dRho_dPsi[b] ); printf( "]\n");
   printf( "BEng:   [" ); for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_E[b]         ); printf( "]\n");
   printf( "InDF:   [" ); for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_IntDFunc[b]  ); printf( "]\n");
   printf( "DisF:   [" ); for (int b=0; b<ENBin; b++)  printf(" %21.14e,", EArray_DFunc[b]     ); printf( "]\n");
//----------------------------------------------------------------------------------------------

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Initializing Par_EquilibriumIC ... done\n" );

} // FUNCTION : initialize



//-------------------------------------------------------------------------------------------------------
// Function    :  constructRadialArray
// Description :  Construct the arrays of the various radial functions
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructRadialArray()
{
// Interval of radial bins
   RArray_dR = Cloud_MaxR/(RNBin-1);

// Array of Radius, R
   for (int b=0; b<RNBin; b++)         RArray_R[b]         = RArray_dR*b;

// Array of Density, Rho(R)
   for (int b=1; b<RNBin; b++)         RArray_Rho[b]       = getDensity( RArray_R[b] );
   RArray_Rho[0]                                           = 2*RArray_Rho[1]-RArray_Rho[2];   // where r=0

// Array of Enclosed Mass, M_Enc(R) = \int_{0}^{R} Rho(r) 4\pi R^2 dR
   RArray_M_Enc[0]                                         = 0;   // where r=0
   for (int b=1; b<RNBin; b++)         RArray_M_Enc[b]     = RArray_M_Enc[b-1] + LinearDensityShellMass( RArray_R[b-1], RArray_R[b], RArray_Rho[b-1], RArray_Rho[b] );

// Array of Gravitational Potential, Phi(R) = Phi(\infty) - \int_{\infty}^{R} g dR, where g = -GM/R^2
   RArray_Phi[RLastIdx]                                    = -NEWTON_G*RArray_M_Enc[RLastIdx]/RArray_R[RLastIdx];
   for (int b=RLastIdx-1; b>0; b--)    RArray_Phi[b]       = RArray_Phi[b+1] + -NEWTON_G*0.5*(RArray_M_Enc[b]/SQR( RArray_R[b] ) + RArray_M_Enc[b+1]/SQR( RArray_R[b+1] ))*RArray_dR;
   RArray_Phi[0]                                           = RArray_Phi[1]   + -NEWTON_G*0.5*(RArray_M_Enc[1]/SQR( RArray_R[1] ))*RArray_dR;

// Adding External Potentil
   if ( AddExtPot_Table )
      for (int b=0; b<RNBin; b++)      RArray_Phi[b]      += ExtendedInterpolatedTable( RArray_R[b], InputTable_ExtPot_nbin, InputTable_ExtPot_radius, InputTable_ExtPot_potential );

   if ( AddExtPot_Analytical )
      for (int b=0; b<RNBin; b++)      RArray_Phi[b]      += AnalyticalExternalPotential( RArray_R[b] );

// Array of dRho/dPsi, dRho/dPsi = -(dRho/dPhi), where Psi = -Phi
   RArray_dRho_dPsi[0]                                     = -(RArray_Rho[1] - RArray_Rho[0])/(RArray_Phi[1] - RArray_Phi[0]);
   for (int b=1; b<RNBin-1; b++)       RArray_dRho_dPsi[b] = -(RArray_Rho[b+1] - RArray_Rho[b-1])/(RArray_Phi[b+1] - RArray_Phi[b-1]); //Slope_LinearRegression( RArray_Phi, RArray_Rho, b-1, 3 );
   RArray_dRho_dPsi[RLastIdx]                              = -(RArray_Rho[RLastIdx] - RArray_Rho[RLastIdx-1])/(RArray_Phi[RLastIdx] - RArray_Phi[RLastIdx-1]);

} // FUNCTION : constructRadialArray



//-------------------------------------------------------------------------------------------------------
// Function    :  constructEnergyArray
// Description :  Construct the energy-space arrays for the distribution function
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructEnergyArray()
{
// Ralative Energy (Binding Energy) ranges from Psi_Min to Psi_Max, where Psi = -Phi is the Relative Potential
   EArray_MinE = -RArray_Phi[RLastIdx];
   EArray_MaxE = -RArray_Phi[0];
   EArray_dE   = (EArray_MaxE-EArray_MinE)/ENBin;

// Array of Relative Energy (Binding Energy), E = -(Phi + 1/2 v^2) = Psi - 1/2 v^2 = 1/2 ( v_{esc}^2 - v^2 )
   for (int b=0; b<ENBin; b++)   EArray_E[b] = EArray_MinE + b*EArray_dE;

// Array of Integrated Distribution Function, IntDFunc(E) = \int_{E_{min}}^{E} (\frac{ 1 }{ \sqrt{E-Psi} }) (\frac{ dRho }{ dPsi }) dPsi
   for (int b=0; b<ENBin; b++)   EArray_IntDFunc[b] = getIntegratedDistributionFunction( EArray_E[b] );

// Array of Distribution Function, DFunc(E) = f(E) = d/dE IntDFunc(E)
   //for (int b=0;       b<2;       b++)   EArray_DFunc[b] = Slope_LinearRegression( EArray_E, EArray_IntDFunc,       0, 5 );
   //for (int b=2;       b<ENBin-2; b++)   EArray_DFunc[b] = Slope_LinearRegression( EArray_E, EArray_IntDFunc,     b-2, 5 );
   //for (int b=ENBin-2; b<ENBin;   b++)   EArray_DFunc[b] = Slope_LinearRegression( EArray_E, EArray_IntDFunc, ENBin-5, 5 );
   EArray_DFunc[0]                                       = (EArray_IntDFunc[1] - EArray_IntDFunc[0])/(EArray_E[1] - EArray_E[0]);
   for (int b=1;       b<ENBin-1; b++)   EArray_DFunc[b] = (EArray_IntDFunc[b+1] - EArray_IntDFunc[b-1])/(EArray_E[b+1] - EArray_E[b-1]);
   EArray_DFunc[ELastIdx]                                = (EArray_IntDFunc[ELastIdx] - EArray_IntDFunc[ELastIdx-1])/(EArray_E[ELastIdx] - EArray_E[ELastIdx-1]);

// check negative distribution function
   for (int b=0; b<ENBin; b++)
      if ( EArray_DFunc[b] < 0 )   EArray_DFunc[b] = 0;

// Smooth the distribution function
   //SmoothArray( EArray_DFunc, 0, ENBin );

} // FUNCTION : constructEnergyArray



//-------------------------------------------------------------------------------------------------------
// Function    :  constructParticles
// Description :  Set the particle's initial conditions (IC) for a cloud that is in equilibrium state
//
// Note        :  1.
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

// Determine the total enclosed mass within the maximum radius
   TotCloudMass = ExtendedInterpolatedTable( Cloud_MaxR, RNBin, RArray_R, RArray_M_Enc );
   ParticleMass = TotCloudMass/Cloud_Par_Num;

   double RandomVectorR[3];
   double RandomVectorV[3];

// Set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)
   {
//    Randomly sample the enclosed mass
      const double RandomSampleM = TotCloudMass*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

//    Ramdomly sample the radius from the enclosed mass profile with linear interpolation
      const double RandomSampleR = Mis_InterpolateFromTable( RNBin, RArray_M_Enc, RArray_R, RandomSampleM );

//    Randomly set the position vector with a given radius
      getRandomVector_GivenLength( RandomSampleR, RandomVectorR );

//    Randomly sample the velocity magnitude from the distribution function
      const double RandomSampleV = getRandomSampleVelocity( RandomSampleR );

//    Randomly set the velocity vector with the given magnitude
      getRandomVector_GivenLength( RandomSampleV, RandomVectorV );

//    Set particle attributes
      Mass_AllRank[p] = ParticleMass;
      for (int d=0; d<3; d++)   Pos_AllRank[d][p] = Cloud_Center[d]  + RandomVectorR[d];
      for (int d=0; d<3; d++)   Vel_AllRank[d][p] = Cloud_BulkVel[d] + RandomVectorV[d];

//    Check periodicity
      for (int d=0; d<3; d++)
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )   Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );

   } // for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing Par_EquilibriumIC ... done\n" );

} // FUNCTION : constructParticle



//-------------------------------------------------------------------------------------------------------
// Function    :  getDensity
// Description :  Get the density of this cloud at radius r
//
// Note        :  1. For CLOUD_MODEL_TABLE, the density is interpolated from the input density profile table
//                2. For other models, the density is from the analytical density profile
//
// Parameter   :  r : radius
//
// Return      :  Density of this cloud at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getDensity( const double r )
{
   double dens = 0.0;

   if      ( Cloud_Model == CLOUD_MODEL_TABLE     )   dens = ExtendedInterpolatedTable   ( r, InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_density );
   else if ( Cloud_Model == CLOUD_MODEL_PLUMMER   )   dens = AnalyticalDensProf_Plummer  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_NFW       )   dens = AnalyticalDensProf_NFW      ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_BURKERT   )   dens = AnalyticalDensProf_Burkert  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_JAFFE     )   dens = AnalyticalDensProf_Jaffe    ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST )   dens = AnalyticalDensProf_Hernquist( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_EINASTO   )   dens = AnalyticalDensProf_Einasto  ( r, Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor );
   else
      Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );

   return dens;

} // FUNCTION : getDensity



//-------------------------------------------------------------------------------------------------------
// Function    :  getAnalEnclosedMass
// Description :  Calculate the enclosed mass of this cloud within radius r
//
// Note        :  1. Get the enclosed mass from the analytical enclosed mass profile
//
// Parameter   :  r : radius
//
// Return      :  Enclosed mass of this cloud within radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getAnalEnclosedMass( const double r )
{
   double enclosed_mass = NULL_REAL;

   if      ( Cloud_Model == CLOUD_MODEL_PLUMMER   )   enclosed_mass = AnalyticalMassProf_Plummer  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_NFW       )   enclosed_mass = AnalyticalMassProf_NFW      ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_BURKERT   )   enclosed_mass = AnalyticalMassProf_Burkert  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_JAFFE     )   enclosed_mass = AnalyticalMassProf_Jaffe    ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST )   enclosed_mass = AnalyticalMassProf_Hernquist( r, Cloud_R0, Cloud_Rho0 );
   else
   {
#     ifdef SUPPORT_GSL
      double abs_error;

//    Arguments for the gsl integration
      const double lower_limit = 0.0;
      const double upper_limit = r;
      const double abs_err_lim = 1e-8;
      const double rel_err_lim = 1e-8;
      const int    integ_size  = 1000;
      const int    integ_rule  = 1;    // 1 = GSL_INTEG_GAUSS15, the 15 point Gauss-Kronrod rule

      gsl_integration_workspace * w = gsl_integration_workspace_alloc( integ_size );

      gsl_function F;

//    Parameters for the integrand
      struct mass_integrand_params_Einasto integrand_params_Einasto = { Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor };
      struct mass_integrand_params_Table   integrand_params_Table   = { InputTable_DensProf_nbin, InputTable_DensProf_radius, InputTable_DensProf_density };

//    Integrand for the integration
      if      ( Cloud_Model == CLOUD_MODEL_EINASTO )
      {
         F.function = &MassIntegrand_Einasto;
         F.params   = &integrand_params_Einasto;
      }
      else if ( Cloud_Model == CLOUD_MODEL_TABLE   )
      {
         F.function = &MassIntegrand_Table;
         F.params   = &integrand_params_Table;
      }
      else
         Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );

//    Integration
      gsl_integration_qag( &F, lower_limit, upper_limit, abs_err_lim, rel_err_lim, integ_size, integ_rule, w, &enclosed_mass, &abs_error );

      gsl_integration_workspace_free( w );
#     endif // #ifdef SUPPORT_GSL
   }

   return enclosed_mass;

} // FUNCTION : getAnalEnclosedMass



//-------------------------------------------------------------------------------------------------------
// Function    :  getRandomSampleVelocity
// Description :  Get the ramdomly sampled magnitude of the velocity of a particle at radius r from the distribution function
//
// Note        :  1. The probability of the magnitude of the velocity with a given radius, p(v|r) \propto v^2*f(E),
//                   where f(E) = DFunc(E) = DFunc(Psi(r)-1/2 v^2) is the distribution function
//                2. Psi = E + 1/2 v^2 -> v = \sqrt{ 2*( Psi - E ) }
//                3. dv = \frac{ 1 }{ \sqrt{ 2*( Psi - E )} } dE
//                4. The state E<0 is unbound, v^2 > v_{esc}^2
//                5. E is in the range  0 < E < Psi(r)
//
// Parameter   :  r : radius
//
// Return      :  Magnitude of velocity for a particle at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getRandomSampleVelocity( const double r )
{
// The relative potential at this radius
   const double Psi = -ExtendedInterpolatedTable( r, RNBin, RArray_R, RArray_Phi );

// CumulativeProbability = \int_{v}^{v_min} v^2 f(E) dv
//                       = \int_{E_min}^{E} (2*(Psi-E)) f(E) \frac{ 1 }{ \sqrt{ 2*( Psi - E )} } dE
//                       = \int_{E_min}^{E} \sqrt{ (2*(Psi-E)) } f(E) dE
   double *CumulativeProbability = new double [ENBin];
   double Probability;

   CumulativeProbability[0] = 0;
   for (int b=1; b<ENBin; b++)
   {
      if ( EArray_E[b] > Psi )   Probability = ( EArray_E[b-1] > Psi ) ? 0 : 0.5*EArray_DFunc[b-1]*sqrt(2*(Psi-EArray_E[b-1]))*(Psi-EArray_E[b-1]);
      else                       Probability = 0.5*( EArray_DFunc[b-1]*sqrt(2*(Psi-EArray_E[b-1])) +
                                                     EArray_DFunc[b  ]*sqrt(2*(Psi-EArray_E[b  ])) )*EArray_dE;

      CumulativeProbability[b] = CumulativeProbability[b-1] + Probability;
   }

   const double TotalProbability        = CumulativeProbability[ELastIdx];
   const double RandomSampleProbability = TotalProbability*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );
   const double RandomSampleE           = Mis_InterpolateFromTable( ENBin, CumulativeProbability, EArray_E, RandomSampleProbability );

// v^2 = 2*(Psi-E)
   const double RandomSampleVelocity = ( RandomSampleE > Psi ) ? 0 : sqrt( 2*(Psi-RandomSampleE) );

   delete [] CumulativeProbability;

   return RandomSampleVelocity;

} // FUNCTION : getRandomSampleVelocity



//-------------------------------------------------------------------------------------------------------
// Function    :  getIntegratedDistributionFunction
// Description :  Get the integrated distribution function
//
// Note        :  1. The integrated deistribution function
//                   = \int_{0}^{\mathcal{E}} \frac{ 1 }{ \sqrt{ \mathcal{E} -\Psi } } \frac{ d\rho }{ d\Psi } d\Psi
//                   = \int_{0}^{\mathcal{E}} -2 \frac{ d\rho }{ d\Psi } d\sqrt{ \mathcal{E} -\Psi } }
//
// Parameter   :  E  : \mathcal{E} in the above equation
//
// Return      :  Integrated distribution function
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getIntegratedDistributionFunction( const double E )
{
   if ( E <= EArray_MinE )   return 0.0;

   const int    N_points = 10000;
   const double dPsi     = ( E - EArray_MinE )/N_points;

   double integral = 0;
   for (int i=0; i<N_points; i++)
   {
      const double Psi             = EArray_MinE + i*dPsi;
      const double dRho_dPsi       = Mis_InterpolateFromTable( RNBin, RArray_Phi, RArray_dRho_dPsi, -(Psi+0.5*dPsi) );
      const double dsqrt_EminusPsi = ( ( (Psi+dPsi) >= E ) ? 0 : sqrt( E-(Psi+dPsi) ) ) - sqrt( E-Psi );

      integral += -2*dRho_dPsi*dsqrt_EminusPsi;
   }

   return integral;

} // FUNCTION : getIntegratedDistributionFunction



//-------------------------------------------------------------------------------------------------------
// Function    :  getRandomVector_GivenLength
// Description :  Get a randomly sampled 3D vector with a given length
//
// Note        :  Uniformly random sample in theta and phi does NOT give a uniformly random sample in 3D space
//                --> Uniformly random sample in a 3D sphere and then normalize all vectors to the given radius
//
// Parameter   :  Length       : Input vector length
//                RandomVector : Array to store the random 3D vector
//
// Return      :  RandomVector
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::getRandomVector_GivenLength( const double Length, double RandomVector[3] )
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
// Function    :  SmoothArray
// Description :  Smooth the input array
//
// Note        :  1.
//
// Parameter   :  array_x     : array of x data
//                index_start : the first index in the array for the smoothing
//                index_end   : the last index in the array for the smoothing
//
// Return      :  array_x
//-------------------------------------------------------------------------------------------------------
void SmoothArray( double* array_x, const int index_start, const int index_end )
{
   int    smoothing_n_elements = 10; // smoothing every "smoothing_n_elements" elements
   double smoothing_criterion  =  3; // smoothing when the ratio is larger than this criterion

// Set the elements as zero if its ratio to other elements is larger than smoothing_criterion
   for (int i=index_start; i<index_end-smoothing_n_elements+1; i++)
   for (int j=i; j<i+smoothing_n_elements; j++)
   for (int k=i; k<i+smoothing_n_elements; k++)
   {
      if ( array_x[k] != 0  &&  fabs( array_x[j]/array_x[k] ) > smoothing_criterion )   array_x[j] = 0;
   }

// Set those zero elements as the average of non-zero elements
   for (int i=index_start; i<index_end-smoothing_n_elements+1; i++)
   {
      double sum_of_nonzero = 0;
      int    num_of_nonzero = 0;
      double ave_of_nonzero = 0;

//    Sum the non-zero elements
      for (int j=i; j<i+smoothing_n_elements; j++)
      {
         if ( array_x[j] != 0 )
         {
            sum_of_nonzero += array_x[j];
            num_of_nonzero ++;
         }
      }

//    Average of non-zero elements
      if ( num_of_nonzero != 0 )   ave_of_nonzero = sum_of_nonzero/num_of_nonzero;

//    Assign the average of non-zero elements to the zero element
      for (int j=i; j<i+smoothing_n_elements; j++)
         if ( array_x[j] == 0 )   array_x[j] = ave_of_nonzero;
   }

} // FUNCTION : SmoothArray



//-------------------------------------------------------------------------------------------------------
// Function    :  ArrayCovariance
// Description :  Get the covariance between two arrays
//
// Note        :  1. if x == y, then the covariance is the variance of x
//
// Parameter   :  array_x     : array of x data
//                array_y     : array of y data
//                index_start : the first index in the array for the linear regression
//                n_elements  : number of elements for the linear regression
//
// Return      :  covariance_xy
//-------------------------------------------------------------------------------------------------------
double ArrayCovariance( const double* array_x, const double* array_y, const int index_start, const int n_elements )
{
   const double normalized_factor = 1.0/n_elements;

// Average
   double average_x = 0.0;
   double average_y = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)
   {
      average_x += array_x[i];
      average_y += array_y[i];
   }

   average_x *= normalized_factor;
   average_y *= normalized_factor;

// Covariance
   double covariance_xy = 0.0;

   for (int i=index_start; i<index_start+n_elements; i++)   covariance_xy += (array_x[i]-average_x)*(array_y[i]-average_y);

   covariance_xy *= normalized_factor;


   return covariance_xy;

} // FUNCTION : ArrayCovariance



//-------------------------------------------------------------------------------------------------------
// Function    :  Slope_LinearRegression
// Description :  Get the slope of y-x using linear regression
//
// Note        :
//
// Parameter   :  array_x     : array of x data
//                array_y     : array of y data
//                index_start : the first index in the array for the linear regression
//                n_elements  : number of elements for the linear regression
//
// Return      :  slope_y
//-------------------------------------------------------------------------------------------------------
double Slope_LinearRegression( const double* array_x, const double* array_y, const int index_start, const int n_elements )
{

   const double variance_x    = ArrayCovariance( array_x, array_x, index_start, n_elements );
   const double covariance_xy = ArrayCovariance( array_x, array_y, index_start, n_elements );

   const double slope_y       = covariance_xy/variance_x;

   return slope_y;

} // FUNCTION : Slope_LinearRegression



//-------------------------------------------------------------------------------------------------------
// Function    :  ExtendedInterpolatedTable
// Description :  Get the interpolated value even x is out of range
//
// Note        :  1. If the x is in the range of Table_x, it call Mis_InterpolateFromTable
//                2. If the x is out of the range of Table_x, it return the values at end points
//
// Parameter   :  x       : Position to get the interpolation
//                N       : Number of bins in the Table
//                Table_x : Table of values to be input
//                Table_y : Table of values to be interpolated
//
// Return      :  Interpolated value
//-------------------------------------------------------------------------------------------------------
double ExtendedInterpolatedTable( const double x, const int N, const double Table_x[], const double Table_y[] )
{
   if      ( x <= Table_x[0]   )   return Table_y[0];
   else if ( x >= Table_x[N-1] )   return Table_y[N-1];
   else                            return Mis_InterpolateFromTable( N, Table_x, Table_y, x );
}



//-------------------------------------------------------------------------------------------------------
// Function    :  LinearDensityShellMass
// Description :  Get the shell mass between to radii according to a linear density profile
//
// Note        :  1. Assume the density profile \rho(r) is linear between two end points, r0 and r1:
//                   \rho(r) = \rho_0 + (\frac{(\rho_1-\rho_0)}{r_1-r_0})(r-r_0)
//                   Then, the integrated shell mass between the two end points is
//                   M_{shell} = \int_{r_0}^{r_1} \rho(r) 4\pi r^2 dr
//                             = 4\pi r_0^2 \Delta r   [ \frac{1}{2} \rho_0 + \frac{1}{2} \rho_1 ] +
//                               4\pi r_0   \Delta r^2 [ \frac{1}{3} \rho_0 + \frac{2}{3} \rho_1 ] +
//                               4\pi       \Delta r^3 [ \frac{1}{12}\rho_0 + \frac{3}{12}\rho_1 ]
//
// Parameter   :  r0   : Radius of the start point
//                r1   : Radius of the end point
//                rho0 : Density of the start point
//                rho1 : Density of the end point
//
// Return      :  Shell mass integrated from the linear density profile
//-------------------------------------------------------------------------------------------------------
double LinearDensityShellMass( const double r0, const double r1, const double rho0, const double rho1 )
{
   const double dr  = r1 - r0;

   return M_PI*dr*( r0*r0*( 6*rho0 + 6*rho1 ) + r0*dr*( 4*rho0 + 8*rho1 ) + dr*dr*( rho0 + 3*rho1 ) )/3.0;

} // FUNCTION : LinearDensityShellMass



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalExternalPotential
// Description :  Analytical external potential
//
// Note        :
//
// Parameter   :  r : input radius
//
// Return      :  external potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalExternalPotential( const double r )
{
   return 0.0;

} // FUNCTION : AnalyticalExternalPotential



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
// Function    :  MassIntegrand_Table
// Description :  Integrand for the enclosed mass profile of the table
//
// Note        :  integrand = 4*\pi*r^2*\rho(r)
//
// Parameter   :  r          : input radius
//                parameters : parameters for the model
//
// Return      :  integrand of mass at the given radius
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



#endif // #ifdef PARTICLE
