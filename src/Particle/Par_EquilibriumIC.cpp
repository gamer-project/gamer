#include "GAMER.h"

#ifdef MASSIVE_PARTICLES

#include "Par_EquilibriumIC.h"



static double ExtendedInterpolatedTable   ( const double x, const int N, const double Table_x[], const double Table_y[] );
static double LinearDensityShellMass      ( const double r0, const double r1, const double rho0, const double rho1 );
static double UserDefAnlaytical_ExtPot    ( const double r );

static double AnalyticalDensProf_Plummer  ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Plummer  ( const double r, const double R0, const double Rho0 );
static double AnalyticalPoteProf_Plummer  ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_NFW      ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_NFW      ( const double r, const double R0, const double Rho0 );
static double AnalyticalPoteProf_NFW      ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Burkert  ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Burkert  ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Jaffe    ( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Jaffe    ( const double r, const double R0, const double Rho0 );
static double AnalyticalPoteProf_Jaffe    ( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Hernquist( const double r, const double R0, const double Rho0 );
static double AnalyticalMassProf_Hernquist( const double r, const double R0, const double Rho0 );
static double AnalyticalPoteProf_Hernquist( const double r, const double R0, const double Rho0 );
static double AnalyticalDensProf_Einasto  ( const double r, const double R0, const double Rho0, const double Einasto_Power_Factor );

static double MassIntegrand_Einasto       ( const double r, void* parameters );
static double MassIntegrand_Table         ( const double r, void* parameters );

// Parameters for the intergration of mass profile
struct mass_integrand_params_Einasto{ double Cloud_R0; double Cloud_Rho0; double Cloud_Einasto_Power_Factor; };
struct mass_integrand_params_Table  { int NBin;        double* Table_R;   double* Table_D;                   };



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_EquilibriumIC
// Description :  Constructor of the class Par_EquilibriumIC
//
// Note        :  1. Cloud_Type is determined during the construction
//
// Parameter   :  Cloud_Type : Type of this particle cloud
//
// Return      :  None
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
   else   Aux_Error( ERROR_INFO, "Unsupported Cloud_Type \"%s\" for Par_EquilibriumIC !!\n", Cloud_Type );

} // FUNCTION : Par_EquilibriumIC



//-------------------------------------------------------------------------------------------------------
// Function    :  ~Par_EquilibriumIC
// Description :  Destructor of the class Par_EquilibriumIC
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
      delete [] InputTable_DensProf_Radius;
      delete [] InputTable_DensProf_Density;
   }

   if ( AddExtPot_Table )
   {
      delete [] InputTable_ExtPot_Radius;
      delete [] InputTable_ExtPot_Potential;
   }

   delete [] RArray_R;
   delete [] RArray_Rho;
   delete [] RArray_M_Enc;
   delete [] RArray_Phi;
   delete [] RArray_dRho_dPsi;

   delete [] EArray_E;
   delete [] EArray_DFunc;
   delete [] EArray_IntDFunc;

   delete Random_Num_Gen;

} // FUNCTION : ~Par_EquilibriumIC



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
// Return      :  Cloud_Center[], Cloud_BulkVel[]
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

} // FUNCTION : setCenterAndBulkVel



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
// Return      :  Cloud_Rho0, Cloud_R0
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setModelParameters( const double Rho0, const double R0 )
{
   Cloud_Rho0 = Rho0;
   Cloud_R0   = R0;

} // FUNCTION : setModelParameters



//-------------------------------------------------------------------------------------------------------
// Function    :  setEinastoPowerFactor
// Description :  Set the power factor in the Einasto model from input parameters outside
//
// Note        :  1. See AnalyticalDensProf_Einasto for details
//
// Parameter   :  EinastoPowerFactor : the power factor in the Einasto density profile
//
// Return      :  Cloud_Einasto_Power_Factor
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setEinastoPowerFactor( const double EinastoPowerFactor )
{
   Cloud_Einasto_Power_Factor = EinastoPowerFactor;

} // FUNCTION : setEinastoPowerFactor



//-------------------------------------------------------------------------------------------------------
// Function    :  setDensProfTableFilename
// Description :  Set the filename for the density profile table from input parameters outside
//
// Note        :  1. The file has two columns: the first is the radius and the second is the density
//
// Parameter   :  DensProfTableFilename : filename for the density profile table
//
// Return      :  DensProf_Table_Name
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setDensProfTableFilename( const char* DensProfTableFilename )
{
   strcpy( DensProf_Table_Name, DensProfTableFilename );

} // FUNCTION : setDensProfTableFilename



//-------------------------------------------------------------------------------------------------------
// Function    :  setParticleParameters
// Description :  Set the parameters related to the construction of the particle cloud from input parameters outside
//
// Note        :  1.
//
// Parameter   :  ParNum : Number of particles of the particle cloud
//                MaxR   : Maximum radius for the scattered particles in this cloud
//                NBin   : Number of bins of radial profiles inside the MaxR
//                RSeed  : Random seed for setting the particle position and velocity
//
// Return      :  Cloud_Par_Num, Cloud_MaxR, Cloud_NBin, Cloud_RSeed
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setParticleParameters( const long ParNum, const double MaxR, const int NBin, const int RSeed )
{
   Cloud_Par_Num  = ParNum;
   Cloud_MaxR     = MaxR;
   Cloud_NBin     = NBin;
   Cloud_RSeed    = RSeed;

   if ( Cloud_Par_Num < 1 )   Aux_Error( ERROR_INFO, "Cloud_Par_Num = %ld is less than 1 !!\n", Cloud_Par_Num );
   if ( Cloud_MaxR <= 0.0 )   Aux_Error( ERROR_INFO, "Cloud_MaxR = % 14.7e is not positive !!\n", Cloud_MaxR );
   if ( Cloud_NBin < 1    )   Aux_Error( ERROR_INFO, "Cloud_NBin = %d is less than 1 !!\n", Cloud_NBin       );

} // FUNCTION : setParticleParameters



//-------------------------------------------------------------------------------------------------------
// Function    :  setExtPotParameters
// Description :  Set the parametes related to the external potential from input parameters outside
//
// Note        :  1. It supports adding either analytical external potential or external potential from a table
//                2. The analytical external potential can be set via getExternalPotential()
//                3. The external potential table has two columns: the first is the radius and the second is the potential
//                4. AddExtPot_Analytical and AddExtPot_Table in Par_EquilibriumIC cannot both be turned on
//
// Parameter   :  AddingExternalPotential_Analytical : Whether adding an analytical external potential
//                AddingExternalPotential_Table      : Whether adding an external potential from a table
//                ExtPotTableFilename                : Filename for the external potential table
//
// Return      :  AddExtPot_Analytical, AddExtPot_Table, ExtPot_Table_Name
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::setExtPotParameters( const int AddingExternalPotential_Analytical, const int AddingExternalPotential_Table, const char* ExtPotTableFilename )
{
   AddExtPot_Analytical = AddingExternalPotential_Analytical;
   AddExtPot_Table      = AddingExternalPotential_Table;

   if ( AddExtPot_Analytical  &&  AddExtPot_Table )
      Aux_Error( ERROR_INFO, "AddExtPot_Analytical and AddExtPot_Table in Par_EquilibriumIC cannot both be turned on !!\n" );

   if ( AddExtPot_Table )   strcpy( ExtPot_Table_Name, ExtPotTableFilename );

} // FUNCTION : setExtPotParameters



//-------------------------------------------------------------------------------------------------------
// Function    :  loadInputDensProfTable
// Description :  Load the density profile from the input table
//
// Note        :  1. Called by constructDistribution()
//
// Parameter   :  None
//
// Return      :  InputTable_DensProf_Radius, InputTable_DensProf_Density, InputTable_DensProf_NBin
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::loadInputDensProfTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ...\n", DensProf_Table_Name );

   const int Col_R[1] = {0};   // target column: radius
   const int Col_D[1] = {1};   // target column: density

   const int NRowR = Aux_LoadTable( InputTable_DensProf_Radius,  DensProf_Table_Name, 1, Col_R, true, true );
   const int NRowD = Aux_LoadTable( InputTable_DensProf_Density, DensProf_Table_Name, 1, Col_D, true, true );

// Check the number of rows is consistent
   if ( NRowR != NRowD )
      Aux_Error( ERROR_INFO, "The number of rows of density (%d) is not equal to the number of rows of radii (%d) in density table %s !!\n",
                             NRowD, NRowR, DensProf_Table_Name );
   else
      InputTable_DensProf_NBin = NRowR;

// Check maximum radius in the density table must be larger than Cloud_MaxR
   if ( InputTable_DensProf_Radius[InputTable_DensProf_NBin-1] < Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Maximum radius (%14.7e) in density table %s is smaller then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_DensProf_Radius[InputTable_DensProf_NBin-1], DensProf_Table_Name, Cloud_MaxR );

// Check minimum radius in the density table must be smaller than Cloud_MaxR
   if ( InputTable_DensProf_Radius[0] > Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Minimum radius (%14.7e) in density table %s is larger then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_DensProf_Radius[0], DensProf_Table_Name, Cloud_MaxR );

// Check the properties of values in the density table
   for (int i=1; i<InputTable_DensProf_NBin; i++)
   {
      if ( InputTable_DensProf_Radius[i] <= InputTable_DensProf_Radius[i-1] )
         Aux_Error( ERROR_INFO, "Radii in density table \"%s\" are not strictly increasing !!\n", DensProf_Table_Name );

      if ( InputTable_DensProf_Density[i] < 0.0 )
         Aux_Error( ERROR_INFO, "Densities in density table \"%s\" have negative value !!\n", DensProf_Table_Name );

      if ( InputTable_DensProf_Density[i] > InputTable_DensProf_Density[i-1] )
         Aux_Error( ERROR_INFO, "Densities in density table \"%s\" are not monotonically decreasing !!\n", DensProf_Table_Name );
   }

// Check Cloud_NBin must be larger than the number of bins inside Cloud_MaxR in the density table
   int InputTable_NBins_InsideMaxR = 0;
   for (int i=0; i<InputTable_DensProf_NBin; i++)
   {
      if ( InputTable_DensProf_Radius[i] > Cloud_MaxR )   break;

      InputTable_NBins_InsideMaxR += 1;
   }
   if ( Cloud_NBin < InputTable_NBins_InsideMaxR  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : Cloud_NBin = %d is smaller than the number of bins (= %d) inside Cloud_MaxR (= %14.7e) in input density table \"%s\" !!\n",
                           Cloud_NBin, InputTable_NBins_InsideMaxR, Cloud_MaxR, DensProf_Table_Name );

// Reset the maximum radius of RArray according to the density table
   RArray_MaxR = MIN( RArray_MaxR, InputTable_DensProf_Radius[InputTable_DensProf_NBin-1] );

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading Density Profile Table: \"%s\" ... done\n", DensProf_Table_Name );

} // FUNCTION : loadInputDensProfTable



//-------------------------------------------------------------------------------------------------------
// Function    :  loadInputExtPotTable
// Description :  Load the external potential profile from the input table
//
// Note        :  1. Called by constructDistribution()
//
// Parameter   :  None
//
// Return      :  InputTable_ExtPot_Radius, InputTable_ExtPot_Potential, InputTable_ExtPot_NBin
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::loadInputExtPotTable()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading ExtPot Profile Table: %s ...\n", ExtPot_Table_Name );

   const int Col_R[1] = {0};   // target column: radius
   const int Col_P[1] = {1};   // target column: potential

   const int NRowR = Aux_LoadTable( InputTable_ExtPot_Radius,    ExtPot_Table_Name, 1, Col_R, true, true );
   const int NRowP = Aux_LoadTable( InputTable_ExtPot_Potential, ExtPot_Table_Name, 1, Col_P, true, true );

// Check the number of rows is consistent
   if ( NRowR != NRowP )
      Aux_Error( ERROR_INFO, "The number of rows of potential (%d) is not equal to the number of rows of radii (%d) in ExtPot table %s !!\n",
                             NRowP, NRowR, ExtPot_Table_Name );
   else
      InputTable_ExtPot_NBin = NRowR;

// Check maximum radius in the external potential table must be larger than Cloud_MaxR
   if ( InputTable_ExtPot_Radius[InputTable_ExtPot_NBin-1] < Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Maximum radius (%14.7e) in external potential table %s is smaller then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_ExtPot_Radius[InputTable_ExtPot_NBin-1], ExtPot_Table_Name, Cloud_MaxR );

// Check minimum radius in the external potential table must be smaller than Cloud_MaxR
   if ( InputTable_ExtPot_Radius[0] > Cloud_MaxR )
      Aux_Error( ERROR_INFO, "Minimum radius (%14.7e) in external potential table %s is larger then Cloud_MaxR (%14.7e) !!\n",
                             InputTable_ExtPot_Radius[0], ExtPot_Table_Name, Cloud_MaxR );

// Check the properties of values in the external potential table
   for (int i=1; i<InputTable_ExtPot_NBin; i++)
   {
      if ( InputTable_ExtPot_Radius[i] <= InputTable_ExtPot_Radius[i-1] )
         Aux_Error( ERROR_INFO, "Radii in external potential table \"%s\" are not strictly increasing !!\n", ExtPot_Table_Name );

      // The potential at r=infinity is assumed to be zero
      if ( InputTable_ExtPot_Potential[i] > 0.0 )
         Aux_Error( ERROR_INFO, "Potential in external potential table \"%s\" have postive value !!\n", ExtPot_Table_Name );

      // Assume there is no outward force
      if ( InputTable_ExtPot_Potential[i] < InputTable_ExtPot_Potential[i-1] )
         Aux_Error( ERROR_INFO, "Potential in external potential table \"%s\" are not monotonically increasing !!\n", ExtPot_Table_Name );
   }

// Reset the maximum radius of RArray according to the external potential table
   RArray_MaxR = MIN( RArray_MaxR, InputTable_ExtPot_Radius[InputTable_ExtPot_NBin-1] );

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "   Loading ExtPot Profile Table: \"%s\" ... done\n", ExtPot_Table_Name );

} // FUNCTION : loadInputExtPotTable



//-------------------------------------------------------------------------------------------------------
// Function    :  constructDistribution
// Description :  Initialization after reading input parameters and before constructing the clouds
//
// Note        :  1. Set the random number generator
//                2. Construct the radial arrays of physical quantities, including radius, mass, density, and gravitational potential
//                3. Construct the distribution function in the energy space
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructDistribution()
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing the distribution in Par_EquilibriumIC ...\n" );

// Maximum radius for the profiles to calculate the distribution function, should be much larger than the Cloud_MaxR
   RArray_MaxR = 1.0e3*Cloud_MaxR;

// Load the input density table
   if ( Cloud_Model == CLOUD_MODEL_TABLE )   loadInputDensProfTable();
   if ( AddExtPot_Table )                    loadInputExtPotTable();

// Construct RArray
   RNPoints         = MIN( 10000000, MAX( 1000, (int)ceil( Cloud_NBin*(RArray_MaxR/Cloud_MaxR)+1 ) ) ); // within a reasonable range, try to keep roughly the same Cloud_NBin inside Cloud_MaxR
   RLastIdx         = RNPoints-1;
   RArray_R         = new double [RNPoints];
   RArray_Rho       = new double [RNPoints];
   RArray_M_Enc     = new double [RNPoints];
   RArray_Phi       = new double [RNPoints];
   RArray_dRho_dPsi = new double [RNPoints];

   constructRadialArray();

// Construct EArray
   ENPoints         = MIN( 10000000, MAX( 1000, (int)ceil( 2*(RArray_Phi[0]/RArray_Phi[RLastIdx])+1 ) ) ); // within a reasonable range, estimated by the potential difference at the center
   ELastIdx         = ENPoints-1;
   EArray_DFunc     = new double [ENPoints];
   EArray_IntDFunc  = new double [ENPoints];
   EArray_E         = new double [ENPoints];

   constructEnergyArray();

#  ifdef GAMER_DEBUG
// Output the arrays for debugging
   printArrays();
#  endif // #ifdef GAMER_DEBUG

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing the distribution in Par_EquilibriumIC ... done\n" );

} // FUNCTION : constructDistribution



//-------------------------------------------------------------------------------------------------------
// Function    :  constructRadialArray
// Description :  Construct the arrays of various radial functions
//
// Note        :  1. Called by constructDistribution()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructRadialArray()
{
// Interval of radial bins
   RArray_dR = RArray_MaxR/(RNPoints-1);

// Array of Radius, R
   for (int i=0; i<RNPoints; i++)     RArray_R[i]         = RArray_dR*i;

// Array of Density, Rho(R)
   for (int i=1; i<RNPoints; i++)     RArray_Rho[i]       = getDensity( RArray_R[i] );
   RArray_Rho[0]                                          = 2.0*RArray_Rho[1]-RArray_Rho[2];   // where r=0

// Array of Enclosed Mass, M_Enc(R) = \int_{0}^{R} Rho(r) 4\pi R^2 dR
   RArray_M_Enc[0]                                        = 0;   // where r=0
   for (int i=1; i<RNPoints; i++)     RArray_M_Enc[i]     = RArray_M_Enc[i-1] + LinearDensityShellMass( RArray_R[i-1], RArray_R[i], RArray_Rho[i-1], RArray_Rho[i] );

// Array of Gravitational Potential, Phi(R) = Phi(\infty) - \int_{\infty}^{R} g dR, where g = -GM/R^2
   RArray_Phi[RLastIdx]                                   = -NEWTON_G*RArray_M_Enc[RLastIdx]/RArray_R[RLastIdx];
   for (int i=RLastIdx-1; i>0; i--)   RArray_Phi[i]       = RArray_Phi[i+1] + -NEWTON_G*0.5*(RArray_M_Enc[i]/SQR( RArray_R[i] ) + RArray_M_Enc[i+1]/SQR( RArray_R[i+1] ))*RArray_dR;
   RArray_Phi[0]                                          = RArray_Phi[1]   + -NEWTON_G*0.5*(RArray_M_Enc[1]/SQR( RArray_R[1] ))*RArray_dR;

// Adding External Potentil
   for (int i=0; i<RNPoints; i++)     RArray_Phi[i]      += getExternalPotential( RArray_R[i] );

// Array of dRho/dPsi, dRho/dPsi = -(dRho/dPhi), where Psi = -Phi
   RArray_dRho_dPsi[0]                                    = -(RArray_Rho[1]        - RArray_Rho[0]         )/(RArray_Phi[1]        - RArray_Phi[0]         );
   for (int i=1; i<RNPoints-1; i++)   RArray_dRho_dPsi[i] = -(RArray_Rho[i+1]      - RArray_Rho[i-1]       )/(RArray_Phi[i+1]      - RArray_Phi[i-1]       );
   RArray_dRho_dPsi[RLastIdx]                             = -(RArray_Rho[RLastIdx] - RArray_Rho[RLastIdx-1])/(RArray_Phi[RLastIdx] - RArray_Phi[RLastIdx-1]);

} // FUNCTION : constructRadialArray



//-------------------------------------------------------------------------------------------------------
// Function    :  constructEnergyArray
// Description :  Construct the energy-space arrays for the distribution function
//
// Note        :  1. Called by constructDistribution()
//                2. Solve the ergodic distribution function from the density profile using the Eddington inversion
//                3. Reference: Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Chapter 4.3 -- Chapter 4.3.1
//                              Eddingtion A. S., 1916, MNRAS, doi:10.1093/mnras/76.7.572
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::constructEnergyArray()
{
// Ralative Energy (Binding Energy) E ranges from Psi_Min to Psi_Max, where Psi = -Phi is the Relative Potential
   EArray_MinE = 0.0;
   EArray_MaxE = -RArray_Phi[0];
   EArray_dE   = (EArray_MaxE-EArray_MinE)/(ENPoints-1);

// Array of Relative Energy (Binding Energy), E = -(Phi + 1/2 v^2) = Psi - 1/2 v^2 = 1/2 ( v_{esc}^2 - v^2 ). See Eq.(4.41) in Binney & Tremaine 2008, Galactic Dynamics
   for (int i=0; i<ENPoints; i++)     EArray_E[i]        = EArray_MinE + i*EArray_dE;

// Array of Integrated Distribution Function, IntDFunc(E) = \frac{1}{\sqrt{8}\pi^2} \int_{0}^{E} (\frac{ 1 }{ \sqrt{E-Psi} }) (\frac{ dRho }{ dPsi }) dPsi. See Eq.(4.46a) in Binney & Tremaine 2008, Galactic Dynamics
   for (int i=0; i<ENPoints; i++)     EArray_IntDFunc[i] = getIntegratedDistributionFunction( EArray_E[i] );

// Array of Distribution Function, DFunc(E) = f(E) = d/dE IntDFunc(E), which is the Eddington's formula Eq.(4.46a) in Binney & Tremaine 2008, Galactic Dynamics
   EArray_DFunc[0]                                       = (EArray_IntDFunc[1]        - EArray_IntDFunc[0]         )/EArray_dE;
   for (int i=1; i<ENPoints-1; i++)   EArray_DFunc[i]    = (EArray_IntDFunc[i+1]      - EArray_IntDFunc[i-1]       )/(2.0*EArray_dE);
   EArray_DFunc[ELastIdx]                                = (EArray_IntDFunc[ELastIdx] - EArray_IntDFunc[ELastIdx-1])/EArray_dE;

// check negative distribution function
   for (int i=0; i<ENPoints; i++)
      if ( EArray_DFunc[i] < 0 )      EArray_DFunc[i]    = 0.0;

} // FUNCTION : constructEnergyArray



//-------------------------------------------------------------------------------------------------------
// Function    :  printArrays
// Description :  Output the RArray and EArray to file for debugging purpose
//
// Note        :  1.
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_EquilibriumIC::printArrays()
{
   if ( MPI_Rank == 0 )
   {
      static int Cloud_ID = 1;

//    RArray
      char Filename_R[MAX_STRING];
      sprintf( Filename_R, "%s%d%s%d", "Record__ParEquilibriumIC_Model_", Cloud_Model, "_RArray_", Cloud_ID );
      FILE *File_R = fopen( Filename_R, "w" );

      fprintf( File_R, "#%21s %21s %21s %21s %21s\n", "R", "Rho", "M_Enc", "Phi", "dRho_dPsi" );

      for (int i=0; i<RNPoints; i++)
         fprintf( File_R, " %21.14e %21.14e %21.14e %21.14e %21.14e\n",
                          RArray_R[i], RArray_Rho[i], RArray_M_Enc[i], RArray_Phi[i], RArray_dRho_dPsi[i] );

      fclose( File_R );

//    EArray
      char Filename_E[MAX_STRING];
      sprintf( Filename_E, "%s%d%s%d", "Record__ParEquilibriumIC_Model_", Cloud_Model, "_EArray_", Cloud_ID );
      FILE *File_E = fopen( Filename_E, "w" );

      fprintf( File_E, "#%21s %21s %21s\n", "E", "IntDFunc", "DFunc" );

      for (int i=0; i<ENPoints; i++)
         fprintf( File_E, " %21.14e %21.14e %21.14e\n", EArray_E[i], EArray_IntDFunc[i], EArray_DFunc[i] );

      fclose( File_E );

      Cloud_ID += 1;
   }

} // FUNCTION : printArrays



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
void Par_EquilibriumIC::constructParticles( real_par *Mass_AllRank, real_par *Pos_AllRank[3], real_par *Vel_AllRank[3], const long Par_Idx0 )
{
   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing the particles in Par_EquilibriumIC ...\n" );

// Determine the total enclosed mass within the maximum radius
   const double TotCloudMass_Analytical = getAnalEnclosedMass( Cloud_MaxR );
   TotCloudMass                         = ExtendedInterpolatedTable( Cloud_MaxR, RNPoints, RArray_R, RArray_M_Enc );
   ParticleMass                         = TotCloudMass/Cloud_Par_Num;
   TotCloudMassError                    = ( TotCloudMass - TotCloudMass_Analytical )/TotCloudMass_Analytical;

// Set random number generator
   Random_Num_Gen = new RandomNumber_t( 1 );
   Random_Num_Gen->SetSeed( 0, Cloud_RSeed );

   double RandomVectorR[3];
   double RandomVectorV[3];

   CumulProbaDistr_GivenRadius = new double [ENPoints];

// Set particle attributes
   for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)
   {
//    Randomly sample the enclosed mass
      const double RandomSampleM = TotCloudMass*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );

//    Ramdomly sample the radius from the enclosed mass profile with linear interpolation
      const double RandomSampleR = Mis_InterpolateFromTable( RNPoints, RArray_M_Enc, RArray_R, RandomSampleM );

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
         if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )   Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real_par)amr->BoxSize[d], (real_par)amr->BoxSize[d] );

   } // for (long p=Par_Idx0; p<Par_Idx0+Cloud_Par_Num; p++)

   delete [] CumulProbaDistr_GivenRadius;

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Constructing the particles in Par_EquilibriumIC ... done\n" );

} // FUNCTION : constructParticles



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

   if      ( Cloud_Model == CLOUD_MODEL_TABLE     )   dens = ExtendedInterpolatedTable   ( r, InputTable_DensProf_NBin, InputTable_DensProf_Radius, InputTable_DensProf_Density );
   else if ( Cloud_Model == CLOUD_MODEL_PLUMMER   )   dens = AnalyticalDensProf_Plummer  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_NFW       )   dens = AnalyticalDensProf_NFW      ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_BURKERT   )   dens = AnalyticalDensProf_Burkert  ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_JAFFE     )   dens = AnalyticalDensProf_Jaffe    ( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_HERNQUIST )   dens = AnalyticalDensProf_Hernquist( r, Cloud_R0, Cloud_Rho0 );
   else if ( Cloud_Model == CLOUD_MODEL_EINASTO   )   dens = AnalyticalDensProf_Einasto  ( r, Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor );
   else   Aux_Error( ERROR_INFO, "Unsupported Cloud_Model = %d !!\n", Cloud_Model );

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
#     ifndef SUPPORT_GSL
      Aux_Error( ERROR_INFO, "Must enable SUPPORT_GSL for integration of enclosed mass in Par_EquilibriumIC !!\n" );
#     endif

#     ifdef SUPPORT_GSL
      double abs_error;

//    Arguments for the gsl integration
      const double lower_limit = 0.0;
      const double upper_limit = r;
      const double abs_err_lim = 1.0e-4;
      const double rel_err_lim = 1.0e-4;
      const int    integ_size  = 1000;
      const int    integ_rule  = 1;    // 1 = GSL_INTEG_GAUSS15, the 15 point Gauss-Kronrod rule

      gsl_integration_workspace * w = gsl_integration_workspace_alloc( integ_size );

      gsl_function F;

//    Parameters for the integrand
      struct mass_integrand_params_Einasto integrand_params_Einasto = { Cloud_R0, Cloud_Rho0, Cloud_Einasto_Power_Factor };
      struct mass_integrand_params_Table   integrand_params_Table   = { InputTable_DensProf_NBin, InputTable_DensProf_Radius, InputTable_DensProf_Density };

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
// Function    :  getExternalPotential
// Description :  Get the external potential at radius r
//
// Note        :  1. Get the external potential from table interpolation or the analytical potential profile
//
// Parameter   :  r : radius
//
// Return      :  external potential at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getExternalPotential( const double r )
{
   double ext_pot;

   if      ( AddExtPot_Table      )   ext_pot = ExtendedInterpolatedTable( r, InputTable_ExtPot_NBin, InputTable_ExtPot_Radius, InputTable_ExtPot_Potential );
   else if ( AddExtPot_Analytical )   ext_pot = UserDefAnlaytical_ExtPot( r );
   else                               ext_pot = 0.0;

   return ext_pot;

} // FUNCTION : getExternalPotential



//-------------------------------------------------------------------------------------------------------
// Function    :  getRandomSampleVelocity
// Description :  Get the ramdomly sampled magnitude of the velocity of a particle at radius r from the distribution function
//
// Note        :  1. The probability of the magnitude of the velocity with a given radius, p(v|r) \propto v^2*f(E),
//                   where f(E) = DFunc(E) = DFunc(Psi(r)-1/2 v^2) is the distribution function
//                2. Psi = E + 1/2 v^2 -> v = \sqrt{ 2*( Psi - E ) }
//                3. dv = - \frac{ 1 }{ \sqrt{ 2*( Psi - E ) } } dE
//                4. The state E<0 is unbound, v^2 > v_{esc}^2
//                5. E is in the range  0 <= E <= Psi(r)
//
// Parameter   :  r : radius
//
// Return      :  Magnitude of velocity for a particle at radius r
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getRandomSampleVelocity( const double r )
{
// The relative potential at this radius
   const double Psi = -ExtendedInterpolatedTable( r, RNPoints, RArray_R, RArray_Phi );

// Cumulative Probability = \int_{v}^{v_max} v^2 f(E) dv
//                        = \int_{E_min}^{E} (2*(Psi-E)) f(E) \frac{ 1 }{ \sqrt{ 2*( Psi - E )} } dE
//                        = \int_{E_min}^{E} \sqrt{ (2*(Psi-E)) } f(E) dE
//                        where Psi - 1/2 v_max^2 = E_min
   double Probability;

   CumulProbaDistr_GivenRadius[0] = 0.0;
   for (int i=1; i<ENPoints; i++)
   {
      if ( EArray_E[i] > Psi )   Probability = ( EArray_E[i-1] > Psi ) ? 0.0 : 0.5*EArray_DFunc[i-1]*sqrt(2.0*(Psi-EArray_E[i-1]))*(Psi-EArray_E[i-1]);
      else                       Probability = 0.5*( EArray_DFunc[i-1]*sqrt(2.0*(Psi-EArray_E[i-1])) +
                                                     EArray_DFunc[i  ]*sqrt(2.0*(Psi-EArray_E[i  ])) )*EArray_dE;

      CumulProbaDistr_GivenRadius[i] = CumulProbaDistr_GivenRadius[i-1] + Probability;
   }

   const double TotalProbability        = CumulProbaDistr_GivenRadius[ELastIdx];
   const double RandomSampleProbability = TotalProbability*Random_Num_Gen->GetValue( 0, 0.0, 1.0 );
   const double RandomSampleE           = Mis_InterpolateFromTable( ENPoints, CumulProbaDistr_GivenRadius, EArray_E, RandomSampleProbability );

// v^2 = 2*(Psi-E)
   const double RandomSampleVelocity = ( RandomSampleE > Psi ) ? 0.0 : sqrt( 2.0*(Psi-RandomSampleE) );

   return RandomSampleVelocity;

} // FUNCTION : getRandomSampleVelocity



//-------------------------------------------------------------------------------------------------------
// Function    :  getIntegratedDistributionFunction
// Description :  Get the integrated distribution function
//
// Note        :  1. The integrated deistribution function
//                   = \frac{1}{\sqrt{8}\pi^2} \int_{0}^{\mathcal{E}} \frac{ 1 }{ \sqrt{ \mathcal{E} -\Psi } } \frac{ d\rho }{ d\Psi } d\Psi
//                   = \frac{1}{\sqrt{8}\pi^2} \int_{0}^{\mathcal{E}} -2 \frac{ d\rho }{ d\Psi } d\sqrt{ \mathcal{E} -\Psi } }
//
// Parameter   :  E  : Ralative energy, \mathcal{E} in the above equation
//
// Return      :  Integrated distribution function
//-------------------------------------------------------------------------------------------------------
double Par_EquilibriumIC::getIntegratedDistributionFunction( const double E )
{
   const int    N_intervals = MIN( (int)ceil(E/EArray_dE), 1000 );
   const double dPsi     = E/N_intervals;

   double dRho_dPsi = 0.0;
   double integral  = 0.0;
   for (int i=0; i<N_intervals; i++)
   {
      const double Psi             = (i+0.5)*dPsi;
      const double dsqrt_EminusPsi = ( ( (Psi+0.5*dPsi) >= E ) ? 0.0 : sqrt( E-(Psi+0.5*dPsi) ) ) - sqrt( E-(Psi-0.5*dPsi) );

      if ( Psi > -RArray_Phi[RLastIdx] )
      {
         dRho_dPsi = Mis_InterpolateFromTable( RNPoints, RArray_Phi, RArray_dRho_dPsi, -(Psi) );
      }
      else
      {
         // dRho_dPsi is extrapolated with a power law for low Psi
         double power_extrapolated = log( RArray_dRho_dPsi[RLastIdx-1]/RArray_dRho_dPsi[RLastIdx] )/log( (-RArray_Phi[RLastIdx-1])/(-RArray_Phi[RLastIdx]) );

         // to avoid the slope of distribution function become negative when the power of Psi is less then 1/2
         power_extrapolated = ( power_extrapolated < 0.5 ) ? 0.5 : power_extrapolated;

         const double dRho_dPsi_Extrapolated = RArray_dRho_dPsi[RLastIdx]*pow( Psi/(-RArray_Phi[RLastIdx]), power_extrapolated );

         dRho_dPsi = ( dRho_dPsi_Extrapolated > 0 ) ? dRho_dPsi_Extrapolated : 0.0;
      }

      integral += -2.0*dRho_dPsi*dsqrt_EminusPsi;
   }

   return integral/(sqrt(8.0)*SQR(M_PI));

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

   return M_PI*dr*( r0*r0*( 6.0*rho0 + 6.0*rho1 ) + r0*dr*( 4.0*rho0 + 8.0*rho1 ) + dr*dr*( rho0 + 3.0*rho1 ) )/3.0;

} // FUNCTION : LinearDensityShellMass



//-------------------------------------------------------------------------------------------------------
// Function    :  UserDefAnlaytical_ExtPot
// Description :  User-defined analytical external potential formula
//
// Note        :  1. As an example, there is a Plummer potential with hardcoded parameters
//
// Parameter   :  r : radius
//
// Return      :  external potential at radius r
//-------------------------------------------------------------------------------------------------------
double UserDefAnlaytical_ExtPot( const double r )
{
// DEFINE YOUR FUNCTION HERE !!!
/*
   return -NEWTON_G*M/r;
*/

// Or use one of the built-in models
/*
   return AnalyticalPoteProf_NFW( r, 0.1, 1.0 );
   return AnalyticalPoteProf_Jaffe( r, 0.1, 1.0 );
   return AnalyticalPoteProf_Hernquist( r, 0.1, 1.0 );
*/

   return AnalyticalPoteProf_Plummer( r, 0.05, 80.0 );

} // FUNCTION : UserDefAnlaytical_ExtPot



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

   return Rho0*pow( 1.0+x*x, -2.5 );

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

   return (4.0/3.0)*M_PI*Rho0*CUBE(r)*pow( 1.0+x*x, -1.5 );

} // FUNCTION : AnalyticalMassProf_Plummer



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalPoteProf_Plummer
// Description :  Analytical gravitational potential of the Plummer model
//
// Note        :  1. \Phi(r) = \frac{ -G*M_0 }{ ( r^2 + a^2 )^{1/2} },
//                    where M_0 = \frac{4\pi}{3} a^3 \rho_0
//
// Parameter   :  r    : input radius
//                R0   : Plummer scale radius, a
//                Rho0 : Plummer scale density, \rho_0
//
// Return      :  gravitational potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalPoteProf_Plummer( const double r, const double R0, const double Rho0  )
{
   return -NEWTON_G*(4.0/3.0)*M_PI*CUBE(R0)*Rho0/sqrt( SQR(r) + SQR(R0) );

} // FUNCTION : AnalyticalPoteProf_Plummer



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

   return Rho0/( x*SQR( 1.0+x ) );

} // FUNCTION : AnalyticalDensProf_NFW



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_NFW
// Description :  Analytical enclosed mass profile of the NFW model
//
// Note        :  1. M(r) = 4\pi \rho_s r_s^3 [ \ln( \frac{ r_s + r }{ r_s } ) - \frac{ r }{ r_s + r } ]
//                2. Reference: Navarro J.~F., Frenk C.~S., White S.~D.~M., 1996, ApJ, doi:10.1086/177173
//                              Li P. et al., 2020, ApJS, doi:10.3847/1538-4365/ab700e
//                              Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.66)
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

   return 4.0*M_PI*Rho0*CUBE(R0)*( log( 1.0+x ) - x/(1.0+x) );

} // FUNCTION : AnalyticalMassProf_NFW



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalPoteProf_NFW
// Description :  Analytical gravitational potential of the NFW model
//
// Note        :  1. \Phi(r) = -4\pi G\rho_s r_s^2 \frac{\ln(1+r/r_s)}{r/r_s}
//                2. Reference: Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.67)
//
// Parameter   :  r    : input radius
//                R0   : NFW scale radius, r_s
//                Rho0 : NFW scale density, \rho_s
//
// Return      :  gravitational potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalPoteProf_NFW( const double r, const double R0, const double Rho0  )
{
   const double x = r/R0;

   return -4.0*M_PI*NEWTON_G*Rho0*SQR(R0)*( log( 1.0+x )/x );

} // FUNCTION : AnalyticalPoteProf_NFW



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

   return Rho0/( (1.0+x)*(1.0+x*x) );

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

   return 2.0*M_PI*Rho0*CUBE(R0)*( 0.5*log( 1.0+x*x ) + log( 1.0+x ) - atan( x ) );

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

   return Rho0/( 4.0*M_PI*SQR(x)*SQR(1.0+x) ); // return Rho0/(x*(1.0+x)); //previous one

} // FUNCTION : AnalyticalDensProf_Jaffe



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Jaffe
// Description :  Analytical enclosed mass profile of the Jaffe model
//
// Note        :  1. M(r) = M_J \frac{ r }{ r_J + r }
//                   ,where M_J = \rho_0*r_J^3 is the total mass
//                2. Reference: Jaffe W., 1983, MNRAS, doi:10.1093/mnras/202.4.995
//                              Ciotti L. and Ziaee Lorzad A., 2018, MNRAS, doi:10.1093/mnras/stx2771
//                              Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.66)
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

   return Rho0*CUBE(R0)*x/(1.0+x);

} // FUNCTION : AnalyticalMassProf_Jaffe



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalPoteProf_Jaffe
// Description :  Analytical gravitational potential of the Jaffe model
//
// Note        :  1. \Phi(r) = -G\rho_0 r_J^2 \ln(1+r_J/r)
//                2. Reference: Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.67)
//
// Parameter   :  r    : input radius
//                R0   : Jaffe scale radius, r_J
//                Rho0 : Jaffe scale density, \rho_0
//
// Return      :  gravitational potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalPoteProf_Jaffe( const double r, const double R0, const double Rho0  )
{
   return -NEWTON_G*Rho0*SQR(R0)*log( 1.0 + R0/r );

} // FUNCTION : AnalyticalPoteProf_Jaffe



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

   return Rho0/( x*CUBE( 1.0+x ) );

} // FUNCTION : AnalyticalDensProf_Hernquist



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalMassProf_Hernquist
// Description :  Analytical enclosed mass profile of the Hernquist model
//
// Note        :  1. M(r) = M_0 \frac{ r^2 }{ (r+a)^2 }
//                   ,where M_0 = 2\pi*\rho_0*a^3 is the total mass
//                2. Reference: Hernquist L., 1990, ApJ, doi:10.1086/168845
//                              Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.66)
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

   return 2.0*M_PI*Rho0*CUBE(R0)*SQR(x)/SQR(1.0+x);

} // FUNCTION : AnalyticalMassProf_Hernquist



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalPoteProf_Hernquist
// Description :  Analytical gravitational potential of the Hernquist model
//
// Note        :  1. \Phi(r) = -4\pi G\rho_0 a^2 \frac{1}{2(1+r/a)}
//                2. Reference: Binney J. & Tremaine S., 2008, Galactic Dynamics (2nd ed.), Eq(2.67)
//
// Parameter   :  r    : input radius
//                R0   : Hernquist scale radius, a
//                Rho0 : Hernquist scale density, \rho_0
//
// Return      :  gravitational potential at the given radius
//-------------------------------------------------------------------------------------------------------
double AnalyticalPoteProf_Hernquist( const double r, const double R0, const double Rho0  )
{
   const double x = r/R0;

   return -2.0*M_PI*NEWTON_G*Rho0*SQR(R0)/(1.0+x);

} // FUNCTION : AnalyticalPoteProf_Hernquist



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

   return 4.0*M_PI*SQR(r)*AnalyticalDensProf_Einasto( r, R0, Rho0, Einasto_Power_Factor );

} // FUNCTION : MassIntegrand_Einasto



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

   return 4.0*M_PI*SQR(r)*ExtendedInterpolatedTable( r, NBin, Table_R, Table_D );

} // FUNCTION : MassIntegrand_Table



#endif // #ifdef MASSIVE_PARTICLE
