#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>
#include <CCfits/CCfits>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "../include/General.h"


#define kTwoPi  (2.0*M_PI)
#define kSpeedOfLight_SI  299792458.0 // m/s
#define kConvertRadiansToDegrees  (180.0/M_PI)



struct Parameters_t
{
   double *energyBins;
   double *numDensPerEnergyArray;
   int     energyBinsSize;

   double  p_plus_2;
   double  scatteredEnergy;
   double  photonEnergy;
   double  gamma_max;
   double (*fun_ptr)( double, double, double );
};

real GammaRay_Leptonic_Emissivity( double energyBins[], double numDensPerEnergyArray[],
                                   int energyBinsSize, real CREngy, real scatteredEnergy,
                                   real gamma_max, real gamma_min, real spectral_index );

void ReadInModel( std::string, std::vector<double>&, std::vector<double>&, std::vector<double>&,
                  std::vector<double>&, std::vector<std::unique_ptr<std::valarray<double>>>& );
void Interpolate3DRadiationField( std::array<double, 3>, std::valarray<double>&, const std::vector<double>&, const std::vector<double>&,
                                  const std::vector<double>&, const std::vector<std::unique_ptr<std::valarray<double>>>& );

double gsl_integrandInner( double x, void *params );
double gsl_integrandOuter( double epsilon, void *params );
double integralInner( double lb, double ub, void *params );
double integralOuter( double lb, double ub, struct Parameters_t *parameters );
double cmb_differetial_number_density( double frequency );
double distributionFun( const double gamma, const double photonEnergy, const double scatteredEnergy );
double Bisection( double xmin, double xmax, double (*fun_ptr)(double,double,double), double tolerence, void *params );
void FindLbUb( double xmin, double xmax, double (*fun)(double,double,double), void *params, double *lb, double *ub );

// gsl arrays (reuse the memory for better performance)
static gsl_integration_workspace *workspaces_inner[NUM_THREADS];
static gsl_integration_workspace *workspaces_outer[NUM_THREADS];
static gsl_interp_accel          *accels          [NUM_THREADS];
static gsl_spline                *splines         [NUM_THREADS];



//-------------------------------------------------------------------------------------------------------
// Function    :  initialize_gsl_arrays
// Description :  Allocate gsl arrays
//
// Note        :  Each MPI_Rank and each openMP thread should have their own array.
//
// Parameter   :  energyBinSize :
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
static void initialize_gsl_arrays( const int energyBinsSize )
{
   for (int w=0; w<NUM_THREADS; w++)
   {
      workspaces_inner[w] = gsl_integration_workspace_alloc( 1000 );
      if ( workspaces_inner[w] == NULL )
      {
         printf( "Failed to allocate workspace %d\n", w );
         exit( -1 );
      }

      workspaces_outer[w] = gsl_integration_workspace_alloc( 1000 );
      if ( workspaces_outer[w] == NULL )
      {
         printf( "Failed to allocate workspace %d\n", w );
         exit( -1 );
      }

      accels[w] = gsl_interp_accel_alloc();
      if ( accels[w] == NULL )
      {
         printf( "Failed to allocate accelerate space %d\n", w );
         exit( -1 );
      }

      splines[w] = gsl_spline_alloc( gsl_interp_cspline, energyBinsSize );
      if ( splines[w] == NULL )
      {
         printf( "Failed to allocate spline space %d\n", w );
         exit( -1 );
      }
   } // for (int w=0; w<NUM_THREADS; w++)
} // FUNCTION : initialize_workspace



//-------------------------------------------------------------------------------------------------------
// Function    :  free_gsl_arrays
// Description :  Free all the allocated gsl arrays
//
// Parameter   :  none
//
// Return      :  none
//-------------------------------------------------------------------------------------------------------
static void free_gsl_arrays()
{
   for (int w=0; w<NUM_THREADS; w++)
   {
      gsl_integration_workspace_free( workspaces_inner[w] );
      gsl_integration_workspace_free( workspaces_outer[w] );
      gsl_interp_accel_free( accels[w] );
      gsl_spline_free( splines[w] );
   } // for (int w=0; w<NUM_THREADS; w++)
} // FUNCTION : free_workspaces



//-------------------------------------------------------------------------------------------------------
// Function    :  Leptonic_3D
// Description :
//
// Note        :
//
// Parameter   :  scatteredEnergy : scattered photon energy              (eV)
//                CREngy          : CR energy density                    (eV/cm**3)
//                gamma_max       : the maximum Lorentz factor of CR
//                gamma_min       : the minimum Lorentz factor of CR
//                spectral_index  : the spectral index of CR
//                numCellXYZ      : number of cells on each side of box
//                BoxSize         : simulation box size
//                **argv          : the name of ISRF table
//
// Return      :  Emissivity      : gamma-ray emissivity                 (1/cm**3/s)
//-------------------------------------------------------------------------------------------------------
void Leptonic_3D( real scatteredEnergy, real ***CREngy, real gamma_max, real gamma_min, real spectral_index,
                  real ***Emissivity, int numCellXYZ[], real BoxSize[], char **argv )
{

   const std::string indexfile = argv[5]; // Filename with .dat suffix

   std::vector<double> r, phi, z, frequency;

   std::vector<std::unique_ptr<std::valarray<double>>> numdenlist;

   MASTER_PRINT( "Reading files ...\n" );

   ReadInModel( indexfile, r, phi, z, frequency, numdenlist );

   MPI_Barrier( MPI_COMM_WORLD );

   MASTER_PRINT( "Reading files ... done\n" );

   int energyBinsSize = frequency.size();

   double *energyBins = (double*)malloc( energyBinsSize*sizeof(double) );

   for (int i=0; i<energyBinsSize; i++)   energyBins[i] = frequency[i]*PLANCK_EV; // eV

// compute frequency difference
   double *diff_frequency = (double *)malloc( energyBinsSize*sizeof(double) );

   for (int i=2; i<energyBinsSize; i++)   diff_frequency[i-1] = 0.5*( frequency[i] - frequency[i-2] );
   diff_frequency[               0] = frequency[               1] - frequency[               0];
   diff_frequency[energyBinsSize-1] = frequency[energyBinsSize-1] - frequency[energyBinsSize-2];

// 0. initialize the gsl arrays
   MASTER_PRINT( "Initializing gsl arrays ...\n" );
   initialize_gsl_arrays( energyBinsSize );
   MASTER_PRINT( "Initializing gsl arrays ... Done\n" );

// 1. define the counts_i
   int counts_i[NRank];
   for (int c=0; c<NRank-1; c++)   counts_i[c] = ( numCellXYZ[0]-numCellXYZ[0]%NRank ) / NRank;

   counts_i[NRank-1] = numCellXYZ[0]-(NRank-1) * counts_i[0];

   if ( MPI_Rank == 1 )
   for (int c=0; c<NRank; c++)   printf( "counts_i[%d]=%d\n", c, counts_i[c] );

// 2. define the counts
   int counts[NRank];
   for (int c=0; c<NRank; c++)
   {
      checkInt32Overflow( numCellXYZ[1],               numCellXYZ[2], '*', __LINE__ );
      checkInt32Overflow( numCellXYZ[1]*numCellXYZ[2], counts_i[c],   '*', __LINE__ );
      counts[c] = numCellXYZ[1] * numCellXYZ[2] * counts_i[c];
   }

// 3. define dis_i
   int dis_i[NRank];
   for (int c=0; c<NRank; c++)   dis_i[c] = c*counts_i[c];


// 4. calculate emissivity
   MASTER_PRINT( "Calculating emissivity ...\n" );
   for (int ii=0; ii<counts_i[MPI_Rank]; ii++)
   {
#     pragma omp parallel for collapse(2) num_threads( NUM_THREADS ) schedule(runtime)
      for (int j=0; j<numCellXYZ[1]; j++)
      {
         //printf("%d/%d\n", ii+dis_i[MPI_Rank], j);
         //fflush(stdout);
// #        pragma omp parallel for num_threads( NUM_THREADS )
         for (int k=0; k<numCellXYZ[2]; k++)
         {
             int i = ii+dis_i[MPI_Rank];

             double Pos[3];
             int Index[3] = { i, j, k };

             Index2Position( Index, Pos, numCellXYZ, BoxSize );

             // interpolate energy density (eV/cm**3)
             std::valarray<double> nd(0., frequency.size());
             std::array<double, 3> Pos_arr = { Pos[0], Pos[1], Pos[2]};
             Interpolate3DRadiationField( Pos_arr, nd, r, phi, z, numdenlist );

             //// compute energy difference
             //double diff_nd[energyBinsSize];

             //for (int ndi=2; ndi<energyBinsSize; ndi++) diff_nd[ndi-1] = 0.5*( nd[ndi] - nd[ndi-2] );
             //diff_nd[               0] = nd[               1] - nd[               0];
             //diff_nd[energyBinsSize-1] = nd[energyBinsSize-1] - nd[energyBinsSize-2];

             //for (int ndi=0; ndi<energyBinsSize; ndi++){
             //  if (diff_nd[ndi] == 0.0) diff_nd[ndi] = 4e-30;
             //}

//           convert energy density (eV cm**-3) to differential number density ( 1 eV**-1 cm**-3 )
             double numDensPerEnergyArray[energyBinsSize];

             for (int f=0; f<energyBinsSize-1; f++)
             {
               //numDensPerEnergyArray[f]  = diff_nd[f]/diff_frequency[f]/energyBins[f]/PLANCK_EV;
               numDensPerEnergyArray[f]  = nd[f] / diff_frequency[f] / energyBins[f] / PLANCK_EV;
               numDensPerEnergyArray[f] += cmb_differetial_number_density( frequency[f] );
             }

             Emissivity[i][j][k] = GammaRay_Leptonic_Emissivity( energyBins, numDensPerEnergyArray, energyBinsSize,
                                                                 CREngy[i][j][k], scatteredEnergy, gamma_max, gamma_min,
                                                                 spectral_index );
         } // for (int k=0; k<numCellXYZ[2]; k++)
      } // for (int j=0; j<numCellXYZ[1]; j++)
   } // for (int ii=0; ii<counts_i[MPI_Rank]; ii++)

   MASTER_PRINT( "Calculating emissivity ... Done\n" );

   MPI_Barrier( MPI_COMM_WORLD );

// 5. free gsl arrays
   MASTER_PRINT( "Freeing gsl arrays ...\n" );
   free_gsl_arrays();
   MASTER_PRINT( "Freeing gsl arrays ... Done\n" );

   MASTER_PRINT( "Before sending ...\n\n" );

   MPI_Barrier( MPI_COMM_WORLD );

   // Declare the buffer and attach it
   int buffer_attached_size = MPI_BSEND_OVERHEAD + counts[MPI_Rank]*sizeof(real);
   char* buffer_attached = (char*)malloc( buffer_attached_size );
   MPI_Buffer_attach( buffer_attached, buffer_attached_size );

   MPI_Bsend( &Emissivity[dis_i[MPI_Rank]][0][0], counts[MPI_Rank], MPI_MYREAL, RootRank, MPI_Rank, MPI_COMM_WORLD );

   // receive Emissivity[][][] from all ranks
   if ( MPI_Rank == RootRank )
   {
     MPI_Status MPI_status;

     for ( int rank=0; rank<NRank; rank++ )
     {
       MPI_Recv( &Emissivity[dis_i[rank]][0][0], counts[rank], MPI_MYREAL, rank, rank, MPI_COMM_WORLD, &MPI_status );
     }
   } // if ( MPI_Rank == RootRank )

   MPI_Barrier( MPI_COMM_WORLD );

   // Detach the buffer. It blocks until all messages stored are sent.
   MPI_Buffer_detach( &buffer_attached, &buffer_attached_size );

   free( buffer_attached );
   free( energyBins );
   free( diff_frequency );

   //if (MPI_Rank == RootRank) OutputBinary(&Emissivity[0][0][0], sizeof(real), numCellXYZ[0]*numCellXYZ[1]*numCellXYZ[2], "leptonic.dat");
} // FUNCTION : Leptonic_3D



//-------------------------------------------------------------------------------------------------------
// Function    :  GammaRay_Leptonic_Emissivity
// Description :
// Note        :
// Parameter   :  energyBins            : 2.734140e+10*PLANCK_EV - 3.287159e+15*PLANCK_EV (eV)
//                numDensPerEnergyArray : differential number density                     ( 1 eV**-1 cm**-3 )
//                energyBinsSize        : 128
//                CREngy                : CR energy density                               (eV/cm**3)
//                scatteredEnergy       : scattered photon energy                         (eV)
//                gamma_max             : the maximum Lorentz factor of CR
//                gamma_min             : the minimum Lorentz factor of CR
//                spectral_index        : the spectral index of CR
// Return      :  Emissivity            : gamma-ray emissivity                            (1/cm**3/s)
//-------------------------------------------------------------------------------------------------------
real GammaRay_Leptonic_Emissivity( double energyBins[], double numDensPerEnergyArray[],
                                   int energyBinsSize, real CREngy, real scatteredEnergy,
                                   real gamma_max, real gamma_min, real spectral_index )
{
  struct Parameters_t parameters;

  parameters.energyBins            = energyBins;
  parameters.numDensPerEnergyArray = numDensPerEnergyArray;
  parameters.energyBinsSize        = energyBinsSize;
  parameters.p_plus_2              = spectral_index + 2.0;
  parameters.scatteredEnergy       = scatteredEnergy;
  parameters.photonEnergy          = NAN;
  parameters.gamma_max             = gamma_max;
  parameters.fun_ptr               = &distributionFun;

  if ( CREngy*EV2ERG < 1e-23 )   return 0.0;

  double NomalizedConst;

  if ( spectral_index == 2.0 )
  {
     NomalizedConst = CREngy / log(gamma_max/gamma_min) / ELECTRON_MASS_ENERGY;
  }
  else
  {
     NomalizedConst  = CREngy * ( 2.0 - spectral_index );
     NomalizedConst /= pow( gamma_max, 2.0 - spectral_index ) - pow( gamma_min, 2.0 - spectral_index );
     NomalizedConst /= ELECTRON_MASS_ENERGY;
  }

  real emissivity;
  emissivity  = integralOuter( energyBins[0], energyBins[energyBinsSize-1], &parameters );
  emissivity *= NomalizedConst;
  emissivity *= 0.75 * THOMSON_CROSS_SECTION;
  emissivity *= SPEED_OF_LIGHT;

  return emissivity;
} // FUNCTION : GammaRay_Leptonic_Emissivity



//-------------------------------------------------------------------------------------------------------
// Function    :  Index2Position
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Index2Position( const int Index[], double Pos[], const int numCellXYZ[], const real BoxSize[] )
{
   for (int d=0; d<3; d++)   Pos[d] = ((real)Index[d]+0.5)*BoxSize[d]/numCellXYZ[d] - (real)0.5*BoxSize[d];
} // FUNCTION : Inde2Position



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadInModel
// Description :
// Note        :  A lot of this comes from the RadiationField class in GALPROP for reading
//                the model files into memory
// Parameter   :  numden : energy density (eV/cm^3)
// Return      :
//-------------------------------------------------------------------------------------------------------
void ReadInModel( std::string name,
                  std::vector<double>& rarray,
                  std::vector<double>& phiarray,
                  std::vector<double>& zarray,
                  std::vector<double>& frequency,
                  std::vector<std::unique_ptr<std::valarray<double>>>& numden)
{
//const double kPlanck_SI = 6.62606876e-34; // J s
//const double e_SI  = 1.602176462e-19; // positron charge in coulomb
  const std::string prefix = name.substr( 0, name.find_last_of("/")+1 );

  std::ifstream rf( name.c_str() );

  std::string line;
  std::getline( rf, line );
  const auto pos = line.find_first_of( "!!3D" );

  if ( pos > line.size() )
  {
    std::cerr << "Error reading 3D header from file" << '\n';
    rf.close();
    exit(-1);
  }

  line.clear();
  std::getline( rf, line );
  std::vector<uint64_t> totals;
  { std::stringstream ss( line ); for (;;) { uint64_t val; ss >> val; if ( ss.fail() ) break; totals.push_back( val ); } }
  auto volumeelements( totals[0] ), wlbins( totals[1] ), numfilters( totals[2] ), numrvals( totals[3] ), numphivals( totals[4] ), numzvals( totals[5] );

  line.clear();
  std::getline(rf, line);
  std::vector<double> rdata;
  { std::stringstream ss( line ); for (;;) { double val; ss >> val; if ( ss.fail() ) break; rdata.push_back( val ); } }
  rarray.resize( rdata.size() );
  std::copy( rdata.begin(), rdata.end(), &rarray[0] );

  line.clear();
  std::getline(rf, line);
  std::vector<double> phidata;
  { std::stringstream ss( line ); for (;;) { double val; ss >> val; if ( ss.fail() ) break; phidata.push_back( val ); } }
  phiarray.resize(phidata.size());
  std::copy( phidata.begin(), phidata.end(), &phiarray[0] );

  line.clear();
  std::getline( rf, line );
  std::vector<double> zdata;
  { std::stringstream ss( line ); for (;;) { double val; ss >> val; if ( ss.fail() ) break; zdata.push_back( val ); } }
  zarray.resize( zdata.size() );
  std::copy( zdata.begin(), zdata.end(), &zarray[0] );

  double luminosity, dustmass;
  rf >> luminosity >> dustmass;

  uint64_t stellarcomponents(0);
  rf >> stellarcomponents;

  for (auto ic(0); ic<stellarcomponents; ++ic)
  {
    double lum(0.);
    std::string cname;
    rf >> cname >> lum;
  }

  uint64_t geometry;
  std::array<double, 6> regiondata = {{ 0. }};
  rf >> geometry;
  rf >> regiondata[0] >> regiondata[1] >> regiondata[2] >> regiondata[3] >> regiondata[4] >> regiondata[5];

  for (auto iv(0); iv<volumeelements; ++iv)
  {
    uint64_t idx;
    double x, y, z, dx, dy, dz;
    rf >> idx;
    rf >> x >> y >> z;
    rf >> dx >> dy >> dz;

    std::string filterfilename, filtercountfilename, directfilename, scatteredfilename, transientfilename, thermalfilename, totalfilename, opticalfilename, infraredfilename, filterfluxfilename, fluxfilename;
    rf >> filterfilename;
    rf >> filtercountfilename;
    rf >> directfilename;
    rf >> scatteredfilename;
    rf >> transientfilename;
    rf >> thermalfilename;
    rf >> totalfilename;
    rf >> opticalfilename;
    rf >> infraredfilename;
    rf >> filterfluxfilename;
    rf >> fluxfilename;

    const std::string fullfluxfilename = prefix + fluxfilename;

    CCfits::FITS fluxfile( fullfluxfilename, CCfits::Read );
    CCfits::ExtHDU& table = fluxfile.currentExtension();
    std::vector<double> wl, total, direct, scattered, transient, thermal;
    table.column(1).read( wl,        1, wlbins );
    table.column(2).read( total,     1, wlbins );
    table.column(3).read( direct,    1, wlbins );
    table.column(4).read( scattered, 1, wlbins );
    table.column(5).read( transient, 1, wlbins );
    table.column(6).read( thermal,   1, wlbins );

    if ( frequency.size() == 0 )
    {
       frequency.resize( wlbins, 0. );
       for (size_t iwl(0); iwl<wl.size(); ++iwl)
       frequency[iwl] = kSpeedOfLight_SI / ( wl[wl.size()-1-iwl]*1e-6 ); // wl is read-in with micron units and ascending. Reverse and convert to frequency
    }

    numden.emplace_back( new std::valarray<double>( 0., wlbins ) );
    std::copy( total.rbegin(), total.rend(), &(*numden.back())[0] ); // Copy is reverse order to make ascending with frequency
    //for (auto i(0); i < frequency.size(); ++i)
    //  (*numden.back())[i] *= std::pow(kPlanck_SI/e_SI*frequency[i], -2.); // Convert to eV^-1 cm^-3
  } // for (auto iv(0); iv<volumeelements; ++iv)
} // FUNCTION : ReadInModel



//-------------------------------------------------------------------------------------------------------
// Function    :  Interpolate3DRadiationField
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Interpolate3DRadiationField( std::array<double, 3> coord,
                                  std::valarray<double>& numberdensity,
                                  const std::vector<double>& rarr,
                                  const std::vector<double>& phiarr,
                                  const std::vector<double>& zarr,
                                  const std::vector<std::unique_ptr<std::valarray<double>>>& numdenlist )
{

   const auto r2( coord[0]*coord[0] + coord[1]*coord[1] ), z( coord[2] );
   auto phi = std::atan2( coord[1], coord[0] );
   phi = ( phi < 0. ? phi + kTwoPi : (phi >= kTwoPi ? phi - kTwoPi : phi) ) * kConvertRadiansToDegrees;

   const auto r = ( r2 > rarr[0]*rarr[0] ? std::sqrt(r2) : rarr[0] );

   if ( r >= rarr[rarr.size()-1]  ||  z < zarr[0]  ||  z >= zarr[zarr.size()-1] )   return;

   const auto ridx = std::upper_bound( rarr.begin(), rarr.end(), r ) - rarr.begin();
   const auto dr   = (rarr[ridx] - rarr[ridx-1]);
   assert( dr > 0. );

   const auto rcoeff = ( rarr[ridx] - r ) / dr;
   const auto phiidx = std::upper_bound( phiarr.begin(), phiarr.end(), phi ) - phiarr.begin();
   const auto dphi   = ( phi < phiarr[phiarr.size()-1] ) ? ( phiarr[phiidx] - phiarr[phiidx-1] ) : ( 360. - phiarr[phiarr.size()-1] );
   assert( dphi > 0. );

   const auto phicoeff = ( phi < phiarr[phiarr.size()-1]  &&  dphi > 0. ) ? ( phiarr[phiidx] - phi ) / dphi : ( 360. - phi ) / dphi;
   const auto zidx     = std::upper_bound( zarr.begin(), zarr.end(), z ) - zarr.begin();
   const auto dz       = ( zarr[zidx] - zarr[zidx-1] );
   assert( dz > 0. );

   const auto zcoeff = ( zarr[zidx] - z ) / dz;

   size_t idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8;

   // This deals with the wrap around in phi coordinate between the last bin in the grid and 360./0. deg
   if ( phi > phiarr[phiarr.size()-1] )
   {
      idx1 = (zidx-1)*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx-1);
      idx2 = (zidx-1)*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (0       );
      idx3 = (zidx-1)*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx-1);
      idx4 = (zidx-1)*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (0       );
      idx5 = (zidx  )*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx-1);
      idx6 = (zidx  )*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (0       );
      idx7 = (zidx  )*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx-1);
      idx8 = (zidx  )*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (0       );
   }
   else
   {
      idx1 = (zidx-1)*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx-1);
      idx2 = (zidx-1)*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx  );
      idx3 = (zidx-1)*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx-1);
      idx4 = (zidx-1)*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx  );
      idx5 = (zidx  )*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx-1);
      idx6 = (zidx  )*rarr.size()*phiarr.size() + (ridx-1)*phiarr.size() + (phiidx  );
      idx7 = (zidx  )*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx-1);
      idx8 = (zidx  )*rarr.size()*phiarr.size() + (ridx  )*phiarr.size() + (phiidx  );
   } // if ( phi > phiarr[phiarr.size()-1] ) ... else ...

   const auto& numberdensity1 = *numdenlist[idx1];
   const auto& numberdensity2 = *numdenlist[idx2];
   const auto& numberdensity3 = *numdenlist[idx3];
   const auto& numberdensity4 = *numdenlist[idx4];
   const auto& numberdensity5 = *numdenlist[idx5];
   const auto& numberdensity6 = *numdenlist[idx6];
   const auto& numberdensity7 = *numdenlist[idx7];
   const auto& numberdensity8 = *numdenlist[idx8];

   numberdensity = numberdensity1 *        zcoeff  *        rcoeff   *        phicoeff   +
                   numberdensity2 *        zcoeff  *        rcoeff   * ( 1. - phicoeff ) +
                   numberdensity3 *        zcoeff  * ( 1. - rcoeff ) *        phicoeff   +
                   numberdensity4 *        zcoeff  * ( 1. - rcoeff ) * ( 1. - phicoeff ) +
                   numberdensity5 * ( 1. - zcoeff )*        rcoeff   *        phicoeff   +
                   numberdensity6 * ( 1. - zcoeff )*        rcoeff   * ( 1. - phicoeff ) +
                   numberdensity7 * ( 1. - zcoeff )* ( 1. - rcoeff ) *        phicoeff   +
                   numberdensity8 * ( 1. - zcoeff )* ( 1. - rcoeff ) * ( 1. - phicoeff );
} // FUNCTION : Interpolate3DRadiationField



//-------------------------------------------------------------------------------------------------------
// Function    :  numDensAtEnergy
// Description :
// Note        :
// Parameter   :  photonEnergy          : incident photon energy                          (eV)
//                energyBins            : 2.734140e+10*PLANCK_EV - 3.287159e+15*PLANCK_EV (eV)
//                numDensPerEnergyArray : differential number density                     ( 1 eV**-1 cm**-3 )
//                energyBinsSize        : 128
//
// Return      : numDensAtEnergy        : number density of ISRF at a specific energy     ( 1 eV**-1 cm**-3 )
//-------------------------------------------------------------------------------------------------------
double numDensAtEnergy( double photonEnergy, double energyBins[], double numDensPerEnergyArray[], int energyBinsSize )
{
   const int tid = omp_get_thread_num();
   double numDensAtEnergy;

   gsl_interp_accel *acc    = accels [tid];
   gsl_spline       *spline = splines[tid];

   gsl_spline_init( spline, energyBins, numDensPerEnergyArray, energyBinsSize );

   numDensAtEnergy = gsl_spline_eval( spline, photonEnergy, acc );

   return numDensAtEnergy;
} // FUNCTION : numDensAtEnergy



//-------------------------------------------------------------------------------------------------------
// Function    :  calloc_3d_array
// Description :
// Note        :
// Parameter   :  gamma           : the Lorentz factor of CR
//                photonEnergy    : the energy of incident photon (eV)
//                scatteredEnergy : the energy of scattered photon (eV)
//
// Return      :  f               : distribution function (dimensionless)
//-------------------------------------------------------------------------------------------------------
double distributionFun( const double gamma, const double photonEnergy, const double scatteredEnergy )
{
   const double cr_engy     = gamma * ELECTRON_MASS_ENERGY;
   const double Gamma       = 4.0 * photonEnergy * gamma / ELECTRON_MASS_ENERGY;
   const double q           = 1.0 / ( cr_engy/scatteredEnergy - 1.0 ) / Gamma;
   const double Gamma_q     = Gamma * q;
   const double one_minus_q = 1.0 - q;

   double f;

   f  = 2.0*q*log(q);
   f += ( 1.0+2.0*q )*one_minus_q;
   f += 0.5 * Gamma_q * Gamma_q * one_minus_q / ( 1.0 + Gamma_q );

   return f;
} // FUNCTION : distributionFun



//-------------------------------------------------------------------------------------------------------
// Function    :  gsl_integrandInner
// Description :
// Note        :
// Parameter   :  gamma : the Lorentz factor of CRe
// Return      :
//-------------------------------------------------------------------------------------------------------
double gsl_integrandInner( double gamma, void *params )
{
   const double p_plus_2        = ((struct Parameters_t*)params)->p_plus_2;
   const double scatteredEnergy = ((struct Parameters_t*)params)->scatteredEnergy;
   const double photonEnergy    = ((struct Parameters_t*)params)->photonEnergy;

   return pow( gamma, -p_plus_2 ) * distributionFun( gamma, photonEnergy, scatteredEnergy );
} // FUNCTION : gsl_integrandInner



//-------------------------------------------------------------------------------------------------------
// Function    :  integralInner
// Description :
// Note        :
// Parameter   :  ub: gamma_max
//                lb: gamma_min
// Return      :
//-------------------------------------------------------------------------------------------------------
double integralInner( double lb, double ub, void *params )
{
   const int tid = omp_get_thread_num();
   gsl_integration_workspace *w = workspaces_inner[tid];
   gsl_function F;

   F.function = &gsl_integrandInner;
   F.params   = params;

   double result, error;
   int status, key;

   key = 1;
   gsl_set_error_handler_off();

   do
   {
      //if ( key>1 ) printf( "Unstable! status: %d key: %d (%s: %d)\n", status, key, __FUNCTION__, __LINE__ );
      status = gsl_integration_qag( &F, lb, ub, 0.0, 1e-7, 1000, key, w, &result, &error );
      key++;
   } while ( status  &&  key <= 6 );

   if ( status )
   {
      printf( "Error! status: %d (%s: %d)\n", status, __FUNCTION__, __LINE__ );
#     ifdef DEBUG
      printf( "scatteredEnergy = %20.16e\n", ((struct Parameters_t*)params)->scatteredEnergy );
      printf( "gamma_max       = %20.16e\n", ((struct Parameters_t*)params)->gamma_max       );
      printf( "photonEnergy    = %20.16e\n", ((struct Parameters_t*)params)->photonEnergy    );
      printf( "lb=%20.16e, ub=%20.16e\n", lb, ub );
#     endif
      exit(0);
   } // if ( status )
   return result;
} // FUNCTION : integralInner



//-------------------------------------------------------------------------------------------------------
// Function    :  gsl_integrandOuter
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
double gsl_integrandOuter( double photonEnergy, void *params )
{
   double *numDensPerEnergyArray           = ((struct Parameters_t*)params)->numDensPerEnergyArray;
   double *energyBins                      = ((struct Parameters_t*)params)->energyBins;
   double  scatteredEnergy                 = ((struct Parameters_t*)params)->scatteredEnergy;
   int     energyBinsSize                  = ((struct Parameters_t*)params)->energyBinsSize;
   double  gamma_max                       = ((struct Parameters_t*)params)->gamma_max;
   double (*fun_ptr)(double,double,double) = ((struct Parameters_t*)params)->fun_ptr;

   ((struct Parameters_t*)params)->photonEnergy = photonEnergy;

   double lb, ub;
   double reduced_scatteredEnergy = scatteredEnergy / ELECTRON_MASS_ENERGY;

   lb  = reduced_scatteredEnergy;
   lb += sqrt( SQR( reduced_scatteredEnergy ) + scatteredEnergy/photonEnergy );
   lb *= 0.5;

   ub = gamma_max;

   if ( ub <= lb )   return 0.0;

   return numDensAtEnergy( photonEnergy, energyBins, numDensPerEnergyArray, energyBinsSize ) *
          ( scatteredEnergy / photonEnergy ) *
          integralInner( lb, gamma_max, params );
} // FUNCTION : gsl_integrandOuter



//-------------------------------------------------------------------------------------------------------
// Function    :  integralOuter
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
double integralOuter( double lb, double ub, struct Parameters_t *parameters )
{
   const int tid = omp_get_thread_num();
   gsl_integration_workspace *w  = workspaces_outer[tid];
   gsl_function F;

   F.function = &gsl_integrandOuter;
   F.params = (void*)parameters;

   double result, error;
   int status, key;
   int key_ini = 6;
   key = key_ini;

   gsl_set_error_handler_off();

   do
   {
      //if ( key>key_ini ) printf("Unstable! status: %d key: %d (%s: %d)\n", status, key, __FUNCTION__, __LINE__);
      status = gsl_integration_qag( &F, lb, ub, 0.0, 1e-7, 1000, key, w, &result, &error );
      key++;
   } while ( status  &&  key <= 6 );

   if ( status )
   {
     printf( "Error! status: %d (%s: %d)\n", status, __FUNCTION__, __LINE__ );
#    ifdef DEBUG
     printf( "scatteredEnergy = %20.16e\n", ((struct Parameters_t*)parameters)->scatteredEnergy );
     printf( "gamma_max       = %20.16e\n", ((struct Parameters_t*)parameters)->gamma_max       );
     printf( "photonEnergy    = %20.16e\n", ((struct Parameters_t*)parameters)->photonEnergy    );
     printf( "lb=%20.16e, ub=%20.16e\n", lb, ub );
#    endif
     exit(0);
   } // if ( status )
   return result;
} // FUNCTION : integralOuter



//-------------------------------------------------------------------------------------------------------
// Function    :  cmb_differential_number_density
// Description :
// Note        :
// Parameter   :  frequency: the frequency of emitted photon  (Hz)
// Return      :  n : differential number density      ( 1 eV**-1 cm**-3 )
//-------------------------------------------------------------------------------------------------------
double cmb_differetial_number_density( double frequency )
{
   const double temperature = 2.72548; // in unit K
   const double energy      = REDUCED_PLANCK_EV * 2.0 * M_PI * frequency; // photon energy
   const double kT          = BOLTZMANN_CONST_EV * temperature;
   double n;

   n  = SQR(energy);
   n /= ( EXP( energy/kT ) - 1.0 );
   n /= M_PI* M_PI * pow( REDUCED_PLANCK_EV*SPEED_OF_LIGHT, 3.0 );

   return n;
} // FUNCTION : cmb_differetial_number_density



// Belows are useless
//-------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------------------
// Function    :  Bisection
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
double Bisection( double xmin, double xmax, double (*fun)(double,double,double), double tolerance, void *params )
{
   double scatteredEnergy = ((struct Parameters_t*)params)->scatteredEnergy;
   double photonEnergy    = ((struct Parameters_t*)params)->photonEnergy;

   double x, fun_x;

   const int itrMax = 100;
   int itr = 1;

   while ( itr <= itrMax )
   {
      x = 0.5*( xmin + xmax );

      fun_x = fun( x, photonEnergy, scatteredEnergy );

      if ( fun_x == 0.0  ||  ( 1.0-xmin/xmax < tolerance  &&  fun_x > 0.0  ) )   return x;

      itr++;

      if ( SIGN(fun_x) == SIGN(fun(xmin, photonEnergy, scatteredEnergy)) )
         xmin = x;
      else
         xmax = x;
   } // while ( itr <= itrMax )

   return x;
} // FUNCTION : Bisection



//-------------------------------------------------------------------------------------------------------
// Function    :  FindLbUb
// Description :
// Note        :
// Parameter   :
// Return      :  none
//-------------------------------------------------------------------------------------------------------
void FindLbUb( double xmin, double xmax, double (*fun)(double,double,double), void *params, double *lb, double *ub )
{
   double scatteredEnergy = ((struct Parameters_t*)params)->scatteredEnergy;
   double photonEnergy    = ((struct Parameters_t*)params)->photonEnergy;

   *lb = xmax;

   do
   {
      *ub = *lb;
      *lb = 0.5*( *lb+xmin );
   } while ( fun( *lb, photonEnergy, scatteredEnergy ) >= 0.0 );
} // FUNCTION : FindLbUb
