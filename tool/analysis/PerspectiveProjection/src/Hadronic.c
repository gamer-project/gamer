#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "../include/General.h"



struct Parameters_t
{
  double spectral_index;
  double scatteredEnergy;
};

real crossSection( real Ep );
real F_gamma ( real x, real Ep );
double Integral( double lb, double ub, struct Parameters_t *parameters );


//-------------------------------------------------------------------------------------------------------
// Function    :  Hadronic_3D
// Description :
// Note        :
// Parameter   :  scatteredEnergy : scattered photon energy              (eV)
//                Density         : gas density                          (g/cm**3)
//                CREngy          : CR energy density                    (eV/cm**3)
//                gamma_max       : the maximum Lorentz factor of CR
//                gamma_min       : the minimum Lorentz factor of CR
//                spectral_index  : the spectral index of CR
//                numCellXYZ      : number of cells on each side of box
//                BoxSize         : simulation box size
//
// Return      :  Emissivity      : gamma-ray emissivity                 (eV/cm**3/s)
//-------------------------------------------------------------------------------------------------------
void Hadronic_3D( real scatteredEnergy, real ***Density, real ***CREngy,
                  real gamma_max, real gamma_min, real spectral_index,
                  real ***Emissivity, int numCellXYZ[], real BoxSize[] )
{
// 1. define the counts_i
   int counts_i[NRank];
   for (int c=0; c<NRank-1; c++)  counts_i[c] = ( numCellXYZ[0]-numCellXYZ[0]%NRank ) / NRank;

   counts_i[NRank-1] = numCellXYZ[0] - (NRank-1) * counts_i[0];

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

   for (int ii=0; ii<counts_i[MPI_Rank]; ii++)
   {
      printf( "%d\n", ii+dis_i[MPI_Rank] );
      fflush( stdout );
      for (int j=0; j<numCellXYZ[1]; j++)
      {
#        pragma omp parallel for  num_threads(NUM_THREADS)
         for (int k=0; k<numCellXYZ[2]; k++)
         {
            const int i = ii + dis_i[MPI_Rank];
            Emissivity[i][j][k] = GammaRay_Hadronic_Emissivity( Density[i][j][k]/MASS_PROTON_GRAM, CREngy[i][j][k],
                                                                scatteredEnergy, gamma_max, gamma_min, spectral_index );
         } // for (int k=0; k<numCellXYZ[2]; k++)
      } // for (int j=0; j<numCellXYZ[1]; j++)
   } // for (int ii=0; ii<counts_i[MPI_Rank]; ii++)

   MPI_Barrier( MPI_COMM_WORLD );

// Declare the buffer and attach it
   int buffer_attached_size = MPI_BSEND_OVERHEAD + counts[MPI_Rank]*sizeof(real);
   char* buffer_attached = (char*)malloc(buffer_attached_size);
   MPI_Buffer_attach( buffer_attached, buffer_attached_size );

   MPI_Bsend( &Emissivity[dis_i[MPI_Rank]][0][0], counts[MPI_Rank], MPI_MYREAL, RootRank, MPI_Rank, MPI_COMM_WORLD );


  // receive Emissivity[][][] from all ranks
  if ( MPI_Rank == RootRank )
  {
    MPI_Status MPI_status;
    for (int rank=0; rank<NRank; rank++)
    {
       MPI_Recv( &Emissivity[dis_i[rank]][0][0], counts[rank], MPI_MYREAL, rank, rank, MPI_COMM_WORLD, &MPI_status );
    }
  } // if ( MPI_Rank == RootRank )

  MPI_Barrier( MPI_COMM_WORLD );

  // Detach the buffer. It blocks until all messages stored are sent.
  MPI_Buffer_detach( &buffer_attached, &buffer_attached_size );
  free( buffer_attached );

  //if (MPI_Rank == RootRank) OutputBinary(&Emissivity[0][0][0], sizeof(real), numCellXYZ[0]*numCellXYZ[1]*numCellXYZ[2], "hadronic.dat");
} // FUNCTION : Hadronic_3D



//-------------------------------------------------------------------------------------------------------
// Function    :  GammaRay_Hadronic_Emissivity
// Description :
// Note        :
// Parameter   :  number_density  : proton number density                (1/cm**3)
//                CREngy          : CR energy density                    (eV/cm**3)
//                scatteredEnergy : scattered photon energy              (eV)
//                gamma_max       : the maximum Lorentz factor of CR
//                gamma_min       : the minimum Lorentz factor of CR
//                spectral_index  : the spectral index of CR
//
// Return      :  Emissivity      : gamma-ray emissivity                 (eV/cm**3/s)
//-------------------------------------------------------------------------------------------------------
real GammaRay_Hadronic_Emissivity( real number_density, real CREngy, real scatteredEnergy,
                                   real gamma_max, real gamma_min, real spectral_index )
{
  struct Parameters_t parameters;

  parameters.spectral_index  = spectral_index;
  parameters.scatteredEnergy = scatteredEnergy;

  double NomalizedConst;

  if ( spectral_index == 2.0 )
  {
     NomalizedConst = CREngy / log(gamma_max/gamma_min) / PROTON_MASS_ENERGY;
  }
  else
  {
     NomalizedConst  = CREngy * ( 2.0 - spectral_index );
     NomalizedConst /= pow( gamma_max, 2.0 - spectral_index ) - pow( gamma_min, 2.0 - spectral_index );
     NomalizedConst /= PROTON_MASS_ENERGY;
  }

  real emissivity;

  emissivity  = Integral( 0.0, 1.0, &parameters );
  emissivity *= SPEED_OF_LIGHT * number_density * spectral_index;
  emissivity *= NomalizedConst;
  emissivity *= pow( scatteredEnergy/PROTON_MASS_ENERGY, -spectral_index );

  return emissivity;
} // FUNCTION : GammaRay_Hadronic_Emissivity



//-------------------------------------------------------------------------------------------------------
// Function    :  gsl_integrand
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
double gsl_integrand( double x, void *params )
{
  double spectral_index   = ((struct Parameters_t*)params)->spectral_index;
  double scatteredEnergy  = ((struct Parameters_t*)params)->scatteredEnergy;
  double Ep               = scatteredEnergy / x;

  return crossSection(Ep) * F_gamma(x, Ep) * pow(x, spectral_index);
} // FUNCTION : gsl_integrand



//-------------------------------------------------------------------------------------------------------
// Function    :  Integral
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
double Integral( double lb, double ub, struct Parameters_t *parameters )
{
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
  gsl_function F;

  F.function = &gsl_integrand;
  F.params   = (void*)parameters;

  double result, error;
  int status, key;
  key = 1;

  gsl_set_error_handler_off();

  do
  {
    //if ( key>1 ) printf("Unstable! status: %d key: %d (%s: %d)\n", status, key, __FUNCTION__, __LINE__);

    status = gsl_integration_qag( &F, lb, ub, 0.0, 1e-7, 1000, key, w, &result, &error );

    key++;

  } while( status  &&  key<=6 );

  if ( status ) { printf( "Error! status: %d (%s: %d)\n", status, __FUNCTION__, __LINE__ ); exit(0); }

  gsl_integration_workspace_free( w );

  return result;
} // FUNCTION : Integral



//-------------------------------------------------------------------------------------------------------
// Function    :  F_gamma
// Description :
// Note        :
// Parameter   :  x  : the ratio of scattered photon energy to proton energy, scatteredEnergy/Ep.
//                Ep : the energy of proton (eV)
//
// Return      :  distribution function : Eq. (58) in Kelner (2006) (dimensionless)
//-------------------------------------------------------------------------------------------------------
real F_gamma( real x, real Ep )
{
  real L, F_gamma;
  real B, beta, k;
  real factor;

  L    = log( Ep / 1e12 ); // L = ln( Ep/1 TeV )
  B    = 1.3 + 0.14*L + 0.011*L*L;
  beta = 1.0 / ( 1.79 + 0.11*L + 0.008*L*L );
  k    = 1.0 / ( 0.801 + 0.049*L + 0.014*L*L );

  real pow_x_beta            = pow( x, beta );
  real one_minus_pow_x_beta  = 1.0 - pow_x_beta;
  real k_times_pow_x_beta    = k*pow_x_beta;
  real ln_x                  = log(x);
  real beta_times_pow_x_beta = beta * pow_x_beta;

  factor  = one_minus_pow_x_beta;
  factor /= 1.0 + k_times_pow_x_beta * one_minus_pow_x_beta;
  factor *= factor;
  factor *= factor;

  F_gamma  = 1.0/ln_x;
  F_gamma -= 4.0*beta_times_pow_x_beta/one_minus_pow_x_beta;
  F_gamma -= 4.0*k*beta_times_pow_x_beta*(1.0-2.0*pow_x_beta)/(1.0+k_times_pow_x_beta*one_minus_pow_x_beta);
  F_gamma *= B*ln_x*factor/x;

  return F_gamma;
} // FUNCTION : F_gamma



//-------------------------------------------------------------------------------------------------------
// Function    :  F_gamma
// Description :
// Note        :
// Parameter   : Ep: the energy of incident proton (eV)
//
// Return      : total cross section of pp-collision (cm^-2), Eq.(73) in Kelner (2006)
//-------------------------------------------------------------------------------------------------------
real crossSection( real Ep )
{
  real Eth, factor, L;

  Eth     = 1.22e9; // Eq(79), the threshold energy (eV) of production of neutral \pi meson

  factor  = Eth / Ep;
  factor *= factor;
  factor *= factor;
  factor  = 1.0 - factor;
  factor *= factor;

  L       = log( Ep / 1e12 ); // L = ln( Ep/1 TeV )

  return ( 34.3 + 1.88*L + 0.25*L*L ) * factor * MILLIBARN_2_CM2;
} // FUNCTION : crossSection
