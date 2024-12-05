#include "../include/General.h"



void ijk2xyz( int i,  int j,  int k,  real *x,  real *y,  real *z, int numCellXYZ[], real BoxSize[] );
double Synchrotron_emissivity( real CREngy, real B, real spectral_index, real observedFreq, real gamma_max, real gamma_min );
double BField( real x,  real y, real z);



//-------------------------------------------------------------------------------------------------------
// Function    :  Synchrotron_3D
// Description :
// Note        :
// Parameter   :  CREngy      : CR energy density                    (eV/cm**3)
//                numCellXYZ  : number of cells on each side of box
//                BoxSize     : simulation box size                  (kpc)
//
// Return      :  Synchrotron : synchrotron emissivity               (erg/cm**3/s)
//-------------------------------------------------------------------------------------------------------
void Synchrotron_3D( real ***CREngy, real spectral_index, real ***Synchrotron, real observedFreq,
                     int numCellXYZ[], real BoxSize[], real gamma_max, real gamma_min )
{
   for (int i=0; i<numCellXYZ[0]; i++)
   {
      for (int j=0; j<numCellXYZ[1]; j++)
      {
#        pragma omp parallel for  num_threads(NUM_THREADS)
         for (int k=0; k<numCellXYZ[2]; k++)
         {
            real x, y, z;
            double B;
            ijk2xyz( i, j, k, &x, &y, &z, numCellXYZ, BoxSize );
            B = BField(x, y, z);
            Synchrotron[i][j][k] = Synchrotron_emissivity( CREngy[i][j][k], B, spectral_index, observedFreq, gamma_max, gamma_min );
         } // for (int k=0; k<numCellXYZ[2]; k++)
      } // for (int j=0; j<numCellXYZ[1]; j++)
   } // for (int i=0; i<numCellXYZ[0]; i++)
} // FUNCTION : Synchrotron_3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Synchrotron_emissivity
// Description :
// Note        :
// Parameter   :  CREngy          : CR energy density                    (eV/cm**3)
//                B               : magnetic field                       (Gauss)
//                spectral_index  : the spectral index of CR
//                observedFreq    : observed frequency                   (Hz)
//
// Return      :  Emissivity      : synchrotron emissivity               (erg/cm**3/s)
//-------------------------------------------------------------------------------------------------------
double Synchrotron_emissivity( real CREngy, real B, real spectral_index, real observedFreq, real gamma_max, real gamma_min )
{
  if ( CREngy*EV2ERG < 1e-23 )   return 0.0;

  double NomalizedConst;

  if ( spectral_index == 2.0 )
  {
     NomalizedConst = CREngy / log(gamma_max/gamma_min) / (ELECTRON_MASS_ENERGY);
  }
  else
  {
     NomalizedConst  = CREngy * ( 2.0 - spectral_index );
     NomalizedConst /= pow( gamma_max, 2.0 - spectral_index ) - pow( gamma_min, 2.0 - spectral_index );
     NomalizedConst /= ELECTRON_MASS_ENERGY;
  }

  double a, emissivity;
  a  =  pow( 2.0, 0.5*(spectral_index - 1.0) ) * 1.732;
  a *=  tgamma( (3.0*spectral_index-1.0)/12.0 );
  a *=  tgamma( (3.0*spectral_index+19.0)/12.0 );
  a *=  tgamma( (    spectral_index+5.0)/ 4.0 );
  a /=  14.179630799142833 * ( spectral_index+1.0 );
  a /= tgamma( (     spectral_index+7.0)/ 4.0 );


  emissivity  = 1e40;  // initialize a huge number to prevent from under flow.
  emissivity *= 4.0 * M_PI * pow( B, 0.5*( 1.0 + spectral_index ));
  emissivity /= ELECTRON_MASS_ENERGY*EV2ERG;
  emissivity *= pow( 3.0*ELETRON_CHARGE / ( 4.0*M_PI*ELECTRON_MASS*SPEED_OF_LIGHT), 0.5*(spectral_index-1.0) );
  emissivity *= pow( observedFreq, -0.5*(spectral_index-1.0) );
  emissivity *= NomalizedConst;
  emissivity *= CUBE( ELETRON_CHARGE );
  emissivity *= a;

  //if ( emissivity != 0.0 ) printf("emissivity=%e\n", emissivity);
  return emissivity;
} // FUNCTION : Synchrotron_emissivity



//-------------------------------------------------------------------------------------------------------
// Function    :  BField
// Description :
// Note        :
// Parameter   :  x: (kpc)
//                y: (kpc)
//                z: (kpc)
// Return      :
//-------------------------------------------------------------------------------------------------------
double BField( real x, real y, real z )
{
  const double R0 = 10.0; // kpc
  const double z0 = 2.0;  // kpc
  const double B0 = 50e-6; // Gauss

  double R = SQRT( x*x + y*y );

  return B0 * exp( -fabs(z)/z0 ) * exp( -R/R0 );
} // FUNCTION : BField



//-------------------------------------------------------------------------------------------------------
// Function    :  ijk2xyz
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void ijk2xyz( int i, int j, int k, real *x, real *y, real *z, int numCellXYZ[], real BoxSize[] )
{
  *x  = (BoxSize[0]*(real)i) / (real)numCellXYZ[0] - BoxSize[0] * 0.5;
  *y  = (BoxSize[1]*(real)j) / (real)numCellXYZ[1] - BoxSize[1] * 0.5;
  *z  = (BoxSize[2]*(real)k) / (real)numCellXYZ[2] - BoxSize[2] * 0.5;
} // FUNCTION : ijk2xyz
