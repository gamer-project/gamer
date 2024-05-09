#ifndef __PAR_EQUILIBRIUM_IC_H__
#define __PAR_EQUILIBRIUM_IC_H__



#include "GAMER.h"

//gsl library
#ifdef SUPPORT_GSL
#include <gsl/gsl_integration.h>
#endif


// Cloud_Model
#define CLOUD_MODEL_TABLE      0
#define CLOUD_MODEL_PLUMMER    1
#define CLOUD_MODEL_NFW        2
#define CLOUD_MODEL_BURKERT    3
#define CLOUD_MODEL_JAFFE      4
#define CLOUD_MODEL_HERNQUIST  5
#define CLOUD_MODEL_EINASTO    6


class Par_EquilibriumIC
{
   public:
      Par_EquilibriumIC( const char* Type );
      virtual ~Par_EquilibriumIC();

      void setCenter( const double Center_X, const double Center_Y, const double Center_Z );
      void setBulkVel( const double BulkVel_X, const double BulkVel_Y, const double BulkVel_Z );
      void setModelParameters( const double Rho0, const double R0 );
      void setEinastoPowerFactor( const double EinastoPowerFactor );
      void setDensProfTableFilename( const char* DensProfTableFilename );
      void setParticleParameters( const long ParNum, const double MaxR, const int Radial_NBin, const int Energy_NBin, const int RSeed );
      void setExternalPotential( const int AddingExternalPotential, const char* ExtPotTableFilename );
      long getParticleNumber( );

      void initialize();
      void constructParticles( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx );

   private:
      char   Cloud_Type[MAX_STRING];
      int    Cloud_Model                     = -1;
      double Cloud_Center[3]                 = { -1, -1, -1 };
      double Cloud_BulkVel[3]                = {  0,  0,  0 };
      double Cloud_Rho0                      = -1;
      double Cloud_R0                        = -1;
      double Cloud_Einasto_Power_Factor      = -1;
      char   DensProf_Table_Name[MAX_STRING];
      long   Cloud_Par_Num                   = -1;
      double Cloud_MaxR                      = -1;
      int    Cloud_RSeed                     = -1;
      int    AddExtPot                       =  0;
      char   ExtPot_Table_Name[MAX_STRING];

      // Get physical attributes for cloud
      double getEnclosedMass        ( const double r );
      double getDensity             ( const double r );
      double getGraviPotential      ( const double r );
      double getRandomSampleVelocity( const double r );

      double getIntegratedDistributionFunction( const double Psi_Min, const double Psi_Max, const int N_points );

      // Auxiliary functions
      void   RandomVector_GivenLength( const double Length, double RandomVector[3] );

      void   SmoothArray( double* array_x, int index_start, int index_end );

      double ArrayCovariance( const double* array_x, const double* array_y,
                              const int index_start, const int n_elements );

      double Slope_LinearRegression( const double* array_x, const double* array_y,
                                     const int index_start, const int n_elements );

      // Input table of density profile
      void    loadInputDensProfTable();
      int     InputTable_DensProf_nbin = -1;
      double *InputTable_DensProf_radius;
      double *InputTable_DensProf_density;
      double *InputTable_DensProf_enclosedmass;

      // Input table of external potential
      void    loadInputExtPotTable();
      int     InputTable_ExtPot_nbin = -1;
      double *InputTable_ExtPot_radius;
      double *InputTable_ExtPot_potential;

      // Arrays of radial distribution of properties of the cloud
      void    constructRadialArray();
      double  RArray_dR = -1;
      double *RArray_R;
      double *RArray_Rho;
      double *RArray_M_Enc;
      double *RArray_dRho_dR;
      double *RArray_dRho_dPsi;
      double *RArray_G;
      double *RArray_Phi;
      int     RLastIdx    = -1;
      int     RNBin       = -1;

      // Arrays in energy space
      void    constructEnergyArray();
      double  EArray_dE   = -1;
      double  EArray_MinE = -1;
      double  EArray_MaxE = -1;
      double *EArray_E;
      double *EArray_DFunc;
      double *EArray_IntDFunc;
      int     ELastIdx    = -1;
      int     ENBin       = -1;

      // Random number generator
      RandomNumber_t *Random_Num_Gen;

}; // class Par_EquilibriumIC



#endif //__PAR_EQUILIBRIUM_IC_H__
