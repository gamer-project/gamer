#ifndef __PAR_EQUILIBRIUM_IC_H__
#define __PAR_EQUILIBRIUM_IC_H__



#include "GAMER.h"
#include "vector"
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

//gsl library
#ifdef SUPPORT_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
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
      void setParticleParameters( const long ParNum, const double MaxR, const int ProfNBin, const int RSeed );
      void setExternalPotential( const int AddingExternalPotential, const char* ExtPotTableFilename );
      long getParticleNumber( );

      void Init();
      void Par_SetEquilibriumIC( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx );


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
      int    Cloud_Rad_NBin                  = -1;
      int    Cloud_RSeed                     = -1;
      int    AddExtPot                       =  0;
      char   ExtPot_Table_Name[MAX_STRING];

      // Get physical attributes for cloud
      double getEnclosedMass        ( const double r );
      double getDensity             ( const double r );
      double getExternalPotential   ( const double r );
      double getGraviPotential      ( const double r );
      double getRandomSampleVelocity( const double r );

      // Initialize arrays in energy space
      void Init_EngArray();

      //  Add External Potential
      void addExternalPotential();

      // Auxiliary functions
      void RandomVector_GivenLength( const double Length, double RandomVector[3] );

      // Solve Eddington's equation
      double getIntegratedDistributionFunction( const double Psi_Min, const double Psi_Max, const int N_points );

      int     Cloud_Eng_NBin             = -1;
      double  EngArray_dBindingEnergy    = -1;
      double  EngArray_MinBindingEnergy  = -1;
      double  EngArray_MaxBindingEnergy  = -1;
      double *EngArray_BindingEnergy;
      double *EngArray_DistriFunc;
      double *EngArray_IntegDistriFunc;
      int     Eng_LastIdx;

      // statistics
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
      double *RadArray_dr;
      double *RadArray_Radius;
      double *RadArray_Density;
      double *RadArray_EnclosedMass;
      double *RadArray_DensitySlope;
      double *RadArray_dRho_dPsi;
      double *RadArray_GraviField;
      double *RadArray_GraviPotential;
      double *RadArray_ExternalPotential;
      int     Rad_LastIdx;

      void setRadArray_Radius();
      void setRadArray_Density();
      void setRadArray_EnclosedMass();
      void setRadArray_DensitySlope();
      void setRadArray_dRho_dPsi();
      void setRadArray_GraviField();
      void setRadArray_GraviPotential();
      void setRadArray_ExternalPotential();

      // Random number generator
      RandomNumber_t *Random_Num_Gen;

}; // class Par_EquilibriumIC



#endif //__PAR_EQUILIBRIUM_IC_H__
