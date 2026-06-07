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
      Par_EquilibriumIC( const char* Cloud_Type );
      virtual ~Par_EquilibriumIC();

//    interface functions to set the parameters from input
      void   setCenterAndBulkVel     ( const double Center_X, const double Center_Y, const double Center_Z,
                                       const double BulkVel_X, const double BulkVel_Y, const double BulkVel_Z );
      void   setModelParameters      ( const double Rho0, const double R0 );
      void   setEinastoPowerFactor   ( const double EinastoPowerFactor );
      void   setDensProfTableFilename( const char* DensProfTableFilename );
      void   setParticleParameters   ( const long ParNum, const double MaxR, const int NBin, const int RSeed );
      void   setExtPotParameters     ( const int AddingExternalPotential_Analytical,
                                       const int AddingExternalPotential_Table, const char* ExtPotTableFilename );

//    main construction process functions to be called
      void   constructDistribution();
      void   constructParticles( real_par *Mass_AllRank, real_par *Pos_AllRank[3], real_par *Vel_AllRank[3], const long Par_Idx );

//    results to be returned
      double TotCloudMass                    = -1;
      double ParticleMass                    = -1;
      double TotCloudMassError               = -1;

   private:
//    parameters for the cloud
      int    Cloud_Model                     = -1;
      double Cloud_Center[3]                 = { -1, -1, -1 };
      double Cloud_BulkVel[3]                = {  0,  0,  0 };
      double Cloud_Rho0                      = -1;
      double Cloud_R0                        = -1;
      double Cloud_Einasto_Power_Factor      = -1;
      char   DensProf_Table_Name[MAX_STRING];
      long   Cloud_Par_Num                   = -1;
      double Cloud_MaxR                      = -1;
      int    Cloud_NBin                      = -1;
      int    Cloud_RSeed                     = -1;
      int    AddExtPot_Analytical            =  0;
      int    AddExtPot_Table                 =  0;
      char   ExtPot_Table_Name[MAX_STRING];

//    functions to get the physical quantities
      double getDensity                       ( const double r );
      double getAnalEnclosedMass              ( const double r );
      double getExternalPotential             ( const double r );
      double getRandomSampleVelocity          ( const double r );
      double getIntegratedDistributionFunction( const double E );
      void   getRandomVector_GivenLength      ( const double Length, double RandomVector[3] );

//    input table of density profile
      void    loadInputDensProfTable();
      int     InputTable_DensProf_NBin    = -1;
      double *InputTable_DensProf_Radius  = NULL;
      double *InputTable_DensProf_Density = NULL;

//    input table of external potential
      void    loadInputExtPotTable();
      int     InputTable_ExtPot_NBin      = -1;
      double *InputTable_ExtPot_Radius    = NULL;
      double *InputTable_ExtPot_Potential = NULL;

//    arrays of radial distribution of properties of the cloud
      void    constructRadialArray();
      int     RNPoints         = -1;
      int     RLastIdx         = -1;
      double  RArray_dR        = -1;
      double  RArray_MaxR      = -1;
      double *RArray_R         = NULL;
      double *RArray_Rho       = NULL;
      double *RArray_M_Enc     = NULL;
      double *RArray_dRho_dPsi = NULL;
      double *RArray_Phi       = NULL;

//    arrays in energy space
      void    constructEnergyArray();
      int     ENPoints        = -1;
      int     ELastIdx        = -1;
      double  EArray_dE       = -1;
      double  EArray_MinE     = -1;
      double  EArray_MaxE     = -1;
      double *EArray_E        = NULL;
      double *EArray_DFunc    = NULL;
      double *EArray_IntDFunc = NULL;

//    output the arrays
      void    printArrays();

//    cumulative probability distribution at the given radius
      double *CumulProbaDistr_GivenRadius = NULL;

//    random number generator
      RandomNumber_t *Random_Num_Gen;

//    flags to record the action status of function calls
      bool    hasSetCenterAndBulkVel      = false;
      bool    hasSetModelParameters       = false;
      bool    hasSetEinastoPowerFactor    = false;
      bool    hasSetDensProfTableFilename = false;
      bool    hasSetParticleParameters    = false;
      bool    hasSetExtPotParameters      = false;
      bool    hasLoadedInputDensProfTable = false;
      bool    hasLoadedInputExtPotTable   = false;
      bool    hasConstructedRadialArray   = false;
      bool    hasConstructedEnergyArray   = false;
      bool    hasConstructedDistribution  = false;
      bool    hasConstructedParticles     = false;

}; // class Par_EquilibriumIC



#endif //__PAR_EQUILIBRIUM_IC_H__
