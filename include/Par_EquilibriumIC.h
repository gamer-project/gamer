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
      void setDensityTableFilename( const char* DensityTableFilename );
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
      char   Density_Table_Name[MAX_STRING];
      long   Cloud_Par_Num                   = -1;
      double Cloud_MaxR                      = -1;
      int    Cloud_TableNBin                 = -1;
      int    Cloud_RSeed                     = -1;
      int    AddExtPot                       =  0;
      char   ExtPot_Table_Name[MAX_STRING];
      // Derive physical attributes for particles
      double Set_Mass( const double r );
      double Set_Density( const double r );
      double Set_Velocity( const double r );

      // Initialize physical parameter tables
      void Init_Mass();
      void Init_Pot();
      void Init_Prob_Dens();

      //  Initialization through Table
      void Init_Mass_Table();
      void Init_Pot_Table();

      //  Add External Potential
      void Add_Ext_Pot();

      // Auxiliary functions
      void RanVec_FixRadius( const double r, double RanVec[] );

      // Solve Eddington's equation
      double potential( const double r );
      double Integration_Eng_base( const double Eng, const int N_points );

      double  delta;
      double  Eng_min;
      double *prob_dens;
      double *int_prob_dens;
      double *psi;

      // statistics
      void   SmoothArray( double* array_x, int index_start, int index_end );

      double ArrayCovariance( const double* array_x, const double* array_y,
                              const int index_start, const int n_elements );
      double Slope_LinearRegression( const double* array_x, const double* array_y,
                                     const int index_start, const int n_elements );

      // Tables of particles' attributes
      double *Table_input_radius;
      double *Table_input_density;
      double *Table_input_enclosedmass;
      int     Table_input_NBin = -1;
      double  Table_dr = 0.0;
      double *Table_Radius;
      double *Table_EnclosedMass;
      double *Table_Density;
      double *Table_dRho_dr;
      double *Table_dRho_dx;
      double *Table_Gravity_Field;
      double *Table_Gravity_Potential;

      // Random number generator
      RandomNumber_t *Random_Num_Gen;

}; // class Par_EquilibriumIC



#endif //__PAR_EQUILIBRIUM_IC_H__
