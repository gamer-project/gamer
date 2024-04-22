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


typedef struct Filename_Parameter{

   int Cloud_Num;
   vector<string> Params_Filenames;

}FP;

typedef struct Physical_Parameter{

   int      Cloud_Num;

   string Params_Filenames;

   char     Cloud_Type[MAX_STRING];
   char     Density_Table_Name[MAX_STRING];
   int      AddExtPot;
   char     ExtPot_Table_Name[MAX_STRING];
   int      Cloud_RSeed;
   double   Cloud_Rho0;
   double   Cloud_R0;
   double   Cloud_MaxR;
   double*  Cloud_Center;
   double*  Cloud_BulkVel;
   double   Cloud_Einasto_Power_Factor;
   long     Cloud_Par_Num;

   int      Cloud_MassProfNBin;

}PhysP;


class Par_EquilibriumIC
{
   public:
      Par_EquilibriumIC();
      virtual ~Par_EquilibriumIC();
      void Read_Filenames( const char *filename_para );
      void Load_Physical_Params( const FP filenames,const int cloud_idx, const long NPar_AllRank );
      void Init();
      void Par_SetEquilibriumIC( real *Mass_AllRank, real *Pos_AllRank[3], real *Vel_AllRank[3], const long Par_Idx );

      PhysP params;
      FP   filenames;

   protected:

   private:
      // Derive physical attributes for particles
      double Set_Mass( const double r );
      double Set_Density( const double r );
      double Set_Velocity(const double x);

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
      int  GetParams( const char *filename, const char *keyword, const int para_num, const char *para_type, vector <string> &container);
      void RanVec_FixRadius( const double r, double RanVec[] );

      // Solve Eddington's equation
      double potential( const double x );
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
      double *Table_r;
      double *Table_Enclosed_Mass;
      double *Table_Density;
      double *Table_dRho_dr;
      double *Table_dRho_dx;
      double *Table_Gravity_Field;
      double *Table_Gravity_Potential;

      // Random number generator
      RandomNumber_t *Random_Num_Gen;

}; // class Par_EquilibriumIC



#endif //__PAR_EQUILIBRIUM_IC_H__
