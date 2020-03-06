#ifndef UNKNOWN_calculator_H
#define UNKNOWN_calculator_H
#include "GAMER.h"
#include<iostream>
#include<fstream>
using namespace std;

/***gsl library***/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define size_UNKNOWN 1000



class UNKNOWN_calculator
{
    public:
        UNKNOWN_calculator();
        virtual ~UNKNOWN_calculator();
        void init(double newton_g,int r_col,int rho_col,const char* Filename);        
        double set_vel(double r);
        RandomNumber_t *RNG ;

        void initialize_mass(int row);
        void initialize_pot(int row);
        void initialize_prob_dens();
        
        //statistics
        double slope(double* a,double* b,int start,int fin);
        void smooth_all(double* x,int start,int fin);
        double set_mass( double x );
        
    protected:
        
    private:
        double MassProf_UNKNOWN( const double r );
        double integration_eng_base_UNKNOWN(double eng);

        double prob_dens[size_UNKNOWN];
        double int_prob_dens[size_UNKNOWN];
        double psi[size_UNKNOWN];
        double delta;
        double eng_min_UNKNOWN;

        //statistics
        double ave(double* a,int start,int fin);
        double var_n(double* a,int start,int fin);
        double cor(double* x,double* y,int start,int fin);
        void mask(double* x,int start,int fin);
        void add_num(double* x,int start,int fin);


};

#endif //UNKNOWN_calculator_H
