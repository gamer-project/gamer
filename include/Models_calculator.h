#ifndef Models_calculator_H
#define Models_calculator_H
#include "GAMER.h"
#include<iostream>
#include<fstream>
using namespace std;

/***gsl library***/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define size_Models 1000



class Models_calculator
{
    public:
        Models_calculator();
        virtual ~Models_calculator();
        void init(double newton_g,int r_col,int rho_col,const char* Filename);        
        double set_vel(double r);
        double set_vel_test(double r)
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
        double MassProf_Models( const double r );
        double integration_eng_base_Models(double eng);

        double prob_dens[size_Models];
        double int_prob_dens[size_Models];
        double psi[size_Models];
        double delta;
        double eng_min_Models;

        //statistics
        double ave(double* a,int start,int fin);
        double var_n(double* a,int start,int fin);
        double cor(double* x,double* y,int start,int fin);
        void mask(double* x,int start,int fin);
        void add_num(double* x,int start,int fin);


};

#endif //Models_calculator_H
