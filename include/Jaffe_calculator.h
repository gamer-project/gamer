#ifndef Jaffe_calculator_H
#define Jaffe_calculator_H
#include "GAMER.h"
#include<iostream>

using namespace std;

/***gsl library***/

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#include<sstream>
#include<fstream>

#define size_Jaffe 100000/***difference***/
#define nbin_Jaffe 1000


class Jaffe_calculator
{
    public:
        //Constructor and Destructor
        Jaffe_calculator();
        virtual ~Jaffe_calculator();
        
        //Principal Functions
        void init(double newton_g,double rho,double r0,double maxr);        
        double set_vel(double r);
        double set_mass(double r);
        
        
    protected:
        
    private:
        double **prob_dens;
        double **int_prob_dens;
        double **psi;
        double delta;

        //Statistics
        double ave(double* a,int start,int fin);
        double var_n(double* a,int start,int fin);
        double cor(double* x,double* y,int start,int fin);
        void mask(double* x,int start,int fin);
        void add_num(double* x,int start,int fin);
        double slope(double* a,double* b,int start,int fin);
        void smooth_all(double* x,int start,int fin);

        //Probability Density
        void initialize_mass();
        void initialize_pot();
        void initialize_prob_dens();
        

        RandomNumber_t *RNG ;
};

#endif //Jaffe_calculator_H
