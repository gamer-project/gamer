#ifndef Burkert_calculator_H
#define Burkert_calculator_H
#include "GAMER.h"
#include<iostream>

using namespace std;

/***gsl library***/

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>

#include<sstream>
#include<fstream>

#define size_Burkert 1000/***difference***/



class Burkert_calculator
{
    public:
        //Constructor and Destructor
        Burkert_calculator();
        virtual ~Burkert_calculator();
        
        //Principal Functions
        void init(double newton_g,double rho,double r);        
        double set_vel(double r);
        
        
    protected:
        
    private:
        double prob_dens[size_Burkert];
        double int_prob_dens[size_Burkert];
        double psi[size_Burkert];
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
        void initialize_prob_dens();

        RandomNumber_t *RNG ;
};

#endif //Burkert_calculator_H
