#ifndef Hernquist_calculator_H
#define Hernquist_calculator_H
#include "GAMER.h"
#include<iostream>
#include<fstream>
using namespace std;

/***gsl library***/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#define size_Hernquist 1000



class Hernquist_calculator
{
    public:
        Hernquist_calculator();
        virtual ~Hernquist_calculator();
        void init(double newton_g,double rho,double r0,int nbin,double rmax,int rseed,bool trunc_flag=false,double trunc_fac=0.7);        
        double set_vel(double r);
        RandomNumber_t *RNG ;

        void initialize_mass();
        void initialize_pot();
        void initialize_prob_dens();
        
        //statistics
        double slope(double* a,double* b,int start,int fin);
        void smooth_all(double* x,int start,int fin);
        double set_mass( double x );
        
    protected:
        
    private:
        double MassProf_Hernquist( const double r );
        double integration_eng_base_Hernquist(double eng);

        double prob_dens[size_Hernquist];
        double int_prob_dens[size_Hernquist];
        double psi[size_Hernquist];
        double delta;

        //truncation
        bool Trunc_Flag;
        double Trunc_Fac;

        //statistics
        double ave(double* a,int start,int fin);
        double var_n(double* a,int start,int fin);
        double cor(double* x,double* y,int start,int fin);
        void mask(double* x,int start,int fin);
        void add_num(double* x,int start,int fin);


};

#endif //Hernquist_calculator_H
