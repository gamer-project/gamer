#ifndef Burkert_calculator_H
#define Burkert_calculator_H
#include "GAMER.h"
<<<<<<< HEAD
#include<iostream>

using namespace std;

/***gsl library***/

#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_vegas.h>
=======
#include<stdio.h>
#include<iostream>
#include<time.h>
using namespace std;



/***gsl library***/

#include <gsl/gsl_integration.h>
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
<<<<<<< HEAD

#include<sstream>
#include<fstream>

#define size_Burkert 1000/***difference***/


=======
#include <random>

extern double Burkert_Rho0;
extern double Burkert_R0;

#define size 1000
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c

class Burkert_calculator
{
    public:
<<<<<<< HEAD
        //Constructor and Destructor
        Burkert_calculator();
        virtual ~Burkert_calculator();
        
        //Principal Functions
        void init(double newton_g,double rho,double r);        
=======
        Burkert_calculator();
        virtual ~Burkert_calculator();
        void init();        
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
        double set_vel(double r);
        
        
    protected:
        
    private:
<<<<<<< HEAD
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
=======
        double f[size];
        double psi[size];
        double delta;


>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
};

#endif //Burkert_calculator_H
