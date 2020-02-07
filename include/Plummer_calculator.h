#ifndef Plummer_calculator_H
#define Plummer_calculator_H
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


#define size_Plummer 1000/***difference***/


class Plummer_calculator
{
    public:
        Plummer_calculator();
        virtual ~Plummer_calculator();
        void init(double newton_g,double rho,double r);      
        double set_vel(double r);
        RandomNumber_t *RNG ;
        
        
    protected:
        
    private:
        double NEWTON_G;
        double Plummer_Rho0;
        double Plummer_R0;


};

#endif //Plummer_calculator_H
