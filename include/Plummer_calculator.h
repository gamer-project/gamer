#ifndef Plummer_calculator_H
#define Plummer_calculator_H
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


#define size_Plummer 1000/***difference***/
=======
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_dawson.h>

extern double Plummer_Rho0;
extern double Plummer_GasMFrac;
extern double Plummer_R0;
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c


class Plummer_calculator
{
    public:
        Plummer_calculator();
        virtual ~Plummer_calculator();
<<<<<<< HEAD
        void init(double newton_g,double rho,double r);      
        double set_vel(double r);
        RandomNumber_t *RNG ;
=======
        double pressure(double r);
        double psi(double r);
        double prob(double v,void* r);
        double max_prob(double r);
        //double distribution(double energy_negative);
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c
        
        
    protected:
        
    private:
<<<<<<< HEAD
        double NEWTON_G;
        double Plummer_Rho0;
        double Plummer_R0;
=======
    
>>>>>>> 59094c0a60c1e0583dd0b96ce0541562dc013a7c


};

#endif //Plummer_calculator_H
