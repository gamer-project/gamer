/***********************************************************************
/
/ Example executable using libgrackle
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
/
/ Modified by: Yuri Oku
************************************************************************/

#include "GAMER.h"

#ifdef SUPPORT_GRACKLE
// specific global variables declared in Init_TestProb_Grackle_Comoving.cpp
extern double              GrackleComoving_InitialMetallicity;
extern grackle_field_data  my_fields;
extern gr_float           *my_temperature;
extern gr_float           *my_gamma;
extern gr_float           *my_cooling_time;



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CoolingCurve
// Description :  Create a cooling table at the specified redshift on the proper coordinate
//
// Note        :  Assuming the cosmic mean density
//
// Parameter   :  z_value : redshift
//                T_min   : minimum temperature in log10(K)
//                T_max   : maximum temperature in log10(K)
//                dT      : temperature interval in log10(K)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_CoolingCurve( const double z_value, const double T_min, const double T_max, const double dT )
{
   double a_value = 1.0 / (1.0 + z_value);

   double rho_crit0 = 3.0 * pow(100.0 * HUBBLE0 * 1e5 / Const_Mpc, 2) / (8.0 * M_PI * Const_NewtonG);
   double rho_mean0 = OMEGA_M0 * rho_crit0;


// no verbose output
   grackle_verbose = 0;


// First, set up the units system.
// These are conversions from code units to cgs.
   code_units my_units;
   my_units.comoving_coordinates = 0;                           // proper coordinate
   my_units.density_units        = rho_mean0 * pow(a_value, -3);
   my_units.length_units         = Const_Mpc;
   my_units.time_units           = Const_Mpc / 1e5;
   my_units.a_units              = 1;                           // units for the expansion factor
   my_units.a_value              = a_value;
   set_velocity_units(&my_units);

   grackle_data->with_radiative_cooling = 0; // turn off cooling for this function


// initialize the chemistry object.
   if ( initialize_chemistry_data(&my_units) == 0 )
      Aux_Error( ERROR_INFO, "Error in initialize_chemistry_data.\n");


// set temperature units
   double temperature_units   = get_temperature_units(&my_units);

   my_fields.density      [0] = 1.0;
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   my_fields.HI_density   [0] = 0.0;
   my_fields.HII_density  [0] = grackle_data->HydrogenFractionByMass * (1.0 - GrackleComoving_InitialMetallicity) * my_fields.density[0];
   my_fields.HeI_density  [0] = 0.0;
   my_fields.HeII_density [0] = 0.0;
   my_fields.HeIII_density[0] = (1.0 - grackle_data->HydrogenFractionByMass) * (1.0 - GrackleComoving_InitialMetallicity) * my_fields.density[0];
   my_fields.e_density    [0] = (my_fields.HII_density[0] + my_fields.HeII_density[0] / 4.0 + 2.0 * my_fields.HeIII_density[0] / 4.0) * Const_me / Const_mp;
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   my_fields.HM_density   [0] = 0.0;
   my_fields.H2I_density  [0] = 0.0;
   my_fields.H2II_density [0] = 0.0;
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   my_fields.DI_density   [0] = 0.0;
   my_fields.DII_density  [0] = grackle_data->DeuteriumToHydrogenRatio * my_fields.HII_density[0];
   my_fields.HDI_density  [0] = 0.0;
   }
   if ( GRACKLE_METAL )
   my_fields.metal_density[0] = GrackleComoving_InitialMetallicity * my_fields.density[0];


   /*********************************************************************
   / Calling the chemistry solver
   / These routines can now be called during the simulation.
   *********************************************************************/

   char ofname[30];
   sprintf(ofname, "coolingrate_z%1.1f.dat", z_value);
   FILE *fp;
   fp = fopen(ofname, "w");

   fprintf(fp, "#%30s%30s\n", "Temperature [K]", "CoolingRate [erg cm^3 s^-1]");


// copy buffer for checking convergence
   double buf_density;
   double buf_HI_density, buf_HII_density, buf_HeI_density, buf_HeII_density, buf_HeIII_density, buf_e_density;
   double buf_HM_density, buf_H2I_density, buf_H2II_density;
   double buf_DI_density, buf_DII_density, buf_HDI_density;

   int maxiter = 0;
   for (double logT=T_min; logT<=T_max; logT+=dT)
   {
//    initial guess for the internal energy
      double mu_init    = 1.0;
      double gamma_init = 5.0/3.0;
      my_fields.internal_energy[0] = pow(10.0, logT) / (mu_init * (gamma_init - 1.) * temperature_units);

      int iter = 0;
      while( true )
      {
//       Copy the fields to check for convergence.
         buf_density       = my_fields.density      [0];
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
         buf_HI_density    = my_fields.HI_density   [0];
         buf_HII_density   = my_fields.HII_density  [0];
         buf_HeI_density   = my_fields.HeI_density  [0];
         buf_HeII_density  = my_fields.HeII_density [0];
         buf_HeIII_density = my_fields.HeIII_density[0];
         buf_e_density     = my_fields.e_density    [0];
         }
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
         buf_HM_density    = my_fields.HM_density   [0];
         buf_H2I_density   = my_fields.H2I_density  [0];
         buf_H2II_density  = my_fields.H2II_density [0];
         }
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
         buf_DI_density    = my_fields.DI_density   [0];
         buf_DII_density   = my_fields.DII_density  [0];
         buf_HDI_density   = my_fields.HDI_density  [0];
         }


//       Calculate cooling time.
         if ( calculate_cooling_time(&my_units, &my_fields, my_cooling_time) == 0 )
            Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n");

         double dt = fabs(0.1 * my_cooling_time[0]);


//       Evolve the chemistry. cooling is disabled by with_radiative_cooling = 0
         if ( solve_chemistry(&my_units, &my_fields, dt) == 0 )
            Aux_Error( ERROR_INFO, "Error in solve_chemistry.\n");


//       reset the internal energy
         if ( calculate_temperature(&my_units, &my_fields, my_temperature) == 0 )
            Aux_Error( ERROR_INFO, "Error in calculate_temperature.\n");

         if ( calculate_gamma(&my_units, &my_fields, my_gamma) == 0 )
            Aux_Error( ERROR_INFO, "Error in calculate_gamma.\n");


         double mu                    = my_temperature[0] / (my_fields.internal_energy[0] * (my_gamma[0] - 1.) * temperature_units);
         my_fields.internal_energy[0] = pow(10.0, logT) / (mu * (my_gamma[0] - 1.) * temperature_units);

         if ( calculate_temperature(&my_units, &my_fields, my_temperature) == 0 )
            Aux_Error( ERROR_INFO, "Error in calculate_temperature.\n");


//       check convergence
         double diff = 0.0;
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
            diff = fmax(diff, fabs(my_fields.HI_density   [0] - buf_HI_density   ) / my_fields.HI_density   [0]);
            diff = fmax(diff, fabs(my_fields.HII_density  [0] - buf_HII_density  ) / my_fields.HII_density  [0]);
            diff = fmax(diff, fabs(my_fields.HeI_density  [0] - buf_HeI_density  ) / my_fields.HeI_density  [0]);
            diff = fmax(diff, fabs(my_fields.HeII_density [0] - buf_HeII_density ) / my_fields.HeII_density [0]);
            diff = fmax(diff, fabs(my_fields.HeIII_density[0] - buf_HeIII_density) / my_fields.HeIII_density[0]);
            diff = fmax(diff, fabs(my_fields.e_density    [0] - buf_e_density    ) / my_fields.e_density    [0]);
         }
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
            diff = fmax(diff, fabs(my_fields.HM_density   [0] - buf_HM_density   ) / my_fields.HM_density   [0]);
            diff = fmax(diff, fabs(my_fields.H2I_density  [0] - buf_H2I_density  ) / my_fields.H2I_density  [0]);
            diff = fmax(diff, fabs(my_fields.H2II_density [0] - buf_H2II_density ) / my_fields.H2II_density [0]);
         }
         if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
            diff = fmax(diff, fabs(my_fields.DI_density   [0] - buf_DI_density   ) / my_fields.DI_density   [0]);
            diff = fmax(diff, fabs(my_fields.DII_density  [0] - buf_DII_density  ) / my_fields.DII_density  [0]);
            diff = fmax(diff, fabs(my_fields.HDI_density  [0] - buf_HDI_density  ) / my_fields.HDI_density  [0]);
         }
         diff = fmax(diff, fabs(my_temperature         [0] - pow(10.0, logT)  ) / my_temperature         [0]);

         iter++;

         if ( diff < 1.e-3 )   break; // converged
      }
      maxiter = iter > maxiter ? iter : maxiter;


//    Calculate cooling time.
      if ( calculate_cooling_time(&my_units, &my_fields, my_cooling_time) == 0 )
         Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n");


      double mu    = my_temperature[0] / (my_fields.internal_energy[0] * (my_gamma[0] - 1.) * temperature_units);
      double nden  = my_fields.density[0] / mu / Const_mp * my_units.density_units;
      double eden  = my_fields.density[0] * my_fields.internal_energy[0] * my_units.density_units * my_units.velocity_units * my_units.velocity_units;
      double lcool = eden / fabs(my_cooling_time[0] * my_units.time_units) / nden / nden;

      fprintf(fp, " %30.8e%30.8e\n", my_temperature[0], lcool);
   }

   fclose(fp);

   return;
}
#endif // SUPPORT_GRACKLE
