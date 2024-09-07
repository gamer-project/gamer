/***********************************************************************
/
/ Example executable using libgrackle
/
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#include <grackle.h>

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16
#define Mpc    3.08567758e24
#define G      6.67259e-8

// create a cooling table at the specified redshift
// assuming the cosmic mean density
int main(int argc, char *argv[])
{

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  if (argc < 3) {
    fprintf(stderr, "Usage: %s z_value log10(T_min) log10(T_max) dlog10T\n", argv[0]);
    return EXIT_FAILURE;
  }

  double z_value = atof(argv[1]);
  double T_min = atof(argv[2]);
  double T_max = atof(argv[3]);
  double dT = atof(argv[4]);

  double a_value = 1.0 / (1.0 + z_value);

  double Metallicity =  1.295e-2;

  double H0 = 70.0;
  double omega_m = 0.3;
  double omega_l = 0.7;

  double rho_crit0 = 3.0 * pow(H0 * 1e5 / Mpc, 2) / (8.0 * M_PI * G);
  double rho_mean0 = omega_m * rho_crit0;


  // Enable output
  grackle_verbose = 0;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = rho_mean0 * pow(a_value, -3);
  my_units.length_units = Mpc;
  my_units.time_units = Mpc / 1e5; 
  my_units.a_units = 1; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = a_value;
  set_velocity_units(&my_units);

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  chemistry_data *my_grackle_data;
  my_grackle_data = malloc(sizeof(chemistry_data));
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }

  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  grackle_data->use_grackle = 1;            // chemistry on
  grackle_data->with_radiative_cooling = 0; // cooling on
  grackle_data->primordial_chemistry = 0;   // molecular network with H, He, D
  grackle_data->dust_chemistry = 0;         // dust processes
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 0;           // UV background on
  grackle_data->grackle_data_file = "../CloudyData_noUVB.h5"; // data file
  // grackle_data->grackle_data_file = "../CloudyData_UVB=HM2012.h5"; // data file
  grackle_data->three_body_rate = 4;        // three-body H2 formation
  grackle_data->cie_cooling = 1;            // collisional ionization equilibrium
  grackle_data->max_iterations = 100000;

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  double tiny_number = 1.e-20;

  // Create struct for storing grackle field data
  grackle_field_data my_fields, my_fields_copy;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_start = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_end = malloc(my_fields.grid_rank * sizeof(int));
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  int i;
  for (i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;

  my_fields.density         = malloc(field_size * sizeof(gr_float));
  my_fields.internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields.x_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields.y_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields.z_velocity      = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields.HI_density      = malloc(field_size * sizeof(gr_float));
  my_fields.HII_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HeI_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HeII_density    = malloc(field_size * sizeof(gr_float));
  my_fields.HeIII_density   = malloc(field_size * sizeof(gr_float));
  my_fields.e_density       = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields.HM_density      = malloc(field_size * sizeof(gr_float));
  my_fields.H2I_density     = malloc(field_size * sizeof(gr_float));
  my_fields.H2II_density    = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields.DI_density      = malloc(field_size * sizeof(gr_float));
  my_fields.DII_density     = malloc(field_size * sizeof(gr_float));
  my_fields.HDI_density     = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields.metal_density   = malloc(field_size * sizeof(gr_float));

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = malloc(field_size * sizeof(gr_float));
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = malloc(field_size * sizeof(gr_float));

  // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
  my_fields.RT_HI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeI_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_HeII_ionization_rate = malloc(field_size * sizeof(gr_float));
  my_fields.RT_H2_dissociation_rate = malloc(field_size * sizeof(gr_float));
  // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  my_fields.RT_heating_rate = malloc(field_size * sizeof(gr_float));

  my_fields_copy.density         = malloc(field_size * sizeof(gr_float));
  my_fields_copy.internal_energy = malloc(field_size * sizeof(gr_float));
  my_fields_copy.x_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields_copy.y_velocity      = malloc(field_size * sizeof(gr_float));
  my_fields_copy.z_velocity      = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 1
  my_fields_copy.HI_density      = malloc(field_size * sizeof(gr_float));
  my_fields_copy.HII_density     = malloc(field_size * sizeof(gr_float));
  my_fields_copy.HeI_density     = malloc(field_size * sizeof(gr_float));
  my_fields_copy.HeII_density    = malloc(field_size * sizeof(gr_float));
  my_fields_copy.HeIII_density   = malloc(field_size * sizeof(gr_float));
  my_fields_copy.e_density       = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 2
  my_fields_copy.HM_density      = malloc(field_size * sizeof(gr_float));
  my_fields_copy.H2I_density     = malloc(field_size * sizeof(gr_float));
  my_fields_copy.H2II_density    = malloc(field_size * sizeof(gr_float));
  // for primordial_chemistry >= 3
  my_fields_copy.DI_density      = malloc(field_size * sizeof(gr_float));
  my_fields_copy.DII_density     = malloc(field_size * sizeof(gr_float));
  my_fields_copy.HDI_density     = malloc(field_size * sizeof(gr_float));
  // for metal_cooling = 1
  my_fields_copy.metal_density   = malloc(field_size * sizeof(gr_float));

  gr_float *temperature;
  temperature = malloc(field_size * sizeof(gr_float));
  gr_float *gamma;
  gamma = malloc(field_size * sizeof(gr_float));
  gr_float *cooling_time;
  cooling_time = malloc(field_size * sizeof(gr_float));

  // set temperature units
  double temperature_units = get_temperature_units(&my_units);

  for (i = 0;i < field_size;i++) {
    my_fields.density[i] = 1.0;
    my_fields.HI_density[i] = grackle_data->HydrogenFractionByMass * my_fields.density[i];
    my_fields.HII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HM_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) *
      my_fields.density[i];
    my_fields.HeII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeIII_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2I_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2II_density[i] = tiny_number * my_fields.density[i];
    my_fields.DI_density[i] = grackle_data->HydrogenFractionByMass * grackle_data->DeuteriumToHydrogenRatio * my_fields.density[i];
    my_fields.DII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HDI_density[i] = tiny_number * my_fields.density[i];
    my_fields.e_density[i] = tiny_number * my_fields.density[i];
    // solar metallicity
    my_fields.metal_density[i] = Metallicity * my_fields.density[i];

    my_fields.x_velocity[i] = 0.0;
    my_fields.y_velocity[i] = 0.0;
    my_fields.z_velocity[i] = 0.0;

    my_fields.internal_energy[i] = pow(10.0, T_min) / temperature_units;

    my_fields.volumetric_heating_rate[i] = 0.0;
    my_fields.specific_heating_rate[i] = 0.0;

    my_fields.RT_HI_ionization_rate[i] = 0.0;
    my_fields.RT_HeI_ionization_rate[i] = 0.0;
    my_fields.RT_HeII_ionization_rate[i] = 0.0;
    my_fields.RT_H2_dissociation_rate[i] = 0.0;
    my_fields.RT_heating_rate[i] = 0.0;
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  char ofname[256];
  sprintf(ofname, "coolingrate_z%1.1f.dat", z_value);
  FILE *fp;
  fp = fopen(ofname, "w");

  fprintf(fp, "#%16s %16s\n", "Temperature [K]", "CoolingRate [erg cm^3 s^-1]");

  int maxiter = 0;
  for (double logT = T_min; logT <= T_max; logT += dT)
  {
    for (i = 0;i < field_size;i++) {
      my_fields.internal_energy[i] = pow(10.0, logT) / temperature_units;
    }

    int converged = 0;
    int iter = 0;

    while(converged == 0)
    {
      // Copy the fields to check for convergence.
      for (i = 0;i < field_size;i++) {
        my_fields_copy.density[i] = my_fields.density[i];
        my_fields_copy.HI_density[i] = my_fields.HI_density[i];
        my_fields_copy.HII_density[i] = my_fields.HII_density[i];
        my_fields_copy.HeI_density[i] = my_fields.HeI_density[i];
        my_fields_copy.HeII_density[i] = my_fields.HeII_density[i];
        my_fields_copy.HeIII_density[i] = my_fields.HeIII_density[i];
        my_fields_copy.HM_density[i] = my_fields.HM_density[i];
        my_fields_copy.H2I_density[i] = my_fields.H2I_density[i];
        my_fields_copy.H2II_density[i] = my_fields.H2II_density[i];
        my_fields_copy.DI_density[i] = my_fields.DI_density[i];
        my_fields_copy.DII_density[i] = my_fields.DII_density[i];
        my_fields_copy.HDI_density[i] = my_fields.HDI_density[i];
        my_fields_copy.e_density[i] = my_fields.e_density[i];
      }

      // Calculate cooling time.
      if (calculate_cooling_time(&my_units, &my_fields, cooling_time) == 0) {
        fprintf(stderr, "Error in calculate_cooling_time.\n");
        return EXIT_FAILURE;
      }

      double dt = fabs(0.1 * cooling_time[0]);
      
      // Evolve the chemistry. cooling is disabled by with_radiative_cooling = 0
      if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
        fprintf(stderr, "Error in solve_chemistry.\n");
        return EXIT_FAILURE;
      }

      // check for convergence
      double diff = 0.0;
      for (i = 0;i < field_size;i++) {
        diff = fmax(diff, fabs(my_fields.HI_density[i] - my_fields_copy.HI_density[i]) / my_fields.HI_density[i]);
        diff = fmax(diff, fabs(my_fields.HII_density[i] - my_fields_copy.HII_density[i]) / my_fields.HII_density[i]);
        diff = fmax(diff, fabs(my_fields.HeI_density[i] - my_fields_copy.HeI_density[i]) / my_fields.HeI_density[i]);
        diff = fmax(diff, fabs(my_fields.HeII_density[i] - my_fields_copy.HeII_density[i]) / my_fields.HeII_density[i]);
        diff = fmax(diff, fabs(my_fields.HeIII_density[i] - my_fields_copy.HeIII_density[i]) / my_fields.HeIII_density[i]);
        diff = fmax(diff, fabs(my_fields.HM_density[i] - my_fields_copy.HM_density[i]) / my_fields.HM_density[i]);
        diff = fmax(diff, fabs(my_fields.H2I_density[i] - my_fields_copy.H2I_density[i]) / my_fields.H2I_density[i]);
        diff = fmax(diff, fabs(my_fields.H2II_density[i] - my_fields_copy.H2II_density[i]) / my_fields.H2II_density[i]);
        diff = fmax(diff, fabs(my_fields.DI_density[i] - my_fields_copy.DI_density[i]) / my_fields.DI_density[i]);
        diff = fmax(diff, fabs(my_fields.DII_density[i] - my_fields_copy.DII_density[i]) / my_fields.DII_density[i]);
        diff = fmax(diff, fabs(my_fields.HDI_density[i] - my_fields_copy.HDI_density[i]) / my_fields.HDI_density[i]);
        diff = fmax(diff, fabs(my_fields.e_density[i] - my_fields_copy.e_density[i]) / my_fields.e_density[i]);
      }

      // reset the internal energy
      if (calculate_temperature(&my_units, &my_fields, temperature) == 0) {
        fprintf(stderr, "Error in calculate_temperature.\n");
        return EXIT_FAILURE;
      }
      if (calculate_gamma(&my_units, &my_fields, gamma) == 0) {
        fprintf(stderr, "Error in calculate_gamma.\n");
        return EXIT_FAILURE;
      }
      double mu = temperature[0] / (my_fields.internal_energy[0] * (gamma[0] - 1.) * temperature_units);
      my_fields.internal_energy[0] = pow(10.0, logT) / (mu * (gamma[0] - 1.) * temperature_units);
      my_fields.metal_density[0] = grackle_data->SolarMetalFractionByMass * my_fields.density[0];


      // Calculate temperature.
      if (calculate_temperature(&my_units, &my_fields, temperature) == 0) {
        fprintf(stderr, "Error in calculate_temperature.\n");
        return EXIT_FAILURE;
      }

      diff = fmax(diff, fabs(temperature[0] - pow(10.0, logT)) / temperature[0]);

      if (diff < 1.e-3)
        converged = 1;
      
      // printf("logT = %e, iter = %d, diff = %e\n", logT, iter, diff);
      // printf("HI_density = %e, HI_density_copy = %e\n", my_fields.HI_density[0]/my_fields.density[0], my_fields_copy.HI_density[0]/my_fields_copy.density[0]);

      iter++;
    }
    maxiter = iter > maxiter ? iter : maxiter;

    // Calculate cooling time.
    if (calculate_cooling_time(&my_units, &my_fields, cooling_time) == 0) {
      fprintf(stderr, "Error in calculate_cooling_time.\n");
      return EXIT_FAILURE;
    }

    double mu = temperature[0] / (my_fields.internal_energy[0] * (gamma[0] - 1.) * temperature_units);

    double nden = my_fields.density[0] / mu / mh * my_units.density_units;
    double eden = my_fields.density[0] * my_fields.internal_energy[0] * my_units.density_units * my_units.velocity_units * my_units.velocity_units;
    double lcool = eden / fabs(cooling_time[0] * my_units.time_units) / nden / nden;

    fprintf(fp, "%16.8e %16.8e\n", temperature[0], lcool);
  }

  printf("maxiter = %d\n", maxiter);

  fclose(fp);

  return EXIT_SUCCESS;
}
