#include <cstring> 
#include "GAMER.h"
#include <tuple>
#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <sstream>
using namespace std;


double Halo_Profile_Param_a;
double Halo_Profile_Param_b;
double Halo_Profile_Param_c;

// bins in the density profile table
const int nbin = 2500;

// NFW density profile
double NFW_dens(double x)
{
   return 1.0 / ((1.0 + x) * (1.0 + x) * x);
}

// Burkert density profile
double Burkert_dens(double x)
{
   return 1.0 / ((1 + x) * (1 + x * x));
}

// Plummer density profile
double Plummer_dens(double x)
{
   return pow(1 + x * x, -2.5);
}

// Double Power Law density profile
double DoublePowerLaw_dens(double x, double alpha, double beta, double gamma)
{
   return pow(x, -gamma) * pow(1 + pow(x, alpha), (gamma - beta) / alpha);
}

// King density profile
double King_dens(double x)
{
   return pow(1 + x * x, -0.5);
}


// return the density base on the inserted model name
double density(double rho_0, double rs, double r, const char* model_name)
{
   double x = r/rs;
   string modelName = string(model_name);
   if ( modelName == "Burkert" )
   {
      return rho_0 * Burkert_dens(x);
   } else if ( modelName == "NFW" )
   {
      return rho_0 * NFW_dens(x);
   } else if ( modelName == "Plummer" )
   {
      return rho_0 * Plummer_dens(x);
   } else if ( modelName == "King" )
   {
      return rho_0 * King_dens(x);
   } else if ( modelName == "DoublePowerLaw" )
   {
      return rho_0 * DoublePowerLaw_dens(x, Halo_Profile_Param_a,Halo_Profile_Param_b,Halo_Profile_Param_c);
   } else
   {
      return 0.0;
   }
}




// Numerical integration using the trapezoidal rule
double integrate(const function<double(double)>& f, double a, double b, int n = 1000)
{
   double h = (b - a) / n; // Step size
   double sum = 0.5 * (f(a) + f(b)); // Start with end points
   for (int i = 1; i < n; i++)
   {
      sum += f(a + i * h);
   }
   return sum * h;
}


// smooth_transition function using sigmoid
double smooth_transition(double r_value, double start, double end, double scale)
{
   if ( r_value < start )
   {
      return 1.0;
   } else if ( r_value > end )
   {
      return exp(-(r_value - start) / scale);
   } else {
      double x = (r_value - start) / (end - start);
      double sigmoid = 1 / (1 + exp(-10 * (x - 0.5)));
      return 1 - sigmoid * (1 - exp(-(r_value - start) / scale));
   }
}

// Function to calculate cluster mass
double clustermass_soft(double rho_0, double rs, double r, const char* model_name, double begin_smooth_r, double end_smooth_r, double GC_ri)
{
   double x = r / rs;
   auto massbase = [rs, &model_name, begin_smooth_r, end_smooth_r, GC_ri](double x) -> double {
      double smooth_transition_factor = smooth_transition(x * rs, begin_smooth_r, end_smooth_r, GC_ri);
      const double prefactor = 4 * M_PI * x * x * pow(rs, 3);
      string modelName = string(model_name);
      if ( modelName == "Burkert" )
      {
         return prefactor * smooth_transition_factor * Burkert_dens(x);
      } else if ( modelName == "NFW" )
      {
         return prefactor * smooth_transition_factor * NFW_dens(x);
      } else if ( modelName == "Plummer" )
      {
         return prefactor * smooth_transition_factor * Plummer_dens(x);
      } else if ( modelName == "King" )
      {
         return prefactor * smooth_transition_factor * King_dens(x);
      } else if ( modelName == "DoublePowerLaw" ) {
         return prefactor * smooth_transition_factor * DoublePowerLaw_dens(x, Halo_Profile_Param_a, Halo_Profile_Param_b, Halo_Profile_Param_c);
      } else
      {
         return 0.0;
      }
   };
 
   double integrated_mass = integrate(massbase, 0, x) * rho_0; 
   return integrated_mass;
}    

// Function to generate even sample in log scale
vector<double> logspace(double start, double end, int n)
{
   vector<double> result(n);
   double delta = (end - start) / (n - 1);
   for(int i = 0; i < n; ++i)
   {
      result[i] = pow(10, start + i * delta);
   }
   return result;
}



// Function to calculate the initial condition
tuple<vector<double>, vector<double>, vector<double>> Calculate_IC( const char* HaloType, const double GC_MASS, const double GC_R,
								    const double Halo_Rho0, const double Halo_Rs, const double Halo_Rt,
								    const char* Table_Name, const bool Pure_Table){

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// Generate the radius bin 
   vector<double> r = logspace(log10(Halo_Rs / 100), log10(Halo_Rt), nbin);
   double dr = r[1] - r[0];

// Calculate the density and mass
   vector<double> dens(nbin);
   vector<double> mass(nbin);
   for (int i = 0; i < nbin; ++i)
   {
// The smooth transition starts from the GC's initial place and ends at the truncation radius
      dens[i] = density(Halo_Rho0, Halo_Rs, r[i], HaloType) * smooth_transition(r[i], GC_R, Halo_Rt, GC_R);
      mass[i] = clustermass_soft(Halo_Rho0, Halo_Rs, r[i], HaloType, GC_R, Halo_Rt, GC_R);
   }

// Write the result to a file
   if ( MPI_Rank == 0)
   {
      char Filename[MAX_STRING];
      sprintf( Filename, "%s", "Profile_Table" );
      FILE *File = fopen( Filename, "a" );
      if (Time[0]==0.0)
      {
         fprintf(File, "%15s\t%15s\t%15s\n", "Radius", "Density", "Enclosed Mass");
      }
      for (int i=0; i<nbin; ++i)
      {
         fprintf(File, "%15.7e\t%15.7e\t%15.7e\n", r[i],dens[i],mass[i]);
      }
      fclose ( File );	
   }


// Start to calculate the initial properties of the GC
// 1. Find the index of radius
   int index = -1;
   double min_diff = abs(r[0] - GC_R);
   for (int i = 1; i < nbin + 1; ++i)
   {
      double diff = abs(r[i] - GC_R);
      if (diff < min_diff)
      {
         min_diff = diff;
         index = i;
      }
   }

// 2. To be more accurate, we interpolate the radius to find the best enclosed mass value
   double m = (mass[index] - mass[index - 1]) / (r[index] - r[index - 1]);
   double b = mass[index - 1] - m * r[index - 1];
   double interpolated_mass = m * GC_R + b;

// 3. calculate the speed of GC and Halo (considering realtive speed)

   double v_initial = sqrt(NEWTON_G * interpolated_mass / GC_R);

   vector<double> GC_position = { amr->BoxSize[0]*0.5 + GC_R, amr->BoxSize[0]*0.5, amr->BoxSize[0]*0.5};
   vector<double> GC_velocity = { 0.0, v_initial * interpolated_mass / ( GC_MASS + interpolated_mass)   , 0.0};
   vector<double> Halo_velocity = {0.0, -v_initial * GC_MASS / ( GC_MASS + interpolated_mass) , 0.0};

// 4. output the Halo velocity to the Input__Profile_Params -> In order the ParticleEquilibriumIC test problem can read the parameter
   ifstream inputFile("Input__Profile_Params");
   vector<string> lines;
   string line;

   while (getline(inputFile, line))
   {
      lines.push_back(line);
   }
   
   inputFile.close();
   
   
   for (auto& currentLine : lines)
   {
      stringstream ss(currentLine);
      string key;
      ss >> key;
      
      if ( key == "Cloud_BulkVelX" )
      {
         currentLine = "Cloud_BulkVelX            " + to_string(Halo_velocity[0]);
      } else if ( key == "Cloud_BulkVelY" )
      {
         currentLine = "Cloud_BulkVelY            " + to_string(Halo_velocity[1]);
      } else if (key == "Cloud_BulkVelZ")
      {
         currentLine = "Cloud_BulkVelZ            " + to_string(Halo_velocity[2]);
      }
   }
   
   ofstream outputFile("Input__Profile_Params", ios::out | ios::trunc); 
   
   for (const auto & modifiedLine : lines)
   {
      outputFile << modifiedLine << endl;
   }
   
   outputFile.close();




   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... Done\n", __FUNCTION__ );
   return make_tuple(GC_position, GC_velocity, Halo_velocity);
}
