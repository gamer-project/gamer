#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-----------------------------------------------------------------------------------------
// Function    :  CPU_GrackleSolver
// Description :  Update the internal energy by the various cooling and heating mechanisms
//                in the original Grackle library
//
// Note        :  1. Currently it is used even when GPU is enabled
//
// Parameter   :  Che_FieldData : Array of Grackle "grackle_field_data" objects
//                Che_Units     : Grackle "code_units" object
//                NPatchGroup   : Number of patch groups to be evaluated
//                dt            : Time interval to advance solution
//-----------------------------------------------------------------------------------------
void CPU_GrackleSolver( grackle_field_data *Che_FieldData, code_units Che_Units, const int NPatchGroup, const real dt )
{

// set grid_dimension, grid_start, and grid_end
   const int OptFac = 16;  // optimization factor
   if ( SQR(PS2)%OptFac != 0 )   Aux_Error( ERROR_INFO, "SQR(PS2) %% OptFac != 0 !!\n" );

   Che_FieldData->grid_dimension[0] = PS2*OptFac;
   Che_FieldData->grid_dimension[1] = 1;
   Che_FieldData->grid_dimension[2] = SQR(PS2)*NPatchGroup/OptFac;

   for (int d=0; d<3; d++)
   {
      Che_FieldData->grid_start[d] = 0;
      Che_FieldData->grid_end  [d] = Che_FieldData->grid_dimension[d] - 1;
   }

   // ====================================================================================================
   // Calculate cooling time.
   if ( MPI_Rank == 0 ){
      Aux_Message( stdout, "=======================================================\n");
      Aux_Message( stdout, "MPI_Rank_0_NPatchGroup = %d\n" , NPatchGroup);
   }

   // Calculate temperature.
   gr_float *temperature;
   temperature = new gr_float[NPatchGroup*SQR(PS2)*PS2];
   if(calculate_temperature(&Che_Units, Che_FieldData,temperature) == 0)
      Aux_Error( ERROR_INFO, "Error in calculate_temperature.\n" );
   if ( MPI_Rank == 0 ){
      Aux_Message( stdout, "temperature = %f\n" , temperature[0]);
      Aux_Message( stdout, "********************************************************\n");
   }


   double dust_dens_max = 1e-300, dust_dens_min = 1e300;
   double metal_dens_max = 1e-300, metal_dens_min = 1e300;
   for (size_t i = 0; i < (size_t)NPatchGroup*SQR(PS2)*PS2; ++i){
         if(Che_FieldData->dust_density[i] > dust_dens_max){
            dust_dens_max = Che_FieldData->dust_density[i];
         }
         if(Che_FieldData->dust_density[i] < dust_dens_min){
            dust_dens_min = Che_FieldData->dust_density[i];
         }
         if(Che_FieldData->metal_density[i] > metal_dens_max){
            metal_dens_max = Che_FieldData->metal_density[i];
         }
         if(Che_FieldData->metal_density[i] < metal_dens_min){
            metal_dens_min = Che_FieldData->metal_density[i];
         }
   }
   if(MPI_Rank == 0){
      Aux_Message( stdout, "dust_dens_max  = %f\n" , dust_dens_max);
      Aux_Message( stdout, "dust_dens_min  = %f\n" , dust_dens_min);
      Aux_Message( stdout, "metal_dens_max = %f\n" , metal_dens_max);
      Aux_Message( stdout, "metal_dens_min = %f\n" , metal_dens_min);
   }


   // double dens_max = 1e-300, dens_min = 1e300;
   // double inter_max = 1e-300, inter_min = 1e300;
   // for (size_t i = 0; i < (size_t)NPatchGroup*SQR(PS2)*PS2; ++i){
   //    if(Che_FieldData->density[i] > dens_max){
   //       dens_max = Che_FieldData->density[i];
   //    }
   //    if(Che_FieldData->density[i] < dens_min){
   //       dens_min = Che_FieldData->density[i];
   //    }
   //    if(Che_FieldData->internal_energy[i] > inter_max){
   //       inter_max = Che_FieldData->internal_energy[i];
   //    }
   //    if(Che_FieldData->internal_energy[i] < inter_min){
   //       inter_min = Che_FieldData->internal_energy[i];
   //    }
   // }
   // if(MPI_Rank == 0){
   //    Aux_Message( stdout, "dens_max = %f\n" , dens_max);
   //    Aux_Message( stdout, "dens_min = %f\n" , dens_min);
   //    Aux_Message( stdout, "internal_energy_max = %f\n" , inter_max);
   //    Aux_Message( stdout, "internal_energy_min = %f\n" , inter_min);
   // }
   // ====================================================================================================

// invoke Grackle
// --> note that we use the OpenMP implementation in Grackle directly, which applies the parallelization to the first two
//     dimensiones of the input grid
// --> this approach is found to be much more efficient than parallelizing different patches or patch groups here

   if (  solve_chemistry( &Che_Units, Che_FieldData, dt ) == 0  )
      Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );

   // gr_float * cooling_time = new gr_float[NPatchGroup*SQR(PS2)*PS2];
   // double dt_remain = dt;
   // double dt_sub;
   // while(dt_remain > 0.0){
   //    if(calculate_cooling_time(&Che_Units, Che_FieldData, cooling_time) == 0)
   //       Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n" );
   //    // double tmin = 1e300;
   //    // for (size_t i = 0; i < (size_t)NPatchGroup*SQR(PS2)*PS2; ++i){
   //    //    if (cooling_time[i] > 0.0 && cooling_time[i] < tmin) tmin = cooling_time[i];
   //    // }
   //    double dt_cool = 0.1 * fabs(cooling_time[0]);
   //    dt_sub = (dt_cool < dt_remain) ? dt_cool : dt_remain;
      
   //    if ( MPI_Rank == 0 ){
   //       // Aux_Message( stdout, "cooling_time   = %.15E \n", dt_cool);
   //       Aux_Message( stdout, "dt_sub_cooling = %.15E \n", dt_sub);
   //    }
   //    if (  solve_chemistry( &Che_Units, Che_FieldData, dt_sub ) == 0  )
   //       Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );
      
   //    dt_remain -= dt_sub;
   // }
   // delete [] cooling_time;
   // delete [] temperature;
} // FUNCTION : CPU_GrackleSolver


#endif // #ifdef SUPPORT_GRACKLE
