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

// invoke Grackle
// --> note that we use the OpenMP implementation in Grackle directly, which applies the parallelization to the first two
//     dimensiones of the input grid
// --> this approach is found to be much more efficient than parallelizing different patches or patch groups here
   if (  solve_chemistry( &Che_Units, Che_FieldData, dt ) == 0  )
      Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );

} // FUNCTION : CPU_GrackleSolver



#endif // #ifdef SUPPORT_GRACKLE
