#include "Copyright.h"
#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-----------------------------------------------------------------------------------------
// Function    :  CPU_GrackleSolver_Original
// Description :  Update the internal energy by the various cooling and heating mechanisms
//                in the original Grackle library
//
// Note        :  1. Work for "GRACKLE_MODE == GRACKLE_MODE_ORI"
//                   --> No matter GPU is enabled or not
//
// Parameter   :  Che_FieldData : Array of Grackle "grackle_field_data" objects
//                Che_Units     : Grackle "code_units" object
//                NPatchGroup   : Number of patch groups to be evaluated
//                dt            : Time interval to advance solution
//-----------------------------------------------------------------------------------------
void CPU_GrackleSolver_Original( grackle_field_data Che_FieldData[], code_units Che_Units, const int NPatchGroup, const real dt )
{

   const int NPatch = NPatchGroup*8;


// loop over all patches
#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<NPatch; p++)
   {
      if (  solve_chemistry( &Che_Units, &Che_FieldData[p], dt ) == 0  )
         Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );
   }

} // FUNCTION : CPU_GrackleSolver_Original



#endif // #ifdef SUPPORT_GRACKLE
