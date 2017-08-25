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

// loop over all patch groups
// --> note that we use the OpenMP implementation in Grackle directly, which applies the parallelization to different cells
//     within the same patch group
// --> this approach is found to be much more efficient than parallelizing different patches or patch groups here
//#  pragma omp parallel for schedule( runtime )
   for (int p=0; p<NPatchGroup; p++)
   {
      if (  solve_chemistry( &Che_Units, &Che_FieldData[p], dt ) == 0  )
         Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );
   }

} // FUNCTION : CPU_GrackleSolver_Original



#endif // #ifdef SUPPORT_GRACKLE
