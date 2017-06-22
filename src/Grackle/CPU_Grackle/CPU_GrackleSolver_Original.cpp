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
// Parameter   :  Che_FieldData :
//                Che_Units     :
//                NPatchGroup   : Number of patch groups to be evaluated
//                dt            : Time interval to advance solution
//-----------------------------------------------------------------------------------------
void CPU_GrackleSolver_Original( grackle_field_data Che_FieldData[], const code_units Che_Units,
                                 const int NPatchGroup, const real dt )
{

   const int NPatch = NPatchGroup*8;


// loop over all patches
#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatch; P++)
   {

   } // for (int P=0; P<NPatch; P++)

} // FUNCTION : CPU_GrackleSolver_Original



#endif // #ifdef SUPPORT_GRACKLE
