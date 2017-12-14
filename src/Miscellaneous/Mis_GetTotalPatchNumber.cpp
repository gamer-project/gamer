#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTotalPatchNumber
// Description :  Get the total number of real patches among all GAMER ranks 
//
// Note        :  None
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Mis_GetTotalPatchNumber( const int lv )
{

   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   MPI_Allreduce( &amr->NPatchComma[lv][1], &NPatchTotal[lv], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

} // FUNCTION : Mis_GetTotalPatchNumber
