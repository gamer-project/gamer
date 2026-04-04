#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_End
// Description :  Free the resources used by Grackle
//
// Note        :  1. Invoked by End_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_End()
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   delete [] Che_FieldData->grid_dimension;
   delete [] Che_FieldData->grid_start;
   delete [] Che_FieldData->grid_end;

   delete Che_FieldData;   Che_FieldData = NULL;   // Che_FieldData is allocated in Grackle_Init_FieldData.cpp
   delete grackle_data;    grackle_data  = NULL;   // grackle_data is allocated in Grackle_Init.cpp


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Grackle_End



#endif // #ifdef SUPPORT_GRACKLE
