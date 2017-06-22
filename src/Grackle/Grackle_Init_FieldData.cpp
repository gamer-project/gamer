#include "Copyright.h"
#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Init_FieldData
// Description :  Initialize the "Che_FieldData" grackle_field_data objects of Grackle
//
// Note        :  1. Must be called AFTER Init_Load_Parameter()
//                2. Invoked by Init_MemAllocate()
//                   --> "Che_FieldData" is freed by End_MemFree()
//                3. Useless if GRACKLE_MODE == GRACKLE_MODE_GAMER
//
// Parameter   :  Che_NPG : Number of patch groups to be evaluated at a time
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Grackle_Init_FieldData( const int Che_NPG )
{

// nothing to do if we are not using the original Grackle
   if ( GRACKLE_MODE != GRACKLE_MODE_ORI )   return;


// allocate memory
   const int Che_NP = 8*Che_NPG;

   Che_FieldData = new grackle_field_data [Che_NP];


// initialization
   const int NDim = 3;

   for (int t=0; t<Che_NP; t++)
   {
      Che_FieldData[t].grid_rank      = NDim;
      Che_FieldData[t].grid_dimension = new int [NDim];
      Che_FieldData[t].grid_start     = new int [NDim];
      Che_FieldData[t].grid_end       = new int [NDim];
      Che_FieldData[t].grid_dx        = NULL_REAL;           // set by Grackle_Prepare()

      for (int d=0; d<NDim; d++)
      {
        Che_FieldData[t].grid_dimension[d] = PS1;
        Che_FieldData[t].grid_start    [d] = 0;
        Che_FieldData[t].grid_end      [d] = PS1 - 1;
      }
   }

} // FUNCTION : Grackle_Init_FieldData



#endif // #ifdef SUPPORT_GRACKLE
