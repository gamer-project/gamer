#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_05
// Description :  Return the number of layers of the NeiFlag_PosList[] and BounFlag_PosList[]
//
// Note        :  Work for Flag_Buffer(), Buf_RecordExchangeDataPatchID() and Buf_RecordBoundaryFlag()
//
// Parameter   :  SibID : Sibling index (0~25)
//-------------------------------------------------------------------------------------------------------
int TABLE_05( const int SibID )
{

   switch ( SibID )
   {
      case 0: case 1: case 2: case 3: case 4: case 5:
         return 2;

      case 6: case 7: case 8: case 9: case 10: case 11: case 12: case 13: case 14: case 15: case 16: case 17:
         return 4;

      case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25:
         return 8;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "SibID", SibID );
         exit(1);

   } // switch ( SibID )

} // FUNCTION : TABLE_05
