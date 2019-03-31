#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_PatchSharingSameEdge
// Description :  Return the three sibling direction indices sharing the same patch edge as the input one
//
// Note        :  1. Results are stored in SibID[]
//                2. Direction of edge indices (EdgeID) are defined in the same way as the sibling indices
//
// Parameter   :  EdgeID : Sibling index corresponding to the target patch edge
//                SibID  : Output array (with the size [3])
//
// Return      :  SibID[]
//-------------------------------------------------------------------------------------------------------
void TABLE_SiblingSharingSameEdge( const int EdgeID, int SibID[] )
{

#  ifdef GAMER_DEBUG
   if ( EdgeID < 6  ||  EdgeID >= 18 )
      Aux_Error( ERROR_INFO, "EdgeID (%d) not in the accepted range [6 ... 17] !!\n", EdgeID );
#  endif

   switch ( EdgeID )
   {
      case 6 :
         SibID[0] = 0;
         SibID[1] = 2;
         SibID[2] = 6;
         break;

      case 7 :
         SibID[0] = 1;
         SibID[1] = 2;
         SibID[2] = 7;
         break;

      case 8 :
         SibID[0] = 0;
         SibID[1] = 3;
         SibID[2] = 8;
         break;

      case 9 :
         SibID[0] = 1;
         SibID[1] = 3;
         SibID[2] = 9;
         break;

      case 10 :
         SibID[0] = 2;
         SibID[1] = 4;
         SibID[2] = 10;
         break;

      case 11 :
         SibID[0] = 3;
         SibID[1] = 4;
         SibID[2] = 11;
         break;

      case 12 :
         SibID[0] = 2;
         SibID[1] = 5;
         SibID[2] = 12;
         break;

      case 13 :
         SibID[0] = 3;
         SibID[1] = 5;
         SibID[2] = 13;
         break;

      case 14 :
         SibID[0] = 0;
         SibID[1] = 4;
         SibID[2] = 14;
         break;

      case 15 :
         SibID[0] = 0;
         SibID[1] = 5;
         SibID[2] = 15;
         break;

      case 16 :
         SibID[0] = 1;
         SibID[1] = 4;
         SibID[2] = 16;
         break;

      case 17 :
         SibID[0] = 1;
         SibID[1] = 5;
         SibID[2] = 17;
         break;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "EdgeID", EdgeID );
         break;

   } // switch ( EdgeID )

} // FUNCTION : TABLE_PatchSharingSameEdge
