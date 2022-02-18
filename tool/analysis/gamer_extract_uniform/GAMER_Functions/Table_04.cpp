#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_04
// Description :  Return the number of patches in each direction (different sibling indices) adjacent to
//                a patch group
//
// Note        :  Work for the functions "Prepare_InputGPUArray", "Buf_RecordExchangeDataPatchID",
//                "Buf_RecordBoundaryPatch", and "SiblingSearch_Base"
//
// Parameter   :  SibID : sibling index (0~25)
//-------------------------------------------------------------------------------------------------------
int TABLE_04( const int SibID )
{

   switch ( SibID )
   {
      case 0: case 1: case 2: case 3: case 4: case 5:
         return 4;

      case 6: case 7: case 8: case 9: case 10: case 11:  case 12: case 13: case 14: case 15: case 16: case 17:
         return 2;

      case 18: case 19: case 20: case 21: case 22: case 23: case 24: case 25:
         return 1;

      default:
         fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "SibID", SibID );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         exit( -1 );
   }

}



