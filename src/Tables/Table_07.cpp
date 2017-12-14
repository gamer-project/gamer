#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_07
// Description :  Return the local patch IDs near the "SibID" boundary
//
// Note        :  a. Must follow the order : x -> y -> z 
//                   (otherwise the function "Buf_FindBoundaryPatch" will fail)
//                b. Similar to the table "TABLE_03"
//                c. Work for the functions "Buf_RecordExchangeDataPatchID" and "Buf_RecordBoundaryPatch"
//
// Parameter   :  SibID : Sibling index (0~25)
//                Count : Patch counter (0~3)
//-------------------------------------------------------------------------------------------------------
int TABLE_07( const int SibID, const int Count )
{

   switch ( SibID )
   {
      case 0:
      {
         switch( Count )
         {
            case 0: return 0;
            case 1: return 2;
            case 2: return 3;
            case 3: return 5;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 1:
      {
         switch ( Count )
         {
            case 0: return 1;
            case 1: return 4;
            case 2: return 6;
            case 3: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 2:
      {
         switch ( Count )
         {
            case 0: return 0;
            case 1: return 1;
            case 2: return 3;
            case 3: return 6;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 3:
      {
         switch ( Count )
         {
            case 0: return 2;
            case 1: return 4;
            case 2: return 5;
            case 3: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 4:
      {
         switch ( Count )
         {
            case 0: return 0;
            case 1: return 1;
            case 2: return 2;
            case 3: return 4;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 5:
      {
         switch ( Count )
         {
            case 0: return 3;
            case 1: return 6;
            case 2: return 5;
            case 3: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 6:
      {
         switch ( Count )
         {
            case 0: return 0;
            case 1: return 3;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 7:
      {
         switch ( Count )
         {
            case 0: return 1;
            case 1: return 6;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 8:
      {
         switch ( Count )
         {
            case 0: return 2;
            case 1: return 5;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 9:
      {
         switch ( Count ) 
         {
            case 0: return 4;
            case 1: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 10:
      {
         switch ( Count )
         {
            case 0: return 0;
            case 1: return 1;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 11:
      {
         switch ( Count )
         {
            case 0: return 2;
            case 1: return 4;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }
         
      case 12:
      {
         switch ( Count )
         {
            case 0: return 3;
            case 1: return 6;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 13:
      {
         switch ( Count )
         {
            case 0: return 5;
            case 1: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 14:
      {
         switch ( Count )
         {
            case 0: return 0;
            case 1: return 2;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 15:
      {
         switch ( Count )
         {
            case 0: return 3;
            case 1: return 5;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 16:
      {
         switch ( Count )
         {
            case 0: return 1;
            case 1: return 4;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 17:
      {
         switch ( Count )
         {
            case 0: return 6;
            case 1: return 7;
            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d, %s = %d\" !!\n", 
                          "SibID", SibID, "Count", Count );
         }
      }

      case 18:    return 0;
      case 19:    return 1;
      case 20:    return 2;
      case 21:    return 4;
      case 22:    return 3;
      case 23:    return 6;
      case 24:    return 5;
      case 25:    return 7;

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "SibID", SibID );
         exit(1);

   } // switch ( SibID )

} // FUNCTION : TABLE_07
