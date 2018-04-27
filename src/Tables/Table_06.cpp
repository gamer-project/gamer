#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_06
// Description :  Return the sibling index used by Flag_Buffer() and Buf_RecordBoundaryFlag()
//
// Parameter   :  SibID     : Sibling index (0~25)
//                FlagLayer : Flag layer (1~8)
//-------------------------------------------------------------------------------------------------------
int TABLE_06( const int SibID, const int FlagLayer )
{

   switch ( SibID )
   {
      case 0: case 1: case 2: case 3: case 4: case 5:
      {
         switch ( FlagLayer )
         {
            case 1:  return SibID;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 6:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 2;
            case 3:  return 6;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 7:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 2;
            case 3:  return 7;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 8:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 3;
            case 3:  return 8;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 9:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 3;
            case 3:  return 9;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 10:
      {
         switch ( FlagLayer )
         {
            case 1:  return 2;
            case 2:  return 4;
            case 3:  return 10;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 11:
      {
         switch ( FlagLayer )
         {
            case 1:  return 3;
            case 2:  return 4;
            case 3:  return 11;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 12:
      {
         switch ( FlagLayer )
         {
            case 1:  return 2;
            case 2:  return 5;
            case 3:  return 12;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 13:
      {
         switch ( FlagLayer )
         {
            case 1:  return 3;
            case 2:  return 5;
            case 3:  return 13;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 14:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 4;
            case 3:  return 14;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 15:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 5;
            case 3:  return 15;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 16:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 4;
            case 3:  return 16;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 17:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 5;
            case 3:  return 17;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 18:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 2;
            case 3:  return 4;
            case 4:  return 6;
            case 5:  return 10;
            case 6:  return 14;
            case 7:  return 18;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 19:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 2;
            case 3:  return 4;
            case 4:  return 7;
            case 5:  return 10;
            case 6:  return 16;
            case 7:  return 19;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 20:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 3;
            case 3:  return 4;
            case 4:  return 8;
            case 5:  return 11;
            case 6:  return 14;
            case 7:  return 20;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 21:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 3;
            case 3:  return 4;
            case 4:  return 9;
            case 5:  return 11;
            case 6:  return 16;
            case 7:  return 21;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 22:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 2;
            case 3:  return 5;
            case 4:  return 6;
            case 5:  return 12;
            case 6:  return 15;
            case 7:  return 22;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 23:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 2;
            case 3:  return 5;
            case 4:  return 7;
            case 5:  return 12;
            case 6:  return 17;
            case 7:  return 23;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 24:
      {
         switch ( FlagLayer )
         {
            case 1:  return 0;
            case 2:  return 3;
            case 3:  return 5;
            case 4:  return 8;
            case 5:  return 13;
            case 6:  return 15;
            case 7:  return 24;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      case 25:
      {
         switch ( FlagLayer )
         {
            case 1:  return 1;
            case 2:  return 3;
            case 3:  return 5;
            case 4:  return 9;
            case 5:  return 13;
            case 6:  return 17;
            case 7:  return 25;
            default: Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "FlagLayer", FlagLayer );
         }
      }

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d\" !!\n", "SibID", SibID );
         exit( EXIT_FAILURE );

   } // switch ( SibID )

} // FUNCTION : TABLE_06
