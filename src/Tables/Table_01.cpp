#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_01
// Description :  Return w0, w1, or w2 according to the inputs "SibIndex" and "dim" 
//
// Note        :  This table is particularly useful when the required values depended on the direction
//                of the sibling patch
//
// Parameter   :  SibIndex : Target sibling direction (0~25)
//                dim      : Target x/y/z direction
//                w0       : Value to be returned if SibIndex belongs to the left   surface
//                w1       : Value to be returned if SibIndex belongs to the middle surface
//                w2       : Value to be returned if SibIndex belongs to the right  surface
//
// Return      :  w0 || w1 || w2
//-------------------------------------------------------------------------------------------------------
int TABLE_01( const int SibIndex, const char dim, const int w0, const int w1, const int w2 )
{

   switch ( dim )
   {
      case 'x':
         switch ( SibIndex )
         {
            case 0: case 6: case 8: case 14: case 15: case 18: case 20: case 22: case 24:
               return w0;
               
            case 2: case 3: case 4: case 5: case 10: case 11: case 12: case 13:
               return w1;
              
            case 1: case 7: case 9: case 16: case 17: case 19: case 21: case 23: case 25:
               return w2;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n", 
                          "dim", dim, "SibIndex", SibIndex );
         }

      case 'y':
         switch ( SibIndex )
         {
            case 2: case 6: case 7: case 10: case 12: case 18: case 19: case 22: case 23:
               return w0;
               
            case 0: case 1: case 4: case 5: case 14: case 15: case 16: case 17:
               return w1;
                 
            case 3: case 8: case 9: case 11: case 13: case 20: case 21: case 24: case 25:
               return w2; 

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n", 
                          "dim", dim, "SibIndex", SibIndex );
         }

      case 'z':
         switch ( SibIndex )
         {
            case 4: case 10: case 11: case 14: case 16: case 18: case 19: case 20: case 21:
               return w0;

            case 0: case 1: case 2: case 3: case 6: case 7: case 8: case 9:
               return w1;

            case 5: case 12: case 13: case 15: case 17: case 22: case 23: case 24: case 25:
               return w2;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n", 
                          "dim", dim, "SibIndex", SibIndex );
         }

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c\" !!\n", "dim", dim );
         exit(1);

   } // switch ( dim )

} // FUNCTION : TABLE_01
