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
template <typename T>
T TABLE_01( const int SibIndex, const char dim, const T w0, const T w1, const T w2 )
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


// explicit template instantiation
template float  TABLE_01 <float > ( const int, const char, const float,  const float,  const float  );
template double TABLE_01 <double> ( const int, const char, const double, const double, const double );
template int    TABLE_01 <int   > ( const int, const char, const int,    const int,    const int    );
template long   TABLE_01 <long  > ( const int, const char, const long,   const long,   const long   );
template bool   TABLE_01 <bool  > ( const int, const char, const bool,   const bool,   const bool   );
