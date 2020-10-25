#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  TABLE_02
// Description :  Return w0, or w1 according to the inputs "LocalID" and "dim"
//
// Note        :  This table is particularly useful when the required values depended on the direction
//                of the local ID of patch within a patch group
//
// Parameter   :  LocalID : Local index within the patch group (0~7)
//                dim     : Target x/y/z direction
//                w0      : Value to be returned if LocalID belongs to the left   surface
//                w1      : Value to be returned if LocalID belongs to the middle surface
//                w2      : Value to be returned if LocalID belongs to the right  surface
//
// Return      :  w0 || w1
//-------------------------------------------------------------------------------------------------------
template <typename T>
T TABLE_02( const int LocalID, const char dim, const T w0, const T w1 )
{

   switch ( dim )
   {
      case 'x':
         switch ( LocalID )
         {
            case 0: case 2: case 3: case 5:
               return w0;

            case 1: case 4: case 6: case 7:
               return w1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n",
                          "dim", dim, "LocalID", LocalID );
         }

      case 'y':
         switch ( LocalID )
         {
            case 0: case 1: case 3: case 6:
               return w0;

            case 2: case 4: case 5: case 7:
               return w1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n",
                          "dim", dim, "LocalID", LocalID );
         }

      case 'z':
         switch ( LocalID )
         {
            case 0: case 1: case 2: case 4:
               return w0;

            case 3: case 5: case 6: case 7:
               return w1;

            default:
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %c, %s = %d\" !!\n",
                          "dim", dim, "LocalID", LocalID );
         }

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %c\" !!\n", "dim", dim );
         exit(1);

   } // switch ( dim )

} // FUNCTION : TABLE_02


// explicit template instantiation
template float  TABLE_02 <float > ( const int, const char, const float,  const float  );
template double TABLE_02 <double> ( const int, const char, const double, const double );
template int    TABLE_02 <int   > ( const int, const char, const int,    const int    );
template long   TABLE_02 <long  > ( const int, const char, const long,   const long   );
template bool   TABLE_02 <bool  > ( const int, const char, const bool,   const bool   );
