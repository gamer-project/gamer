#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_BinarySearch
// Description :  Use binary search to locate the position of an input number in a sorted array
//
// Note        :  1. "Array" must be sorted in advance in ascending numerical order
//                2. If there are multiple elements matching Key, the return index can be any of them
//                3. An overloaded function for the "long" Array and Key is also created
//                4. Overloaded with different types
//                   --> Explicit template instantiation is put in the end of this file
//
// Return      :     match --> array index
//                no match --> -1
//
// Parameter   :  Array : Sorted look-up integer array (in ascending numerical order)
//                Min   : Minimum array index for searching
//                Max   : Maximum array index for searching
//                Key   : Integer number to search for
//-------------------------------------------------------------------------------------------------------
template <typename U, typename T>
U Mis_BinarySearch( const T Array[], U Min, U Max, const T Key )
{

   U Mid = -1;

   while ( Min <= Max )
   {
      Mid = ( Min + Max ) / 2;

      if      ( Array[Mid] > Key )  Max = Mid-1;
      else if ( Array[Mid] < Key )  Min = Mid+1;
      else                          return Mid;
   }

   return -1;

} // FUNCTION : Mis_BinarySearch



// explicit template instantiation
template int  Mis_BinarySearch <int,int>    ( const int   Array[], int  Min, int  Max, const int   Key );
template int  Mis_BinarySearch <int,long>   ( const long  Array[], int  Min, int  Max, const long  Key );
template int  Mis_BinarySearch <int,ulong>  ( const ulong Array[], int  Min, int  Max, const ulong Key );

template long Mis_BinarySearch <long,int>   ( const int   Array[], long Min, long Max, const int   Key );
template long Mis_BinarySearch <long,long>  ( const long  Array[], long Min, long Max, const long  Key );
template long Mis_BinarySearch <long,ulong> ( const ulong Array[], long Min, long Max, const ulong Key );
