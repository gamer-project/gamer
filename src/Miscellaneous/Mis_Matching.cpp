#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Matching_char/int
// Description :  Use binary search (Mis_BinarySearch) to locate the positions of the integer numbers
//                in "Key" from the sorted integer array "Array"
//
// Note        :  1. "Array" and "Key" must be sorted in advance in ascending numerical order
//                2. If there are multiple elements matching Key, the return index can be any of them
//                3.    match --> record "1/array index" in the character/integer array "Match"
//                   no match --> record "0/-1" in the character/integer array "Match"
//                4. Overloaded with different types
//                   --> Explicit template instantiation is put in the end of this file
//                5. "Key" can contain multiple elements to be searched, and the minimum array index
//                   for searching (Min) will be set to the array index of the previous target in order
//                   to enhance the searching efficiency
//
// Return      :  Number of matching elements
//
// Parameter   :  N     : Size of the array "Array"
//                Array : Sorted look-up integer array
//                M     : Size of the array "Key"
//                Key   : Sorted integer array to search for
//                Match : Array storing the matching results
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// for char Match[]
//-------------------------------------------------------------------------------------------------------
template <typename T>
int Mis_Matching_char( const int N, const T Array[], const int M, const T Key[], char Match[] )
{

   int MatchIdx, Min = 0, NMatch = 0;

   for (int t=0; t<M; t++)
   {
      MatchIdx = Mis_BinarySearch( Array, Min, N-1, Key[t] );

      if ( MatchIdx != -1 )
      {
         Min      = MatchIdx;
         Match[t] = 1;
         NMatch ++;
      }
      else
         Match[t] = 0;
   }

   return NMatch;

} // FUNCTION : Mis_Matching_char



//-------------------------------------------------------------------------------------------------------
// for int Match[]
//-------------------------------------------------------------------------------------------------------
template <typename T>
int Mis_Matching_int( const int N, const T Array[], const int M, const T Key[], int Match[] )
{

   int Min = 0, NMatch = 0;

   for (int t=0; t<M; t++)
   {
      Match[t] = Mis_BinarySearch( Array, Min, N-1, Key[t] );

      if ( Match[t] != -1 )
      {
         Min = Match[t];
         NMatch ++;
      }
   }

   return NMatch;

} // FUNCTION : Mis_Matching_int



// explicit template instantiation
template int Mis_Matching_char <int>   ( const int N, const int   Array[], const int M, const int   Key[], char Match[] );
template int Mis_Matching_char <long>  ( const int N, const long  Array[], const int M, const long  Key[], char Match[] );
template int Mis_Matching_char <ulong> ( const int N, const ulong Array[], const int M, const ulong Key[], char Match[] );
template int Mis_Matching_int  <int>   ( const int N, const int   Array[], const int M, const int   Key[], int  Match[] );
template int Mis_Matching_int  <long>  ( const int N, const long  Array[], const int M, const long  Key[], int  Match[] );
template int Mis_Matching_int  <ulong> ( const int N, const ulong Array[], const int M, const ulong Key[], int  Match[] );
