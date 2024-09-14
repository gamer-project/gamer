#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_BinarySearch_Real
// Description :  Use binary search to find the proper array index "Idx" in the input "Array" satisfying
//
//                   Array[Idx] <= Key < Array[Idx+1]
//
// Note        :  1. "Array" must be sorted in advance in ascending numerical order
//                2. If there are multiple elements equal to Key, the return index will be the maximum index among them
//                3. Typical returned index should be within the range "Min <= Idx < Max", However,
//                   if Key <  Array[Min] (minimum value) --> return Min-1
//                   if Key >= Array[Max] (maximum value) --> return Max
//                4. We must have "Max > Min" in the input parameters
//                5. Overloaded with different types
//                6. Explicit template instantiation is put in the end of this file
//
// Parameter   :  Array : Sorted look-up array (in ascending numerical order)
//                Min   : Minimum array index for searching
//                Max   : Maximum array index for searching
//                Key   : Target value to search for
//
// Return      :  Idx   : if target is found
//                Min-1 : if Key <  Array[Min]
//                Max   : if Key >= Array[Max]
//-------------------------------------------------------------------------------------------------------
int BinarySearch( const real Array[], int Min, int Max, const real Key )
{

// check whether the input key lies outside the target range
   if ( Key <  Array[Min] )   return Min-1;
   if ( Key >= Array[Max] )   return Max;


// binary search
   int Idx = -2;

   while (  ( Idx=(Min+Max)/2 ) != Min  )
   {
      if   ( Array[Idx] > Key )  Max = Idx;
      else                       Min = Idx;
   }

   return Idx;

} // FUNCTION : Mis_BinarySearch_Real



int double_BinarySearch( const double Array[], int Min, int Max, const real Key )
{

// check whether the input key lies outside the target range
   if ( Key <  Array[Min] )   return Min-1;
   if ( Key >= Array[Max] )   return Max;


// binary search
   int Idx = -2;

   while (  ( Idx=(Min+Max)/2 ) != Min  )
   {
      if   ( Array[Idx] > Key )  Max = Idx;
      else                       Min = Idx;
   }

   return Idx;

} // FUNCTION : Mis_BinarySearch_Real
