#include "GAMER.h"




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
template <typename T>
int Mis_BinarySearch_Real( const T Array[], int Min, int Max, const T Key )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( Min < 0 )    Aux_Error( ERROR_INFO, "incorrect input parameter \"Min (%d) < 0\" !!\n", Min );
   if ( Max <= Min ) Aux_Error( ERROR_INFO, "incorrect input parameters \"Max (%d) <= Min (%d)\" !!\n", Max, Min );
#  endif


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


// check whether the found array index is correct
#  ifdef GAMER_DEBUG
   if ( Idx < Min  ||  Idx >= Max )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %d, Min %d, Max%d) !!\n", Idx, Min, Max );

   if (  Array[Idx] > Key  ||  Array[Idx+1] <= Key )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %d, ValueL %14.7e, ValueR %14.7e, Key %14.7e) !!\n",
                 Idx, Array[Idx], Array[Idx+1], Key );
#  endif


   return Idx;

} // FUNCTION : Mis_BinarySearch_Real



// explicit template instantiation
template int Mis_BinarySearch_Real <float>  ( const float  Array[], int Min, int Max, const float  Key );
template int Mis_BinarySearch_Real <double> ( const double Array[], int Min, int Max, const double Key );
