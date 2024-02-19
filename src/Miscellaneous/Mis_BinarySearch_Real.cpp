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
template <typename U, typename T>
U Mis_BinarySearch_Real( const T Array[], U Min, U Max, const T Key )
{

// initial check; force print format as long since typename U can be "int" or "long"
#  ifdef GAMER_DEBUG
   if ( Min < 0 )    Aux_Error( ERROR_INFO, "incorrect input parameter \"Min (%ld) < 0\" !!\n", (long)Min );
   if ( Max <= Min ) Aux_Error( ERROR_INFO, "incorrect input parameters \"Max (%ld) <= Min (%ld)\" !!\n", (long)Max, (long)Min );
#  endif


// check whether the input key lies outside the target range
   if ( Key <  Array[Min] )   return Min-1;
   if ( Key >= Array[Max] )   return Max;


// binary search
   U Idx = -2;

   while (  ( Idx=(Min+Max)/2 ) != Min  )
   {
      if   ( Array[Idx] > Key )  Max = Idx;
      else                       Min = Idx;
   }


// check whether the found array index is correct; force print format as long since typename U can be "int" or "long"
#  ifdef GAMER_DEBUG
   if ( Idx < Min  ||  Idx >= Max )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %ld, Min %ld, Max%ld) !!\n", (long)Idx, (long)Min, (long)Max );

   if (  Array[Idx] > Key  ||  Array[Idx+1] <= Key )
      Aux_Error( ERROR_INFO, "incorrect output index (Idx %ld, ValueL %14.7e, ValueR %14.7e, Key %14.7e) !!\n",
                 (long)Idx, Array[Idx], Array[Idx+1], Key );
#  endif


   return Idx;

} // FUNCTION : Mis_BinarySearch_Real



// explicit template instantiation
template int  Mis_BinarySearch_Real <int,float>   ( const float  Array[], int  Min, int  Max, const float  Key );
template int  Mis_BinarySearch_Real <int,double>  ( const double Array[], int  Min, int  Max, const double Key );

template long Mis_BinarySearch_Real <long,float>  ( const float  Array[], long Min, long Max, const float  Key );
template long Mis_BinarySearch_Real <long,double> ( const double Array[], long Min, long Max, const double Key );
