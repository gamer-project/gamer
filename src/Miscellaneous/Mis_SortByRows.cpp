#include "GAMER.h"

template <typename T>
static void SortByRows( T const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_SortByRows
// Description :  Indirectly sort the columns of a 2D array using a sequence of sorting rows.
//
// Note        :  1. Invoked by Par_MassAssignment() and FB_AdvanceDt()
//                2. Example:
//
//                   int  NOrder         = 2;
//                   int  SortOrder[2]   = { 0, 1 }; // Sort by Array[0] first, then by Array[1]
//                   long NSort          = 7;
//                   int Array[2][NSort] = { { 1, 5, 1, 4, 3, 4, 4 },
//                                           { 9, 4, 0, 4, 0, 2, 1 } };
//                   long IdxTable[NSort];
//
//                   Mis_SortByRows( Array, IdxTable, NSort, SortOrder, NOrder );
//                   // => SortIdxTable becomes { 2, 0, 4, 6, 5, 3, 1 };
//
//                2. It serves as a gate function between the local function SortByRows() and external functions.
//                   --> It's sole purpose is to initialize IdxTable[]
//
// Parameter   :  Array     : A 2D m by n array to be sorted, Array[m][n], which can be interpreted as m rows of 1D arrays with n columns.
//                            The array will not be changed in this function since it is an indirect sort.
//                IdxTable  : Table of indices of the columns to be sorted, with a size of NSort.
//                            This will be changed and returned to sort the array.
//                NSort     : Number of columns (n) to be sorted in Array[m][n].
//                SortOrder : Array with a size of NOrder specifies which rows to compare first, second, etc.
//                NOrder    : Size of SortOrder (must be <= m).
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_SortByRows( T const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder )
{

//  initialize the IdxTable
    for (long i=0; i<NSort; i++)   IdxTable[i] = i;

    SortByRows<T>( Array, IdxTable, NSort, SortOrder, NOrder );

} // FUNCTION : Mis_SortByRow



//-------------------------------------------------------------------------------------------------------
// Function    :  SortByRows
// Description :  Local function called by Mis_SortByRow()
//
// Note        :  1. Recursive function
//
// Parameter   :  See Mis_SortByRow()
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
template <typename T>
void SortByRows( T const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder )
{

// check
   if ( NSort  < 0L )   Aux_Error( ERROR_INFO, "NSort < 0 !!\n" );
   if ( NOrder < 1  )   Aux_Error( ERROR_INFO, "NOrder < 1 !!\n" );


   T    *Array_Sorted = new T    [NSort];
   long *Idx_Sorted   = new long [NSort];
   long *IdxTable_old = new long [NSort];

// 0. back up the data to be sorted
   for (long i=0; i<NSort; i++)
   {
      Array_Sorted[i] = Array[ SortOrder[0] ][ IdxTable[i] ];
      IdxTable_old[i] = IdxTable[i];
   }

// 1. sort the data
   Mis_Heapsort( NSort, Array_Sorted, Idx_Sorted );

// 2. update the index table
   for (long i=0; i<NSort; i++)   IdxTable[i] = IdxTable_old[ Idx_Sorted[i] ];

// 3. check the same value
   for (long i=0; i<NSort-1L; i++)
   {
      long NSameVal = 1L; // 1 --> itself

      while ( i+NSameVal < NSort  &&  Array_Sorted[i] == Array_Sorted[i+NSameVal] )    NSameVal++;

      if ( NSameVal == 1L )   continue;
      else if ( NOrder == 1 )
      {
         Aux_Message( stderr, "WARNING : Cannot sort the exact same value !!\n" );
         break;
      }

//    4. sort the same values by the next row
      long *IdxTable_same = new long [NSameVal];
      for (long j=0; j<NSameVal; j++)   IdxTable_same[j] = IdxTable[i+j];

      SortByRows<T>( Array, IdxTable_same, NSameVal, SortOrder+1, NOrder-1 );

//    5. store the result
      for (long j=0; j<NSameVal; j++)   IdxTable[i+j] = IdxTable_same[j];

      delete [] IdxTable_same;

      i += NSameVal-1L;
   } // for (long i=0; i<NSort-1L; i++)

   delete [] Array_Sorted;
   delete [] IdxTable_old;
   delete [] Idx_Sorted;

} // FUNCTION : SortByRows



// explicit template instantiation
template void Mis_SortByRows <int   > ( int    const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByRows <long  > ( long   const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByRows <ulong > ( ulong  const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByRows <float > ( float  const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByRows <double> ( double const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
