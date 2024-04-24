#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_SortByMultiField
// Description :  Sort the data by the input fields
//
// Note        :  1. Sorting by velocity may be necessary for STAR_FORMATION, where the new star particles
//                   created at different time but the same position may still have the same position for a
//                   while if velocity*dt is on the order of round-off errors
//                   --> Not supported yet since we may not have the velocity information (e.g., when adopting
//                       UseInputMassPos in Par_MassAssignment())
//                2. Invoked by Par_MassAssignment() and FB_AdvanceDt()
//
// Parameter   :  Array     : The array to be sorted. DO NOT change the array in this function!
//                IdxTable  : Index table to be returned has a size of NSort.
//                NSort     : The size of field to be sorted.
//                SortOrder : The order of sort field has a size of NOrder.
//                NOrder    : The size of the sort field order.
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_SortByMultiField( T **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder )
{
   T    *Array_Sorted = new T    [NSort];
   long *Idx_Sorted   = new long [NSort];
   long *IdxTable_old = new long [NSort];

// 0. back up the field
   for (long i=0; i<NSort; i++)
   {
      Array_Sorted[i] = Array[SortOrder[0]][IdxTable[i]];
      IdxTable_old[i] = IdxTable[i];
   }

// 1. sort by the field
   Mis_Heapsort( NSort, Array_Sorted, Idx_Sorted );

// 2. update the index table
   for (long i=0; i<NSort; i++)   IdxTable[i] = IdxTable_old[Idx_Sorted[i]];

// 3. check the same value
   long NSameVal;
   for (long i=0; i<NSort-1; i++)
   {
      NSameVal = 1; // 1 --> itself

      while ( i+NSameVal < NSort  &&  Array_Sorted[i] == Array_Sorted[i+NSameVal] )   NSameVal++;

      if ( NSameVal == 1 )   continue;
      if ( NOrder == 1 )
      {
         Aux_Message( stderr, "WARNING : Can not sort the exact same value.\n" );
         break;
      }

//    4. Sort the same values again
      long *IdxTable_same  = new long [NSameVal];
      for (long j=0; j<NSameVal; j++)   IdxTable_same[j]  = IdxTable[i+j];

      Mis_SortByMultiField( Array, IdxTable_same, NSameVal, SortOrder+1, NOrder-1 );

//    5. store the result
      for (long j=0; j<NSameVal; j++)   IdxTable[i+j] = IdxTable_same[j];

      delete [] IdxTable_same;

      i += NSameVal-1;
   } // for (long i=0; i<FieldSize-1; i++)

   delete [] Array_Sorted;
   delete [] IdxTable_old;
   delete [] Idx_Sorted;

} // FUNCTION :  Mis_SortByMultiField



// explicit template instantiation
template void Mis_SortByMultiField <int   > ( int    **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByMultiField <long  > ( long   **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByMultiField <ulong > ( ulong  **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByMultiField <float > ( float  **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
template void Mis_SortByMultiField <double> ( double **Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
