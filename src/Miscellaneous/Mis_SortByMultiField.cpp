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
// Parameter   :  NField    : Number of the fields
//                FieldSize : Size of a field
//                Array     : The array to be sort has size of [NField*FieldSize]. DO NOT change the array in this function!
//                IdxTable  : Index table to be returned
//                start_idx : The index start to sort
//                SortField : The field to be sorted
//                NSort     : The size of field to be sorted
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_SortByMultiField( const int NField, const long FieldSize, T **Array, long *IdxTable,
                           const long start_idx, const int SortField, const long NSort )
{
   if ( SortField == NField )
   {
      Aux_Message( stderr, "WARNING : Can not sort the exact same value.\n" );
      return;
   }

// check the inputs for the first call
   if ( SortField == 0 )
   {
      if ( IdxTable  == NULL  )   Aux_Error( ERROR_INFO, "NULL IdxTable.\n" );
      if ( FieldSize != NSort )   Aux_Error( ERROR_INFO, "FieldSize != NSort for the first call.\n" );
      if ( start_idx != 0     )   Aux_Error( ERROR_INFO, "start_idx != 0 for the first call.\n" );

      for (long i=0; i<FieldSize; i++)   IdxTable[i] = i;
   }

// 0. back up the field and the table
   T *Array_Sorted     = new T [NSort];
   for (long i=0; i<NSort; i++)       Array_Sorted[i]  = Array[SortField][IdxTable[start_idx+i]];
   long *IdxTable_copy = new long [FieldSize];
   for (long i=0; i<FieldSize; i++)   IdxTable_copy[i] = IdxTable[i];

// 1. sort by the field
   long *Idx_Sorted = new long [NSort];
   Mis_Heapsort( NSort, Array_Sorted, Idx_Sorted );

// 2. update the index table
   for (long i=0; i<NSort; i++)   IdxTable_copy[start_idx+i] = IdxTable[start_idx+Idx_Sorted[i]];

// 3. check the same value
   long NSameVal;
   for (long i=0; i<NSort-1; i++)
   {
      NSameVal = 1; // 1 --> itself

      while ( i+NSameVal < NSort  &&  Array_Sorted[i] == Array_Sorted[i+NSameVal] )   NSameVal++;

      if ( NSameVal == 1 )   continue;

      // 4. Sort the same values again
      Mis_SortByMultiField( NField, FieldSize, Array, IdxTable_copy, start_idx+i, SortField+1, NSameVal );

      i += NSameVal-1;
   } // for (long i=0; i<FieldSize-1; i++)

// 5. store the result
   for (long i=start_idx; i<start_idx+NSort; i++)   IdxTable[i] = IdxTable_copy[i];

   delete [] Array_Sorted;
   delete [] IdxTable_copy;
   delete [] Idx_Sorted;

} // FUNCTION :  Mis_SortByMultiField



// explicit template instantiation
template void Mis_SortByMultiField <int   > ( const int NField, const long FieldSize, int    **Array, long *IdxTable, const long start_idx, const int SortField, const long NSort );
template void Mis_SortByMultiField <long  > ( const int NField, const long FieldSize, long   **Array, long *IdxTable, const long start_idx, const int SortField, const long NSort );
template void Mis_SortByMultiField <ulong > ( const int NField, const long FieldSize, ulong  **Array, long *IdxTable, const long start_idx, const int SortField, const long NSort );
template void Mis_SortByMultiField <float > ( const int NField, const long FieldSize, float  **Array, long *IdxTable, const long start_idx, const int SortField, const long NSort );
template void Mis_SortByMultiField <double> ( const int NField, const long FieldSize, double **Array, long *IdxTable, const long start_idx, const int SortField, const long NSort );
