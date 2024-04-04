#include "GAMER.h"


template <typename T>

//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_SortByMultiField
// Description :  Sort the data by the input fields
//
// Note        :  1.
//                2.
//                3. Invoked by Par_MassAssignment() and FB_AdvanceDt()
//
// Parameter   :
//
// Return      :  IdxTable
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_SortByMultiField( const int NField, const long FieldSize, const T **Array, long *IdxTable,
                           const long start_idx, const int SortField, const long NSort )
{
   if ( SortField == NField-1 )
   {
//    TODO : print the warning message if the data is exactally the same
      printf("WARNING : Can not sort the exact same value.\n");
      return;
   }

// check the size of the **Array and IdxTable and FieldToBeSorted
   if ( SortField == 0 ) {
       if (IdxTable == NULL)
       {
           printf("NULL IdxTable\n");
           return;
       }
       if (FieldSize != NSort)
       {
           printf("FieldSize != NSort\n");
           return;
       }
       if (start_idx != 0)
       {
           printf("start_idx != 0\n");
           return;
       }
       for (long i=0; i<FieldSize; i++) IdxTable[i] = i;
   }

// 0. back up the field
   T *Array_Sorted     = new T [NSort];
   for (long i=0; i<NSort; i++)
   {
       Array_Sorted[i]  = Array[SortField][start_idx+i];
   }
   long *IdxTable_copy = new long [NSort];
   for (long i=0; i<FieldSize; i++)
   {
       IdxTable_copy[i] = IdxTable[i];
   }
   
// 1. sort by the field
   long *Idx_Sorted = new long [NSort];
   Mis_Heapsort( NSort, Array_Sorted, Idx_Sorted );

// 2. update the index table
   for (long i=0; i<NSort; i++)
   {
      IdxTable_copy[start_idx+i] = IdxTable[start_idx+Idx_Sorted[i]];
   }

// 3. check the same value
   long NSameVal;
   for (long i=0; i<NSort-1; i++)
   {
       NSameVal = 1; // 1 --> itself

       while ( i+NSameVal < NSort  &&  Array_Sorted[i] == Array_Sorted[i+NSameVal] )   NSameVal ++;

       if ( NSameVal == 1 ) continue;

       Mis_SortByMultiField( NField, FieldSize, Array, IdxTable_copy, start_idx+i, SortField+1, NSameVal );

       i += NSameVal - 1;
   } // for (long i=0; i<FieldSize-1; i++)

// 4. store the result
   for (long i=0; i<NSort; i++)
   {
      IdxTable[i] = IdxTable_copy[i];
   }

   delete [] Array_Sorted;
   delete [] IdxTable_copy;
   delete [] Idx_Sorted;


} // FUNCTION :  Mis_SortByMultiField



// explicit template instantiation
template void Mis_SortByMultiField <int>    ( const int NField, const long FieldSize, const int    **Array );
template void Mis_SortByMultiField <long>   ( const int NField, const long FieldSize, const long   **Array );
template void Mis_SortByMultiField <ulong>  ( const int NField, const long FieldSize, const ulong  **Array );
template void Mis_SortByMultiField <float>  ( const int NField, const long FieldSize, const float  **Array );
template void Mis_SortByMultiField <double> ( const int NField, const long FieldSize, const double **Array );
