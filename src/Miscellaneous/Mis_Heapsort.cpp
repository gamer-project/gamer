#include "GAMER.h"

template <typename U, typename T>
static void Heapsort_SiftDown( const U L, const U R, T Array[], U IdxTable[] );




//OPTIMIZATION : (1) quick sort  (2) try the "qsort" library
//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Heapsort
// Description :  Use the Heapsort algorithm to sort the input array into ascending numerical order
//                --> An index table will also be constructed if "IdxTable != NULL"
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//                2. Overloaded with different types
//                   --> Explicit template instantiation is put in the end of this file
//
// Parameter   :  N        :  Size of Array
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
template <typename U, typename T>
void Mis_Heapsort( const U N, T Array[], U IdxTable[] )
{

// initialize the IdxTable
   if ( IdxTable != NULL )
      for (U t=0; t<N; t++)    IdxTable[t] = t;

// heap creation
   for (U L=N/2-1; L>=0; L--)  Heapsort_SiftDown<U,T>( L, N-1, Array, IdxTable );

// retirement-and-promotion
   T Buf;
   for (U R=N-1; R>0; R--)
   {
      Buf      = Array[R];
      Array[R] = Array[0];
      Array[0] = Buf;

      if ( IdxTable != NULL )
      {
         Buf         = IdxTable[R];
         IdxTable[R] = IdxTable[0];
         IdxTable[0] = Buf;
      }

      Heapsort_SiftDown<U,T>( 0, R-1, Array, IdxTable );
   }

} // FUNCTION : Mis_Heapsort



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort_SiftDown
// Description :  Sift-down process for the Heapsort algorithm
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//                2. Overloaded with different types
//
// Parameter   :  L        :  Left  range of the sift-down
//                R        :  Right range of the sift-down
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
template <typename U, typename T>
void Heapsort_SiftDown( const U L, const U R, T Array[], U IdxTable[] )
{

   U  Idx_up    = L;
   U  Idx_down  = 2*Idx_up + 1;
   T  Target    = Array[Idx_up];
   U  TargetIdx = ( IdxTable == NULL ) ? -1 : IdxTable[Idx_up];

   while ( Idx_down <= R )
   {
//    find the better employee
      if ( Idx_down < R  &&  Array[Idx_down+1] > Array[Idx_down] )   Idx_down ++;

//    terminate the sift-down process if the target (supervisor) is better than both its employees
      if ( Target >= Array[Idx_down] )    break;

//    otherwise, promote the better employee
      Array[Idx_up] = Array[Idx_down];
      if ( IdxTable != NULL )    IdxTable[Idx_up] = IdxTable[Idx_down];

//    prepare the next sift-down operation
      Idx_up   = Idx_down;
      Idx_down = 2*Idx_up + 1;
   }

// put target at its best position
   Array[Idx_up] = Target;
   if ( IdxTable != NULL )    IdxTable[Idx_up] = TargetIdx;

} // FUNCTION : Heapsort_SiftDown



// explicit template instantiation
template void Mis_Heapsort <int,int>     ( const int  N, int    Array[], int  IdxTable[] );
template void Mis_Heapsort <int,long>    ( const int  N, long   Array[], int  IdxTable[] );
template void Mis_Heapsort <int,ulong>   ( const int  N, ulong  Array[], int  IdxTable[] );
template void Mis_Heapsort <int,float>   ( const int  N, float  Array[], int  IdxTable[] );
template void Mis_Heapsort <int,double>  ( const int  N, double Array[], int  IdxTable[] );

template void Mis_Heapsort <long,int>    ( const long N, int    Array[], long IdxTable[] );
template void Mis_Heapsort <long,long>   ( const long N, long   Array[], long IdxTable[] );
template void Mis_Heapsort <long,ulong>  ( const long N, ulong  Array[], long IdxTable[] );
template void Mis_Heapsort <long,float>  ( const long N, float  Array[], long IdxTable[] );
template void Mis_Heapsort <long,double> ( const long N, double Array[], long IdxTable[] );

