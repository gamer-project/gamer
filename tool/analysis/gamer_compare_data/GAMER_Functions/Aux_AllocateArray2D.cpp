#include "CompareData.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_AllocateArray2D
// Description :  Allocate a continuous 2D array with size [J][I]
//
// Note        :  1. Call-by-reference
//                2. Free memory by Aux_DeallocateArray2D
//
// Parameter   :  Array : Pointer to be allocated
//                J/I   : Array dimensions
//-------------------------------------------------------------------------------------------------------
template<typename T>
void Aux_AllocateArray2D( T** &Array, const int J, const int I )
{

   if ( J < 0  ||  I < 0 )    Aux_Error( ERROR_INFO, "incorrect array size (J = %d, I = %d) !!\n", J, I );

   if ( J == 0  ||  I == 0 )
   {
      Array = NULL;
      return;
   }

   Array    = new T* [J  ];
   Array[0] = new T  [J*I];

   for (int j=1; j<J; j++)    Array[j] = Array[j-1] + I;

} // FUNCTION : Aux_AllocateArray2D



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeallocateArray2D
// Description :  Free an array previously allocated by Aux_AllocateArray2D
//
// Note        :  1. Call-by-reference
//                2. Pointer is reset to NULL
//
// Parameter   :  Array : Pointer to be deallocated
//-------------------------------------------------------------------------------------------------------
template<typename T>
void Aux_DeallocateArray2D( T** &Array )
{

   if ( Array == NULL )    return;

   delete [] Array[0];
   delete [] Array;

   Array = NULL;

} // FUNCTION : Aux_DeallocateArray2D


// explicit template instantiation
template void Aux_AllocateArray2D <int>    ( int**    &Array, const int J, const int I );
template void Aux_AllocateArray2D <long>   ( long**   &Array, const int J, const int I );
template void Aux_AllocateArray2D <float>  ( float**  &Array, const int J, const int I );
template void Aux_AllocateArray2D <double> ( double** &Array, const int J, const int I );

template void Aux_DeallocateArray2D <int>    ( int**    &Array );
template void Aux_DeallocateArray2D <long>   ( long**   &Array );
template void Aux_DeallocateArray2D <float>  ( float**  &Array );
template void Aux_DeallocateArray2D <double> ( double** &Array );
