#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_AllocateArray2D
// Description :  Allocate a continuous 2D array with size [J][I]
//
// Note        :  1. Overloaded with different types
//                2. Call-by-reference
//                3. Free memory by Aux_DeallocateArray2D
//                4. Explicit template instantiation is put in the end of this file
//                5. One must be careful about using the Array pointer returned from this function,
//                   which is set to NULL if either I or J is 0
//                   --> Some functions may NOT accept a NULL pointer even when I or J is zero
//                   --> Moreover, Array[0...J-1] will become ill-defined when I == 0
//                       --> So accessing Array[0] will be illegal when I == 0 even though J > 0
//
// Parameter   :  Array : Pointer to be allocated
//                J/I   : Array dimensions
//-------------------------------------------------------------------------------------------------------
template <typename T>
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
// Function    :  Aux_AllocateArray3D
// Description :  Allocate a continuous 3D array with size [K][J][I]
//
// Note        :  1. Overloaded with different types
//                2. Call-by-reference
//                3. Free memory by Aux_DeallocateArray3D
//                4. Explicit template instantiation is put in the end of this file
//                5. One must be careful about using the Array pointer returned from this function,
//                   which is set to NULL if either I or J is 0
//                   --> Some functions may NOT accept a NULL pointer even when I or J is zero
//                   --> Moreover, Array[0...J-1] will become ill-defined when I == 0
//                       --> So accessing Array[0] will be illegal when I == 0 even though J > 0
//
// Parameter   :  Array : Pointer to be allocated
//                K/J/I : Array dimensions
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Aux_AllocateArray3D( T*** &Array, const int K, const int J, const int I )
{

   if ( K < 0  ||  J < 0  ||  I < 0 )  Aux_Error( ERROR_INFO, "incorrect array size (K = %d, J = %d, I = %d) !!\n", K, J, I );

   if ( K == 0  ||  J == 0  ||  I == 0 )
   {
      Array = NULL;
      return;
   }

   Array       = new T** [K    ];
   Array[0]    = new T*  [K*J  ];
   Array[0][0] = new T   [K*J*I];

   for (int k=1; k<K; k++)   Array[k]    = Array[0  ]    + k*J;

   for (int k=1; k<K; k++)   Array[k][0] = Array[k-1][0] + J*I;

   for (int k=0; k<K; k++)
   for (int j=1; j<J; j++)   Array[k][j] = Array[k  ][0] + j*I;

} // FUNCTION : Aux_AllocateArray3D



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeallocateArray2D
// Description :  Free an array previously allocated by Aux_AllocateArray2D
//
// Note        :  1. Overloaded with different types
//                2. Call-by-reference
//                3. Pointer is reset to NULL
//                4. Explicit template instantiation is put in the end of this file
//
// Parameter   :  Array : Pointer to be deallocated
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Aux_DeallocateArray2D( T** &Array )
{

// do NOT call delete if Array is NULL
   if ( Array == NULL )    return;

   delete [] Array[0];
   delete [] Array;

   Array = NULL;

} // FUNCTION : Aux_DeallocateArray2D



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_DeallocateArray3D
// Description :  Free an array previously allocated by Aux_AllocateArray3D
//
// Note        :  1. Overloaded with different types
//                2. Call-by-reference
//                3. Pointer is reset to NULL
//                4. Explicit template instantiation is put in the end of this file
//
// Parameter   :  Array : Pointer to be deallocated
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Aux_DeallocateArray3D( T*** &Array )
{

// do NOT call delete if Array is NULL
   if ( Array == NULL )    return;

   delete [] Array[0][0];
   delete [] Array[0];
   delete [] Array;

   Array = NULL;

} // FUNCTION : Aux_DeallocateArray3D



// explicit template instantiation
template void Aux_AllocateArray2D <float>  ( float**  &Array, const int J, const int I );
template void Aux_AllocateArray2D <double> ( double** &Array, const int J, const int I );
template void Aux_AllocateArray2D <int>    ( int**    &Array, const int J, const int I );
template void Aux_AllocateArray2D <long>   ( long**   &Array, const int J, const int I );

template void Aux_AllocateArray3D <float>  ( float***  &Array, const int K, const int J, const int I );
template void Aux_AllocateArray3D <double> ( double*** &Array, const int K, const int J, const int I );
template void Aux_AllocateArray3D <int>    ( int***    &Array, const int K, const int J, const int I );
template void Aux_AllocateArray3D <long>   ( long***   &Array, const int K, const int J, const int I );

template void Aux_DeallocateArray2D <float>  ( float**  &Array );
template void Aux_DeallocateArray2D <double> ( double** &Array );
template void Aux_DeallocateArray2D <int>    ( int**    &Array );
template void Aux_DeallocateArray2D <long>   ( long**   &Array );

template void Aux_DeallocateArray3D <float>  ( float***  &Array );
template void Aux_DeallocateArray3D <double> ( double*** &Array );
template void Aux_DeallocateArray3D <int>    ( int***    &Array );
template void Aux_DeallocateArray3D <long>   ( long***   &Array );

