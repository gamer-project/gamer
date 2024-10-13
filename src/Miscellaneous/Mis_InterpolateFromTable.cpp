#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFromTable
// Description :  Assuming y=y(x), return the interpolated value of y for a given point x
//
// Note        :  1. Interpolation table Table_x must be sorted into ascending numerical order in advance
//                2. Target point x must lie in the range Table_x[0] <= x < Table_x[N-1]
//                   --> Otherwise the function returns NULL_REAL
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put in the end of this file
//
// Parameter   :  N        : Number of elements in the interpolation tables Table_x and Table_y
//                           --> Must be >= 2
//                Table_x  : Interpolation table x
//                Table_y  : Interpolation table y
//                x        : Target point x for interpolation
//
// Return      :  y(x)      if x lies in the range Table_x[0] <= x < Table_x[N-1]
//                NULL_REAL if x lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFromTable( const int N, const T Table_x[], const T Table_y[], const T x )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( N <= 1 )           Aux_Error( ERROR_INFO, "incorrect input parameter \"N (%d) <= 1\" !!\n", N );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_y == NULL )  Aux_Error( ERROR_INFO, "Table_y == NULL !!\n" );
#  endif


// check whether the target x lies within the accepted range
   if ( x < Table_x[0]  ||  x >= Table_x[N-1] )    return NULL_REAL;


// binary search
   int IdxL, IdxR;
   T   xL, xR, yL, yR, y;

   IdxL = Mis_BinarySearch_Real( Table_x, 0, N-1, x );

#  ifdef GAMER_DEBUG
   if ( IdxL < 0  ||  IdxL >= N-1 )
      Aux_Error( ERROR_INFO, "IdxL (%d) lies outside the accepted range [%d ... %d] !!\n", IdxL, 0, N-2 );
#  endif

   IdxR = IdxL + 1;
   xL   = Table_x[IdxL];
   xR   = Table_x[IdxR];
   yL   = Table_y[IdxL];
   yR   = Table_y[IdxR];


// linear interpolation
   y = yL + (yR-yL)/(xR-xL)*(x-xL);

   return y;

} // FUNCTION : Mis_InterpolateFromTable



// explicit template instantiation
template float  Mis_InterpolateFromTable <float>  ( const int N, const float  Table_x[], const float  Table_y[], const float  x );
template double Mis_InterpolateFromTable <double> ( const int N, const double Table_x[], const double Table_y[], const double x );
