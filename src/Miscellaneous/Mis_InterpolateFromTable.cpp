#include "GAMER.h"

template <typename T>
static int GetIdxL_From1DTable( const int N, const T Table_x[], const T x );
static int GetIdx_corner_nDim_Table( const int Idx_corner_local, const int nDim, const int N_x[], const int IdxL[] );



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
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_y == NULL )  Aux_Error( ERROR_INFO, "Table_y == NULL !!\n" );
#  endif


// get index
   int IdxL = GetIdxL_From1DTable( N, Table_x, x );

   if ( IdxL == NULL_INT )    return NULL_REAL;


// linear interpolation
   T y = Mis_LinearInterpolate( Table_x[IdxL], Table_x[IdxL+1], Table_y[IdxL], Table_y[IdxL+1], x );

   return y;

} // FUNCTION : Mis_InterpolateFromTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom_nDim_Table
// Description :  Assuming f=f(\vec{x}), return the interpolated value of f for a given point
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in N-dimension space
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. Target coordinate x_d in each dimension must lie in the range Table_x[d][0] <= x_d < Table_x[d][N_x[d]-1]
//                   --> Otherwise the function returns NULL_REAL
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put in the end of this file
//
// Parameter   :  nDim     : Number of dimensions of the Table_f
//                N_x      : Number of elements in the interpolation tables of each dimension Table_x[d]
//                           --> it is an array with size = nDim
//                           --> Must be >= 2 for each dimension
//                Table_x  : Interpolation table x of each dimension
//                           --> it is an array of arrays
//                           --> Table_x[d] has size = N_x[d], where d is the dimension from 0 to nDim-1
//                Table_f  : Interpolation table f
//                           --> it is an array with size = N_x[0]*N_x[1]*...*N_x[nDim-1]
//                x        : Target point \vec{x} for interpolation
//                           --> it is an array with size = nDim
//
// Return      :  f(\vec{x}) if \vec{x} lies in the range Table_x[d][0] <= x[d] < Table_x[d][N_x[d]-1] for d from 0 to nDim-1
//                NULL_REAL  if \vec{x} lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom_nDim_Table( const int nDim, const int N_x[], T const* const* Table_x, const T Table_f[], const T x[] )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim <=0 !!\n" );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
   if ( x == NULL )        Aux_Error( ERROR_INFO, "x == NULL !!\n" );
   for (int d=0; d<nDim; d++)
      if ( Table_x[d] == NULL )  Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
#  endif


// get index
   int *IdxL = new int [nDim];
   for (int d=0; d<nDim; d++)
   {
      IdxL[d] = GetIdxL_From1DTable( N_x[d], Table_x[d], x[d] );

      if ( IdxL[d] == NULL_INT )    return NULL_REAL;
   }

// linear interpolation
   T y = Mis_InterpolateFrom_nDim_Table_withIdxL( nDim, N_x, Table_x, Table_f, x, IdxL );

   delete [] IdxL;

   return y;

} // FUNCTION : Mis_InterpolateFrom_nDim_Table




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom_nDim_Table_withIdxL
// Description :  Assuming f=f(\vec{x}), return the interpolated value of f for a given point
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in N-dimension space
//                with the index on the left in the table given
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. Target coordinate x_d in each dimension must lie in the range Table_x[d][IdxL[d]] <= x_d < Table_x[d][IdxL[d]+1]
//                   --> Otherwise the function returns NULL_REAL
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put in the end of this file
//
// Parameter   :  nDim     : Number of dimensions of the Table_f
//                N_x      : Number of elements in the interpolation tables of each dimension Table_x[d]
//                           --> it is an array with size = nDim
//                           --> Must be >= 2 for each dimension
//                Table_x  : Interpolation table x of each dimension
//                           --> it is an array of arrays
//                           --> Table_x[d] has size = N_x[d], where d is the dimension from 0 to nDim-1
//                Table_f  : Interpolation table f
//                           --> it is an array with size = N_x[0]*N_x[1]*...*N_x[nDim-1]
//                x        : Target point \vec{x} for interpolation
//                           --> it is an array with size = nDim
//                IdxL     : Index of Table_x which Table_x[d][IdxL[d]] lies on the left side of x[d] for each dimension
//                           --> it is an array with size = nDim
//
// Return      :  f(\vec{x}) if \vec{x} lies in the range Table_x[d][IdxL[d]] <= x[d] < Table_x[d][IdxL[d]+1] for d from 0 to nDim-1
//                NULL_REAL  if \vec{x} lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom_nDim_Table_withIdxL( const int nDim, const int N_x[], T const* const* Table_x, const T Table_f[], const T x[], const int IdxL[] )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim <=0 !!\n" );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
   if ( x == NULL )        Aux_Error( ERROR_INFO, "x == NULL !!\n" );
   if ( IdxL == NULL )     Aux_Error( ERROR_INFO, "IdxL == NULL !!\n" );
   for (int d=0; d<nDim; d++)
      if ( Table_x[d] == NULL )  Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
   for (int d=0; d<nDim; d++)
      if ( IdxL[d] >= N_x[d]-1 )  Aux_Error( ERROR_INFO, "IdxL >= N_x-1 in d=%d !!\n", d, IdxL[d], N_x[d]-1 );
#  endif


   T *xL = new T [nDim];
   T *xR = new T [nDim];
   for (int d=0; d<nDim; d++)
   {
      xL[d] = Table_x[d][IdxL[d]  ];
      xR[d] = Table_x[d][IdxL[d]+1];
      if ( x[d] < xL[d]  ||  x[d] >= xR[d] )    return NULL_REAL;
   }

   int Ncorners = (int)POW( 2, nDim );
   T *fC = new T [Ncorners];
   for (int i=0; i<Ncorners; i++)
   {
      const int Idx_f_corner = GetIdx_corner_nDim_Table( i, nDim, N_x, IdxL );
      fC[i]  = Table_f[ Idx_f_corner ];
   }

   T y = Mis_MultilinearInterpolate( nDim, xL, xR, fC, x );

   delete [] xL;
   delete [] xR;
   delete [] fC;

   return y;

} // FUNCTION : Mis_InterpolateFromNDTable_withIdx




//-------------------------------------------------------------------------------------------------------
// Function    :  GetIdxL_From1DTable
// Description :  Assuming a table of x, return the left-side index, IdxL, in the table for a given point x
//                such that Table_x[IdxL] <= x < Table_x[IdxL+1]
//
// Note        :  1. Table_x must be sorted into ascending numerical order in advance
//                2. Target point x must lie in the range Table_x[0] <= x < Table_x[N-1]
//                   --> Otherwise the function returns NULL_INT
//                3. Overloaded with different types
//
// Parameter   :  N        : Number of elements in the table Table_x
//                           --> Must be >= 2
//                Table_x  : Table of x
//                x        : Target point x for interpolation
//
// Return      :  IdxL      if x lies in the range Table_x[0] <= x < Table_x[N-1]
//                NULL_INT  if x lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
int GetIdxL_From1DTable( const int N, const T Table_x[], const T x )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( N <= 1 )           Aux_Error( ERROR_INFO, "incorrect input parameter \"N (%d) <= 1\" !!\n", N );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
#  endif


// check whether the target x lies within the accepted range
   if ( x < Table_x[0]  ||  x >= Table_x[N-1] )    return NULL_INT;


// binary search
   int IdxL;

   IdxL = Mis_BinarySearch_Real( Table_x, 0, N-1, x );

#  ifdef GAMER_DEBUG
   if ( IdxL < 0  ||  IdxL >= N-1 )
      Aux_Error( ERROR_INFO, "IdxL (%d) lies outside the accepted range [%d ... %d] !!\n", IdxL, 0, N-2 );
#  endif

   return IdxL;

} // FUNCTION : GetIdxL_From1DTable



//-------------------------------------------------------------------------------------------------------
// Function    :  GetIdx_corner_nDim_Table
// Description :  Transform the 1D left-side index, IdxL, in each dimension and
//                the target local index of target corner into the index of the corners in a nDim Table
//
// Note        :  1. Table_x must be sorted into ascending numerical order in advance
//                2. Target point x must lie in the range Table_x[0] <= x < Table_x[N-1]
//                   --> Otherwise the function returns NULL_INT
//                3. Overloaded with different types
//
// Parameter   :  Idx_corner_local : Index of the corner on a n-dimensional hypercube
//                                   range from 0 to 2^(nDim)-1,
//                                   the row-major order following the order of dimensions in N_x nad IdxL
//                                   3D Example:
//                                      0 --> x=xL, y=yL, z=zL
//                                      1 --> x=xR, y=yL, z=zL
//                                      2 --> x=xL, y=yR, z=zL
//                                      3 --> x=xR, y=yR, z=zL
//                                      4 --> x=xL, y=yL, z=zR
//                                      5 --> x=xR, y=yL, z=zR
//                                      6 --> x=xL, y=yR, z=zR
//                                      7 --> x=xR, y=yR, z=zR
//                x                : Target point x for interpolation
//                nDim             : Number of dimensions of the the table
//                N_x              : Number of elements in each dimension of the table
//                                   --> it is an array with size = nDim
//                IdxL             : Index of the leftmost corner in each dimension
//                                   --> it is an array with size = nDim
//
// Return      :  IdxL      if x lies in the range Table_x[0] <= x < Table_x[N-1]
//                NULL_INT  if x lies outside the above range
//-------------------------------------------------------------------------------------------------------
int GetIdx_corner_nDim_Table( const int Idx_corner_local, const int nDim, const int N_x[], const int IdxL[] )
{
   int *LorR = new int [nDim];

   int Idx_temp = Idx_corner_local;
   for (int d=0; d<nDim; d++)
   {
      LorR[d] = ( Idx_temp%2 == 0 ) ? 0 : 1;
      Idx_temp /= 2;
   }

   int Idx_corner_inTable = 0;
   for (int d=nDim-1; d>=0; d++)
   {
      Idx_corner_inTable *= N_x[d];
      Idx_corner_inTable += IdxL[d] + LorR[d];
   }

   delete [] LorR;

   return Idx_corner_inTable;

} // FUNCTION : GetIdx_corner_nDim_Table



// explicit template instantiation
template float  Mis_InterpolateFromTable <float>  ( const int N, const float  Table_x[], const float  Table_y[], const float  x );
template double Mis_InterpolateFromTable <double> ( const int N, const double Table_x[], const double Table_y[], const double x );
template float  Mis_InterpolateFrom_nDim_Table <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[] );
template double Mis_InterpolateFrom_nDim_Table <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[] );
template float  Mis_InterpolateFrom_nDim_Table_withIdxL <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[], const int IdxL[] );
template double Mis_InterpolateFrom_nDim_Table_withIdxL <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[], const int IdxL[] );
