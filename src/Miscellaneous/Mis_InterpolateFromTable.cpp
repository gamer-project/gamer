#include "GAMER.h"

template <typename T>
static int GetIdxL_From1DCoordinateTable( const int N, const T Table_x[], const T x );
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
   int IdxL = GetIdxL_From1DCoordinateTable( N, Table_x, x );

   if ( IdxL == NULL_INT )    return NULL_REAL;


// linear interpolation
   T y = Mis_LinearInterpolate( x, Table_x[IdxL], Table_x[IdxL+1], Table_y[IdxL], Table_y[IdxL+1] );

   return y;

} // FUNCTION : Mis_InterpolateFromTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom_nDim_Table
// Description :  Assuming f=f(\vec{x}), return the interpolated value of f for a given point
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in n-dimensional table
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. For interpolation, target coordinate x_d in each dimension should lie in the range Table_x[d][0] <= x_d < Table_x[d][N_x[d]-1]
//                   --> Otherwise the OutsideMethod is used
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put in the end of this file
//
// Parameter   :  nDim          : Number of dimensions of the Table_f
//                N_x           : Number of elements in the interpolation tables of each dimension Table_x[d]
//                                --> it is an array with size = nDim
//                                --> Must be >= 2 for each dimension
//                Table_x       : Interpolation table x of each dimension
//                                --> it is an array of arrays
//                                --> Table_x[d] has size = N_x[d], where d is the dimension from 0 to nDim-1
//                Table_f       : Interpolation table f
//                                --> it is an array with size = N_x[0]*N_x[1]*...*N_x[nDim-1]
//                x             : Target point \vec{x} for interpolation
//                                --> it is an array with size = nDim
//                OutsideMethod : Method to use when the given point is outside the range of table
//                                0 = return NULL_REAL
//                                1 = extend with the values on the boundary
//                                2 = extrapolate from the nearest points
//
// Return      :  f(\vec{x}) if \vec{x} lies in the range Table_x[d][0] <= x[d] < Table_x[d][N_x[d]-1] for d from 0 to nDim-1
//                NULL_REAL or extended values or extrapolation  if \vec{x} lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom_nDim_Table( const int nDim, const int N_x[], T const* const* Table_x, const T Table_f[], const T x[], const int OutsideMethod )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim <=0 !!\n" );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
   if ( x == NULL )        Aux_Error( ERROR_INFO, "x == NULL !!\n" );
   for (int d=0; d<nDim; d++)
   {
      if ( N_x[d] <= 1 )          Aux_Error( ERROR_INFO, "N_x[%d] <= 1 !!\n", d );
      if ( Table_x[d] == NULL )   Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
      for (int i=0; i<N_x[d]-1; i++)
      {
         if ( Table_x[d][i] >= Table_x[d][i+1] )
            Aux_Error( ERROR_INFO, "Table_x[%d][%d] >= Table_x[%d][%d] !!\n", d, i, d, i+1 );
      }
   }
#  endif


// get index
   bool ReturnNull = false;
   int *IdxL = new int [nDim];
   for (int d=0; d<nDim; d++)
   {
      IdxL[d] = GetIdxL_From1DCoordinateTable( N_x[d], Table_x[d], x[d] );

      if ( IdxL[d] == NULL_INT )
      {
         if      ( OutsideMethod == 0 )  ReturnNull = true;
         else if ( OutsideMethod == 1 )  IdxL[d] = ( x[d] < Table_x[d][0] ) ? -1 : N_x[d]-1; // set it to invalid index for boundary extension
         else if ( OutsideMethod == 2 )  IdxL[d] = ( x[d] < Table_x[d][0] ) ?  0 : N_x[d]-2; // set it to the boundary index for two-point extrapolation
         else                            Aux_Error( ERROR_INFO, "Unknown OutsideMethod !!\n" );
      }
   }


// interpolation
   T f_interpolated;
   if ( ReturnNull )   f_interpolated = NULL_REAL;
   else                f_interpolated = Mis_InterpolateFrom_nDim_Table_withIdxL( nDim, N_x, Table_x, Table_f, x, IdxL );


// free memory
   delete [] IdxL;

   return f_interpolated;

} // FUNCTION : Mis_InterpolateFrom_nDim_Table




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom_nDim_Table_withIdxL
// Description :  Assuming f=f(\vec{x}), return the interpolated value of f for a given point
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in n-dimensional table with the left-side index in the table given
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. Target coordinate x_d is not necessarily in the range Table_x[d][IdxL[d]] <= x_d < Table_x[d][IdxL[d]+1]
//                   Since IdxL is given, the corresponding corners will be used anyway.
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
//                IdxL     : Index of Table_x which Table_x[d][IdxL[d]] is the left-side corner in each dimension
//                           --> it is an array with size = nDim
//                           if IdxL[d] < 0 or IdxL[d] >= N_x[d]-1, the boundary value of the table will be extended
//
// Return      :  f_interpolated : Interpolated value of f
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
   {
      if ( N_x[d] <= 1 )          Aux_Error( ERROR_INFO, "N_x[%d] <= 1 !!\n", d );
      if ( Table_x[d] == NULL )   Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
      for (int i=0; i<N_x[d]-1; i++)
      {
         if ( Table_x[d][i] >= Table_x[d][i+1] )
            Aux_Error( ERROR_INFO, "Table_x[%d][%d] >= Table_x[%d][%d] !!\n", d, i, d, i+1 );
      }
   }
#  endif


// prepare the coordinates of the corners
   T *xL = new T [nDim];
   T *xR = new T [nDim];
   for (int d=0; d<nDim; d++)
   {
      if ( IdxL[d] < 0 )
      {
         xL[d] = Table_x[d][0] - ( Table_x[d][1] - Table_x[d][0] );                      // arbitrary x outside the table for boundary extension
         xR[d] = Table_x[d][0];                                                          // left boundary of the table
      }
      else if ( IdxL[d] >= N_x[d]-1 )
      {
         xL[d] = Table_x[d][N_x[d]-1];                                                   // right boundary of the table
         xR[d] = Table_x[d][N_x[d]-1] + ( Table_x[d][N_x[d]-1] - Table_x[d][N_x[d]-2] ); // arbitrary x outside the table for boundary extension
      }
      else
      {
         xL[d] = Table_x[d][IdxL[d]  ];
         xR[d] = Table_x[d][IdxL[d]+1];
      }
   }


// prepare the values at the corners
   int Ncorners = (int)POW( 2, nDim );
   T *fC = new T [Ncorners];
   for (int i=0; i<Ncorners; i++)
   {
      const int Idx_f_corner = GetIdx_corner_nDim_Table( i, nDim, N_x, IdxL );
      fC[i]  = Table_f[ Idx_f_corner ];
   }


// linear interpolation
   T f_interpolated = Mis_MultilinearInterpolate( nDim, x, xL, xR, fC );


// free memory
   delete [] xL;
   delete [] xR;
   delete [] fC;

   return f_interpolated;

} // FUNCTION : Mis_InterpolateFrom_nDim_Table_withIdxL



//-------------------------------------------------------------------------------------------------------
// Function    :  GetIdxL_From1DCoordinateTable
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
int GetIdxL_From1DCoordinateTable( const int N, const T Table_x[], const T x )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( N <= 1 )           Aux_Error( ERROR_INFO, "incorrect input parameter \"N (%d) <= 1\" !!\n", N );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
#  endif


// check whether the target x lies within the accepted range
   if ( x < Table_x[0]  ||  x >= Table_x[N-1] )    return NULL_INT;


// binary search
   int IdxL = Mis_BinarySearch_Real( Table_x, 0, N-1, x );

#  ifdef GAMER_DEBUG
   if ( IdxL < 0  ||  IdxL >= N-1 )
      Aux_Error( ERROR_INFO, "IdxL (%d) lies outside the accepted range [%d ... %d] !!\n", IdxL, 0, N-2 );
#  endif

   return IdxL;

} // FUNCTION : GetIdxL_From1DCoordinateTable



//-------------------------------------------------------------------------------------------------------
// Function    :  GetIdx_corner_nDim_Table
// Description :  Return the "1D" index of the corners in a n-dimensional table from
//                the "1D" local index of corner and the n-dimensional index of the leftmost corner
//
// Note        :  1.
//
// Parameter   :  Idx_corner_local : Index of the corner on a n-dimensional hypercube
//                                   range from 0 to 2^(nDim)-1,
//                                   the row-major order following the order of dimensions in N_x and IdxL
//                                   3D Example:
//                                      0 --> x=xL, y=yL, z=zL
//                                      1 --> x=xR, y=yL, z=zL
//                                      2 --> x=xL, y=yR, z=zL
//                                      3 --> x=xR, y=yR, z=zL
//                                      4 --> x=xL, y=yL, z=zR
//                                      5 --> x=xR, y=yL, z=zR
//                                      6 --> x=xL, y=yR, z=zR
//                                      7 --> x=xR, y=yR, z=zR
//                nDim             : Number of dimensions of the table
//                N_x              : Number of elements in each dimension of the table
//                                   --> it is an array with size = nDim
//                IdxL             : Index of the leftmost corner in each dimension
//                                   --> it is an array with size = nDim
//                                   if IdxL[d] < 0 or IdxL[d] >= N_x[d]-1, the index of the boundary will be used
//
// Return      :  Idx_corner_inTable : The corresponding "1D" index of corner in the n-dimensional table
//-------------------------------------------------------------------------------------------------------
int GetIdx_corner_nDim_Table( const int Idx_corner_local, const int nDim, const int N_x[], const int IdxL[] )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim <=0 !!\n" );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( IdxL == NULL )     Aux_Error( ERROR_INFO, "IdxL == NULL !!\n" );
   for (int d=0; d<nDim; d++)
      if ( N_x[d] <= 1 )   Aux_Error( ERROR_INFO, "N_x[%d] <= 1 !!\n", d );
   if ( Idx_corner_local < 0  ||  Idx_corner_local >= POW( 2, nDim ) )
      Aux_Error( ERROR_INFO, "Idx_corner_local=%d is outside of the range !!\n", Idx_corner_local );
#  endif


// get the corresponding "1D" index in the n-dimensional table
   int Idx_corner_inTable = 0;
   for (int d=nDim-1; d>=0; d--)
   {
      Idx_corner_inTable *= N_x[d];

      if      ( IdxL[d] <  0        )   Idx_corner_inTable += 0;        // index of the left  boundary, no matter what LorR
      else if ( IdxL[d] >= N_x[d]-1 )   Idx_corner_inTable += N_x[d]-1; // index of the right boundary, no matter what LorR
      else                              Idx_corner_inTable += IdxL[d] + (bool)( Idx_corner_local & (int)POW( 2, d ) ); // if right corner, +1 to get IdxR
   }

   return Idx_corner_inTable;

} // FUNCTION : GetIdx_corner_nDim_Table



//-------------------------------------------------------------------------------------------------------
// Function    :  UnitTest_Mis_InterpolateFromTable
// Description :  Test the functionality in the file
//
// Note        :
//
// Parameter   :
//
// Return      :  0     : All the tests finish
//-------------------------------------------------------------------------------------------------------
int UnitTest_Mis_InterpolateFromTable()
{

// Test 1. Mis_LinearInterpolate
   double Result_1 = Mis_LinearInterpolate( 2.5, 2.0, 3.0, 4.0, -2.0 );
   double Answer_1 = 1.0;
   if ( !Mis_CompareRealValue( Result_1, Answer_1, NULL, false ) )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Fail in %s Test 1 !!\n", __FUNCTION__ );
   }
   else
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Pass in %s Test 1 !!\n", __FUNCTION__ );
   }


// Test 2. Mis_BilinearInterpolate
// f(x,y) = x^2 + y^2
   double Result_2 = Mis_BilinearInterpolate( 2.5, 1.5,
                                             -2.0, 4.0, 0.0, 5.0,
                                              4.0, 16.0, 29.0, 41.0 );
   double Answer_2 = 20.5;
   if ( !Mis_CompareRealValue( Result_2, Answer_2, NULL, false ) )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Fail in %s Test 2 !!\n", __FUNCTION__ );
   }
   else
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "Pass in %s Test 2 !!\n", __FUNCTION__ );
   }


   return 0;

} // FUNCTION : UnitTest_Mis_InterpolateFromTable



// explicit template instantiation
template float  Mis_InterpolateFromTable <float>  ( const int N, const float  Table_x[], const float  Table_y[], const float  x );
template double Mis_InterpolateFromTable <double> ( const int N, const double Table_x[], const double Table_y[], const double x );
template float  Mis_InterpolateFrom_nDim_Table <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[] );
template double Mis_InterpolateFrom_nDim_Table <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[] );
template float  Mis_InterpolateFrom_nDim_Table_withIdxL <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[], const int IdxL[] );
template double Mis_InterpolateFrom_nDim_Table_withIdxL <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[], const int IdxL[] );
