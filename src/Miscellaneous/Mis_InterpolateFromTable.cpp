#include "GAMER.h"

template <typename T>
static int GetIdxL_From1DCoordinateTable( const int N, const T Table_x[], const T x );
static long GetIdx_corner_nDim_Table( const int Idx_corner_local, const int nDim, const int N_x[], const int IdxL[] );



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
//                6. See the unit tests at the bottom of this file for example usage
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
   const int IdxL = GetIdxL_From1DCoordinateTable( N, Table_x, x );

   if ( IdxL == NULL_INT )    return NULL_REAL;


// linear interpolation
   T y = Mis_LinearInterpolate( x, Table_x[IdxL], Table_x[IdxL+1], Table_y[IdxL], Table_y[IdxL+1] );

   return y;

} // FUNCTION : Mis_InterpolateFromTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom2DTable
// Description :  Assuming f=f(x,y), return the interpolated value of f for a given point (x,y)
//
// Note        :  1. Interpolation table Table_x and Table_y must be sorted into ascending numerical order in advance
//                2. Target point x must lie in the range Table_x[0] <= x < Table_x[N_x-1]
//                   Target point y must lie in the range Table_y[0] <= y < Table_y[N_y-1]
//                   --> Otherwise the function returns NULL_REAL
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put at the end of this file
//                6. See the unit tests at the bottom of this file for example usage
//
// Parameter   :  N_x      : Number of elements in the interpolation table Table_x
//                           --> Must be >= 2
//                N_y      : Number of elements in the interpolation table Table_y
//                           --> Must be >= 2
//                Table_x  : Interpolation table x, size = N_x
//                Table_y  : Interpolation table y, size = N_y
//                Table_f  : Interpolation table f, size = N_y*N_x
//                x        : Target point x for interpolation
//                y        : Target point y for interpolation
//
// Return      :  f(x,y)    if x lies in the range Table_x[0] <= x < Table_x[N_x-1] and
//                             y lies in the range Table_y[0] <= y < Table_y[N_y-1]
//                NULL_REAL if (x,y) lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom2DTable( const int N_x, const int N_y,
                              const T Table_x[], const T Table_y[], const T Table_f[],
                              const T x, const T y )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_y == NULL )  Aux_Error( ERROR_INFO, "Table_y == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
#  endif


// get index
   int IdxL_x = GetIdxL_From1DCoordinateTable( N_x, Table_x, x );
   int IdxL_y = GetIdxL_From1DCoordinateTable( N_y, Table_y, y );

   if ( IdxL_x == NULL_INT  ||  IdxL_y == NULL_INT )    return NULL_REAL;


// linear interpolation
   T f = Mis_BilinearInterpolate( x, y,
                                  Table_x[IdxL_x], Table_x[IdxL_x+1],
                                  Table_y[IdxL_y], Table_y[IdxL_y+1],
                                  Table_f[ IdxL_y   *N_x +  IdxL_x   ],
                                  Table_f[ IdxL_y   *N_x + (IdxL_x+1)],
                                  Table_f[(IdxL_y+1)*N_x +  IdxL_x   ],
                                  Table_f[(IdxL_y+1)*N_x + (IdxL_x+1)] );

   return f;

} // FUNCTION : Mis_InterpolateFrom2DTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom3DTable
// Description :  Assuming f=f(x,y,z), return the interpolated value of f for a given point (x,y,z)
//
// Note        :  1. Interpolation table Table_x, Table_y, and Table_z must be sorted into ascending numerical order in advance
//                2. Target point x must lie in the range Table_x[0] <= x < Table_x[N_x-1]
//                   Target point y must lie in the range Table_y[0] <= y < Table_y[N_y-1]
//                   Target point z must lie in the range Table_z[0] <= z < Table_z[N_z-1]
//                   --> Otherwise the function returns NULL_REAL
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put at the end of this file
//                6. See the unit tests at the bottom of this file for example usage
//
// Parameter   :  N_x      : Number of elements in the interpolation table Table_x
//                           --> Must be >= 2
//                N_y      : Number of elements in the interpolation table Table_y
//                           --> Must be >= 2
//                N_z      : Number of elements in the interpolation table Table_z
//                           --> Must be >= 2
//                Table_x  : Interpolation table x, size = N_x
//                Table_y  : Interpolation table y, size = N_y
//                Table_z  : Interpolation table z, size = N_z
//                Table_f  : Interpolation table f, size = N_z*N_y*N_x
//                x        : Target point x for interpolation
//                y        : Target point y for interpolation
//                z        : Target point z for interpolation
//
// Return      :  f(x,y,z)  if x lies in the range Table_x[0] <= x < Table_x[N_x-1] and
//                             y lies in the range Table_y[0] <= y < Table_y[N_y-1] and
//                             z lies in the range Table_z[0] <= z < Table_z[N_z-1]
//                NULL_REAL if (x,y,z) lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom3DTable( const int N_x, const int N_y, const int N_z,
                              const T Table_x[], const T Table_y[], const T Table_z[], const T Table_f[],
                              const T x, const T y, const T z )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_y == NULL )  Aux_Error( ERROR_INFO, "Table_y == NULL !!\n" );
   if ( Table_z == NULL )  Aux_Error( ERROR_INFO, "Table_z == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
#  endif


// get index
   int IdxL_x = GetIdxL_From1DCoordinateTable( N_x, Table_x, x );
   int IdxL_y = GetIdxL_From1DCoordinateTable( N_y, Table_y, y );
   int IdxL_z = GetIdxL_From1DCoordinateTable( N_z, Table_z, z );

   if ( IdxL_x == NULL_INT  ||  IdxL_y == NULL_INT  ||  IdxL_z == NULL_INT )    return NULL_REAL;


// linear interpolation
   T f = Mis_TrilinearInterpolate( x, y, z,
                                   Table_x[IdxL_x], Table_x[IdxL_x+1],
                                   Table_y[IdxL_y], Table_y[IdxL_y+1],
                                   Table_z[IdxL_z], Table_z[IdxL_z+1],
                                   Table_f[IDX321( IdxL_x,   IdxL_y,   IdxL_z,   (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x+1, IdxL_y,   IdxL_z,   (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x,   IdxL_y+1, IdxL_z,   (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x+1, IdxL_y+1, IdxL_z,   (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x,   IdxL_y,   IdxL_z+1, (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x+1, IdxL_y,   IdxL_z+1, (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x,   IdxL_y+1, IdxL_z+1, (long)N_x, (long)N_y )],
                                   Table_f[IDX321( IdxL_x+1, IdxL_y+1, IdxL_z+1, (long)N_x, (long)N_y )] );

   return f;

} // FUNCTION : Mis_InterpolateFrom3DTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_InterpolateFrom_nDim_Table
// Description :  Assuming f=f(\vec{x}), return the interpolated value of f for a given point
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in an n-dimensional table
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. For interpolation, target coordinate x_d in each dimension should lie in the range Table_x[d][0] <= x_d < Table_x[d][N_x[d]-1]
//                   --> Otherwise the OutsideMethod is used
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put at the end of this file
//                6. See the unit tests at the bottom of this file for example usage
//
// Parameter   :  nDim          : Number of dimensions of Table_f
//                N_x           : Number of elements in the interpolation tables of each dimension Table_x[d]
//                                --> It is an array with size = nDim
//                                --> Must be >= 2 for each dimension
//                Table_x       : Interpolation table x of each dimension
//                                --> It is an array of arrays
//                                --> Table_x[d] has size = N_x[d], where d is the dimension from 0 to nDim-1
//                Table_f       : Interpolation table f
//                                --> It is an array with size = N_x[0]*N_x[1]*...*N_x[nDim-1]
//                x             : Target point \vec{x} for interpolation
//                                --> It is an array with size = nDim
//                OutsideMethod : Method to use when the given point is outside the range of table
//                                0 = return NULL_REAL
//                                1 = extend with the values on the boundary
//                                2 = extrapolate from the nearest points
//
// Return      :  f(\vec{x}) if \vec{x} lies in the range Table_x[d][0] <= x[d] < Table_x[d][N_x[d]-1] for d from 0 to nDim-1
//                NULL_REAL, an extended value, or an extrapolated value, depending on OutsideMethod, if \vec{x} lies outside the above range
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom_nDim_Table( const int nDim, const int N_x[], T const* const* Table_x, const T Table_f[], const T x[], const int OutsideMethod )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim (%d) <=0 !!\n", nDim );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
   if ( x == NULL )        Aux_Error( ERROR_INFO, "x == NULL !!\n" );
   for (int d=0; d<nDim; d++)
   {
      if ( N_x[d] <= 1 )          Aux_Error( ERROR_INFO, "N_x[%d] (%d) <= 1 !!\n", d, N_x[d] );
      if ( Table_x[d] == NULL )   Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
      for (int i=0; i<N_x[d]-1; i++)
      {
         if ( Table_x[d][i] >= Table_x[d][i+1] )
            Aux_Error( ERROR_INFO, "Table_x[%d][%d] (%16.8e) >= Table_x[%d][%d] (%16.8e) !!\n", d, i, Table_x[d][i], d, i+1, Table_x[d][i+1] );
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
         else                            Aux_Error( ERROR_INFO, "Unknown OutsideMethod = %d !!\n", OutsideMethod );
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
//                \vec{x} = (x_0, x_1, ..., x_{n-1}) in an n-dimensional table with the left-side indices in the table given
//
// Note        :  1. Interpolation table Table_x of each dimension must be sorted into ascending numerical order in advance
//                2. Target coordinate x_d is not necessarily in the range Table_x[d][IdxL[d]] <= x_d < Table_x[d][IdxL[d]+1]
//                   Since IdxL is given, the corresponding corners will be used anyway.
//                3. Currently the function only supports linear interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put at the end of this file
//                6. See the unit tests at the bottom of this file for example usage
//
// Parameter   :  nDim    : Number of dimensions of Table_f
//                N_x     : Number of elements in the interpolation tables of each dimension Table_x[d]
//                          --> It is an array with size = nDim
//                          --> Must be >= 2 for each dimension
//                Table_x : Interpolation table x of each dimension
//                          --> It is an array of arrays
//                          --> Table_x[d] has size = N_x[d], where d is the dimension from 0 to nDim-1
//                Table_f : Interpolation table f
//                          --> It is an array with size = N_x[0]*N_x[1]*...*N_x[nDim-1]
//                x       : Target point \vec{x} for interpolation
//                          --> It is an array with size = nDim
//                IdxL    : Index of Table_x such that Table_x[d][IdxL[d]] is the left-side corner in each dimension
//                          --> It is an array with size = nDim
//                          If IdxL[d] < 0 or IdxL[d] >= N_x[d]-1, the boundary value of the table will be extended
//
// Return      :  f_interpolated : Interpolated value of f
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_InterpolateFrom_nDim_Table_withIdxL( const int nDim, const int N_x[], T const* const* Table_x, const T Table_f[], const T x[], const int IdxL[] )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim (%d) <=0 !!\n", nDim );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( Table_x == NULL )  Aux_Error( ERROR_INFO, "Table_x == NULL !!\n" );
   if ( Table_f == NULL )  Aux_Error( ERROR_INFO, "Table_f == NULL !!\n" );
   if ( x == NULL )        Aux_Error( ERROR_INFO, "x == NULL !!\n" );
   if ( IdxL == NULL )     Aux_Error( ERROR_INFO, "IdxL == NULL !!\n" );
   for (int d=0; d<nDim; d++)
   {
      if ( N_x[d] <= 1 )          Aux_Error( ERROR_INFO, "N_x[%d] (%d) <= 1 !!\n", d, N_x[d] );
      if ( Table_x[d] == NULL )   Aux_Error( ERROR_INFO, "Table_x[%d] == NULL !!\n", d );
      for (int i=0; i<N_x[d]-1; i++)
      {
         if ( Table_x[d][i] >= Table_x[d][i+1] )
            Aux_Error( ERROR_INFO, "Table_x[%d][%d] (%16.8e) >= Table_x[%d][%d] (%16.8e) !!\n", d, i, Table_x[d][i], d, i+1, Table_x[d][i+1] );
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
   const int Ncorners = ( 1 << nDim );
   T *fC = new T [Ncorners];
   for (int i=0; i<Ncorners; i++)
   {
      const long Idx_f_corner = GetIdx_corner_nDim_Table( i, nDim, N_x, IdxL );
      fC[i] = Table_f[ Idx_f_corner ];
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
// Parameter   :  N       : Number of elements in the table Table_x
//                          --> Must be >= 2
//                Table_x : Table of x
//                x       : Target point x for interpolation
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
// Description :  Return the "1D" index of the corners in an n-dimensional table from
//                the "1D" local index of corner and the n-dimensional index of the leftmost corner
//
// Note        :  1. Invoked by Mis_InterpolateFrom_nDim_Table_withIdxL()
//
// Parameter   :  Idx_corner_local : Local index of the corner on an n-dimensional hypercube, ranging from 0 to 2^(nDim)-1.
//                                   The row-major order follows the order of dimensions in N_x and IdxL.
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
//                                   --> It is an array with size = nDim
//                IdxL             : Index of the leftmost corner in each dimension
//                                   --> It is an array with size = nDim
//                                   If IdxL[d] < 0 or IdxL[d] >= N_x[d]-1, the index of the boundary will be used
//
// Return      :  Idx_corner_inTable : The corresponding "1D" index of corner in the n-dimensional table
//-------------------------------------------------------------------------------------------------------
long GetIdx_corner_nDim_Table( const int Idx_corner_local, const int nDim, const int N_x[], const int IdxL[] )
{

// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )        Aux_Error( ERROR_INFO, "nDim (%d) <=0 !!\n", nDim );
   if ( N_x == NULL )      Aux_Error( ERROR_INFO, "N_x == NULL !!\n" );
   if ( IdxL == NULL )     Aux_Error( ERROR_INFO, "IdxL == NULL !!\n" );
   for (int d=0; d<nDim; d++)
      if ( N_x[d] <= 1 )   Aux_Error( ERROR_INFO, "N_x[%d] (%d) <= 1 !!\n", d, N_x[d] );
   if ( Idx_corner_local < 0  ||  Idx_corner_local >= ( 1L << nDim ) )
      Aux_Error( ERROR_INFO, "Idx_corner_local=%d is outside of the range !!\n", Idx_corner_local );
#  endif


// get the corresponding "1D" index in the n-dimensional table
   long Idx_corner_inTable = 0;
   for (int d=nDim-1; d>=0; d--)
   {
      Idx_corner_inTable *= N_x[d];

      if      ( IdxL[d] <  0        )   Idx_corner_inTable += 0;        // index of the left  boundary
      else if ( IdxL[d] >= N_x[d]-1 )   Idx_corner_inTable += N_x[d]-1; // index of the right boundary
      else                              Idx_corner_inTable += IdxL[d] + (bool)( Idx_corner_local & ( 1 << d ) ); // for right corner, +1 to get IdxR
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
// Return      :  0 : All the tests finish
//-------------------------------------------------------------------------------------------------------
int UnitTest_Mis_InterpolateFromTable()
{

// Example 1D table
   int    N_1D      = 10;
   double Table_X_1D[10] = { 1.0, 2.0, 3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0,  10.0 };
   double Table_f_1D[10] = { 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0 };
   int         N_x_1D[1] = { N_1D };
   double *Table_x_1D[1] = { Table_X_1D };


// Example 2D table
   int     N_X_2D       = 4;
   int     N_Y_2D       = 3;
   int     N_x_2D[2]    = { N_X_2D, N_Y_2D };

   float *Table_X_2D  = new float [N_X_2D];
   for (int i=0; i<N_X_2D; i++)
      Table_X_2D[i] = 1.0+i;   // { 1.0, 2.0, 3.0, 4.0 }

   float *Table_Y_2D  = new float [N_Y_2D];
   for (int j=0; j<N_Y_2D; j++)
      Table_Y_2D[j] = 1.0+j*j; // { 1.0, 2.0, 5.0 }

   float **Table_x_2D = new float* [2];
   Table_x_2D[0] = Table_X_2D;
   Table_x_2D[1] = Table_Y_2D;

   float *Table_f_2D = new float [N_Y_2D*N_X_2D];
   for (int j=0; j<N_Y_2D; j++)
   for (int i=0; i<N_X_2D; i++)
      Table_f_2D[ j*N_X_2D + i ] = Table_X_2D[i]*Table_Y_2D[j];
                                   // { 1.0,  2.0,  3.0,  4.0,
                                   //   2.0,  4.0,  6.0,  8.0,
                                   //   5.0, 10.0, 15.0, 20.0 }

// Example 3D table
   int     N_X_3D    = 4;
   int     N_Y_3D    = 2;
   int     N_Z_3D    = 3;
   int     N_x_3D[3] = { N_X_3D, N_Y_3D, N_Z_3D };

   double *Table_X_3D = new double [N_X_3D];
   for (int i=0; i<N_X_3D; i++)   Table_X_3D[i] = 1.0+i;   // { 1.0, 2.0, 3.0, 4.0 }

   double *Table_Y_3D = new double [N_Y_3D];
   for (int j=0; j<N_Y_3D; j++)   Table_Y_3D[j] = 2.0+j;   // { 2.0, 3.0 }

   double *Table_Z_3D = new double [N_Z_3D];
   for (int k=0; k<N_Z_3D; k++)   Table_Z_3D[k] = 4.0+k*k; // { 4.0, 5.0, 8.0 }

   double **Table_x_3D = new double* [3];
   Table_x_3D[0] = Table_X_3D;
   Table_x_3D[1] = Table_Y_3D;
   Table_x_3D[2] = Table_Z_3D;

   double *Table_f_3D = new double [N_Z_3D*N_Y_3D*N_X_3D];
   for (int k=0; k<N_Z_3D; k++)
   for (int j=0; j<N_Y_3D; j++)
   for (int i=0; i<N_X_3D; i++)
      Table_f_3D[ ( k*N_Y_3D + j )*N_X_3D + i ] = Table_X_3D[i] * SQR( Table_Y_3D[j] ) * Table_Z_3D[k];
                                                  // { 16.0,  32.0,   48.0,   64.0,
                                                  //   36.0,  72.0,  108.0,  144.0,
                                                  //
                                                  //   20.0,  40.0,   60.0,   80.0,
                                                  //   45.0,  90.0,  135.0,  180.0,
                                                  //
                                                  //   32.0,  64.0,   96.0,  128.0,
                                                  //   72.0, 144.0,  216.0,  288.0 }


// Test 1. Mis_InterpolateFromTable, normal case
   double Result_1 = Mis_InterpolateFromTable( N_1D, Table_X_1D, Table_f_1D, 5.5 );
   double Answer_1 = 30.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_1, Answer_1, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 1 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 1 !!\n", __FUNCTION__ );


// Test 2. Mis_InterpolateFromTable, x is outside the table
   double Result_2 = Mis_InterpolateFromTable( N_1D, Table_X_1D, Table_f_1D, 10.5 );
   double Answer_2 = NULL_REAL;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_2, Answer_2, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 2 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 2 !!\n", __FUNCTION__ );


// Test 3. Mis_InterpolateFrom2DTable, normal case
   float Result_3 = Mis_InterpolateFrom2DTable( N_X_2D, N_Y_2D,
                                                Table_X_2D, Table_Y_2D, Table_f_2D,
                                                (float)2.5, (float)4.0 );
   float Answer_3 = 10.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_3, Answer_3, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 3 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 3 !!\n", __FUNCTION__ );


// Test 4. Mis_InterpolateFrom3DTable, normal case
   double Result_4 = Mis_InterpolateFrom3DTable( N_X_3D, N_Y_3D, N_Z_3D,
                                                 Table_X_3D, Table_Y_3D, Table_Z_3D, Table_f_3D,
                                                 3.5, 2.2, 4.6 );
   double Answer_4 = 80.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_4, Answer_4, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 4 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 4 !!\n", __FUNCTION__ );


// Test 5. Mis_InterpolateFrom_nDim_Table, 1D normal case
   double x_5[1] = { 5.5 };
   double Result_5 = Mis_InterpolateFrom_nDim_Table( 1, N_x_1D, Table_x_1D, Table_f_1D, x_5, 0 );
   double Answer_5 = 30.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_5, Answer_5, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 5 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 5 !!\n", __FUNCTION__ );


// Test 6. Mis_InterpolateFrom_nDim_Table, 2D normal case
   float x_6[2] = { 2.5, 4.0 };
   float Result_6 = Mis_InterpolateFrom_nDim_Table( 2, N_x_2D, Table_x_2D, Table_f_2D, x_6, 0 );
   float Answer_6 = 10.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_6, Answer_6, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 6 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 6 !!\n", __FUNCTION__ );


// Test 7. Mis_InterpolateFrom_nDim_Table, 3D normal case
   double x_7[3] = { 3.5, 2.2, 4.6 };
   double Result_7 = Mis_InterpolateFrom_nDim_Table( 3, N_x_3D, Table_x_3D, Table_f_3D, x_7, 0 );
   double Answer_7 = 80.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_7, Answer_7, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 7 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 7 !!\n", __FUNCTION__ );


// Test 8. Mis_InterpolateFrom_nDim_Table, 2D and use the extended value for the point outside the table
   float x_8[2] = { 5.0, 4.0 };
   float Result_8 = Mis_InterpolateFrom_nDim_Table( 2, N_x_2D, Table_x_2D, Table_f_2D, x_8, 1 );
   float Answer_8 = 16.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_8, Answer_8, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 8 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 8 !!\n", __FUNCTION__ );


// Test 9. Mis_InterpolateFrom_nDim_Table, 2D and use the extrapolated value for the point outside the table
   float x_9[2] = { 6.0, 0.5 };
   float Result_9 = Mis_InterpolateFrom_nDim_Table( 2, N_x_2D, Table_x_2D, Table_f_2D, x_9, 2 );
   float Answer_9 = 3.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_9, Answer_9, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 9 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 9 !!\n", __FUNCTION__ );


// Test 10. Mis_InterpolateFrom_nDim_Table, 3D and use the extrapolated value for the point outside the table
   double x_10[3] = { 1.5, 4.0, 9.0 };
   double Result_10 = Mis_InterpolateFrom_nDim_Table( 3, N_x_3D, Table_x_3D, Table_f_3D, x_10, 2 );
   double Answer_10 = 189.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_10, Answer_10, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 10 !!\n", __FUNCTION__ );
   else                                                                Aux_Message( stdout, "Pass in %s Test 10 !!\n", __FUNCTION__ );


// Test 11. Mis_InterpolateFrom_nDim_Table, 2D and IdxL is provided
   float  x_11[2] = { 4.0, 5.0 };
   int IdxL_11[2] = { 0 , 0 };
   float Result_11 = Mis_InterpolateFrom_nDim_Table_withIdxL( 2, N_x_2D, Table_x_2D, Table_f_2D, x_11, IdxL_11 );
   float Answer_11 = 20.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_11, Answer_11, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 11 !!\n", __FUNCTION__ );
   else                                                                Aux_Message( stdout, "Pass in %s Test 11 !!\n", __FUNCTION__ );


// Test 12. Mis_InterpolateFrom_nDim_Table, 3D and IdxL is provided
   double x_12[3] = { 3.0, 4.0, 8.0 };
   int IdxL_12[3] = { 1, 0, 1 };
   double Result_12 = Mis_InterpolateFrom_nDim_Table_withIdxL( 3, N_x_3D, Table_x_3D, Table_f_3D, x_12, IdxL_12 );
   double Answer_12 = 336.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_12, Answer_12, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 12 !!\n", __FUNCTION__ );
   else                                                                Aux_Message( stdout, "Pass in %s Test 12 !!\n", __FUNCTION__ );


// Test 13. Mis_InterpolateFrom_nDim_Table, 3D and use the extended value when an out-of-range IdxL is provided
   double x_13[3] = { 1.5, 2.5, 8.5 };
   int IdxL_13[3] = { -1, -2, 3 };
   double Result_13 = Mis_InterpolateFrom_nDim_Table_withIdxL( 3, N_x_3D, Table_x_3D, Table_f_3D, x_13, IdxL_13 );
   double Answer_13 = 32.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_13, Answer_13, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 13 !!\n", __FUNCTION__ );
   else                                                                Aux_Message( stdout, "Pass in %s Test 13 !!\n", __FUNCTION__ );


// Test 14. GetIdx_corner_nDim_Table, 3D normal case
   int IdxL_14[3] = { 1, 1, 2 };
   int Result_14 = GetIdx_corner_nDim_Table( 1, 3, N_x_3D, IdxL_14 );
   int Answer_14 = 22;
   if ( MPI_Rank == 0 )
   if ( Result_14 != Answer_14 )   Aux_Message( stdout, "Fail in %s Test 14 !!\n", __FUNCTION__ );
   else                            Aux_Message( stdout, "Pass in %s Test 14 !!\n", __FUNCTION__ );


// free memory
   delete [] Table_X_2D;
   delete [] Table_Y_2D;
   delete [] Table_x_2D;
   delete [] Table_f_2D;

   delete [] Table_X_3D;
   delete [] Table_Y_3D;
   delete [] Table_Z_3D;
   delete [] Table_x_3D;
   delete [] Table_f_3D;


   return 0;

} // FUNCTION : UnitTest_Mis_InterpolateFromTable



// explicit template instantiation
template float  Mis_InterpolateFromTable <float>  ( const int N, const float  Table_x[], const float  Table_y[], const float  x );
template double Mis_InterpolateFromTable <double> ( const int N, const double Table_x[], const double Table_y[], const double x );
template float  Mis_InterpolateFrom2DTable <float>  ( const int N_x, const int N_y, const float  Table_x[], const float  Table_y[], const float  Table_f[], const float  x, const float  y );
template double Mis_InterpolateFrom2DTable <double> ( const int N_x, const int N_y, const double Table_x[], const double Table_y[], const double Table_f[], const double x, const double y );
template float  Mis_InterpolateFrom3DTable <float>  ( const int N_x, const int N_y, const int N_z, const float  Table_x[], const float  Table_y[], const float  Table_z[], const float  Table_f[], const float  x, const float  y, const float  z );
template double Mis_InterpolateFrom3DTable <double> ( const int N_x, const int N_y, const int N_z, const double Table_x[], const double Table_y[], const double Table_z[], const double Table_f[], const double x, const double y, const double z );
template float  Mis_InterpolateFrom_nDim_Table <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[], const int OutsideMethod );
template double Mis_InterpolateFrom_nDim_Table <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[], const int OutsideMethod );
template float  Mis_InterpolateFrom_nDim_Table_withIdxL <float>  ( const int nDim, const int N_x[], float  const* const* Table_x, const float  Table_f[], const float  x[], const int IdxL[] );
template double Mis_InterpolateFrom_nDim_Table_withIdxL <double> ( const int nDim, const int N_x[], double const* const* Table_x, const double Table_f[], const double x[], const int IdxL[] );
