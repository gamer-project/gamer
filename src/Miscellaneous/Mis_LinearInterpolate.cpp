#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_LinearInterpolate
// Description :  Return the interpolated value of f=f(x) for a given point x in 1-dimensional grid
//                with the coordinates and values at the 2 corners provided
//
// Note        :  1. Assume xL <= x <= xR; otherwise, it is an extrapolation
//                2. Overloaded with different types
//                3. Explicit template instantiation is put at the end of this file
//
// Parameter   :  x    : Target point x coordinate at which to evaluate the interpolated value
//                xL   : Left corner x coordinate
//                xR   : Right corner x coordinate
//                f_xL : Value of f( xL )
//                f_xR : Value of f( xR )
//
// Return      :  f_x : Interpolated value f( x )
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_LinearInterpolate( const T x, const T xL, const T xR, const T f_xL, const T f_xR )
{
#  ifdef GAMER_DEBUG
   if ( xL >= xR )   Aux_Error( ERROR_INFO, "xL >= xR !!\n" );
#  endif

   const T f_x = f_xL + ( f_xR - f_xL )/( xR - xL )*( x - xL );

   return f_x;

} // FUNCTION : Mis_LinearInterpolate



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_BilinearInterpolate
// Description :  Return the interpolated value of f=f(x,y) for a given point (x,y) on a 2-dimensional rectilinear grid
//                with the coordinates and values at the 4 corners provided
//
// Note        :  1. Assume xL <= x <= xR, yL <= y <= yR; Otherwise, it is an extrapolation
//                2. Perform linear interpolation by multiplying weighting
//                3. Reference: https://en.wikipedia.org/wiki/Bilinear_interpolation
//                3. Overloaded with different types
//                4. Explicit template instantiation is put at the end of this file
//
// Parameter   :  x       : Target point x coordinate at which to evaluate the interpolated value
//                y       : Target point y coordinate at which to evaluate the interpolated value
//                xL      : Left corner x coordinate
//                xR      : Right corner x coordinate
//                yL      : Left corner y coordinate
//                yR      : Right corner y coordinate
//                f_xL_yL : Value of f( xL, yL )
//                f_xR_yL : Value of f( xR, yL )
//                f_xL_yR : Value of f( xL, yR )
//                f_xR_yR : Value of f( xR, yR )
//
// Return      :  f_x_y   : Interpolated value f( x, y )
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_BilinearInterpolate( const T x, const T y,
                           const T xL, const T xR, const T yL, const T yR,
                           const T f_xL_yL, const T f_xR_yL, const T f_xL_yR, const T f_xR_yR )
{
#  ifdef GAMER_DEBUG
   if ( xL >= xR )   Aux_Error( ERROR_INFO, "xL >= xR !!\n" );
   if ( yL >= yR )   Aux_Error( ERROR_INFO, "yL >= yR !!\n" );
#  endif

   const T weight_xR  = ( x  - xL )/( xR - xL );
   const T weight_xL  = ( xR - x  )/( xR - xL );
   const T weight_yR  = ( y  - yL )/( yR - yL );
   const T weight_yL  = ( yR - y  )/( yR - yL );

   const T f_x_y = f_xL_yL * weight_xL * weight_yL +
                   f_xR_yL * weight_xR * weight_yL +
                   f_xL_yR * weight_xL * weight_yR +
                   f_xR_yR * weight_xR * weight_yR;

/*
// Alternatively, call linear interpolation repeatedly
   const T f_xL_y = Mis_LinearInterpolate( y, yL, yR, f_xL_yL, f_xL_yR );
   const T f_xR_y = Mis_LinearInterpolate( y, yL, yR, f_xR_yL, f_xR_yR );

   const T f_x_y  = Mis_LinearInterpolate( x, xL, xR, f_xL_y,  f_xR_y  );
*/

   return f_x_y;

} // FUNCTION : Mis_BilinearInterpolate



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_TrilinearInterpolate
// Description :  Return the interpolated value of f=f(x,y,z) for a given point (x,y,z) on a 3-dimensional rectilinear grid
//                with the coordinates and values at the 8 corners provided
//
// Note        :  1. Assume xL <= x <= xR, yL <= y <= yR, zL <= z <= zR; Otherwise, it is an extrapolation
//                2. Perform linear interpolation by multiplying weighting
//                3. Reference: https://en.wikipedia.org/wiki/Trilinear_interpolation
//                4. Overloaded with different types
//                5. Explicit template instantiation is put at the end of this file
//
// Parameter   :  x          : Target point x coordinate at which to evaluate the interpolated value
//                y          : Target point y coordinate at which to evaluate the interpolated value
//                z          : Target point z coordinate at which to evaluate the interpolated value
//                xL         : Left corner x coordinate
//                xR         : Right corner x coordinate
//                yL         : Left corner y coordinate
//                yR         : Right corner y coordinate
//                zL         : Left corner z coordinate
//                zR         : Right corner z coordinate
//                f_xL_yL_zL : Value of f( xL, yL, zL )
//                f_xR_yL_zL : Value of f( xR, yL, zL )
//                f_xL_yR_zL : Value of f( xL, yR, zL )
//                f_xR_yR_zL : Value of f( xR, yR, zL )
//                f_xL_yL_zR : Value of f( xL, yL, zR )
//                f_xR_yL_zR : Value of f( xR, yL, zR )
//                f_xL_yR_zR : Value of f( xL, yR, zR )
//                f_xR_yR_zR : Value of f( xR, yR, zR )
//
// Return      :  f_x_y_z    : Interpolated value f( x, y, z )
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_TrilinearInterpolate( const T x, const T y, const T z,
                            const T xL, const T xR, const T yL, const T yR, const T zL, const T zR,
                            const T f_xL_yL_zL, const T f_xR_yL_zL, const T f_xL_yR_zL, const T f_xR_yR_zL,
                            const T f_xL_yL_zR, const T f_xR_yL_zR, const T f_xL_yR_zR, const T f_xR_yR_zR )
{
#  ifdef GAMER_DEBUG
   if ( xL >= xR )   Aux_Error( ERROR_INFO, "xL >= xR !!\n" );
   if ( yL >= yR )   Aux_Error( ERROR_INFO, "yL >= yR !!\n" );
   if ( zL >= zR )   Aux_Error( ERROR_INFO, "zL >= zR !!\n" );
#  endif

   const T weight_xR  = ( x  - xL )/( xR - xL );
   const T weight_xL  = ( xR - x  )/( xR - xL );
   const T weight_yR  = ( y  - yL )/( yR - yL );
   const T weight_yL  = ( yR - y  )/( yR - yL );
   const T weight_zR  = ( z  - zL )/( zR - zL );
   const T weight_zL  = ( zR - z  )/( zR - zL );

   const T f_z_y_x = f_xL_yL_zL * weight_xL * weight_yL * weight_zL +
                     f_xR_yL_zL * weight_xR * weight_yL * weight_zL +
                     f_xL_yR_zL * weight_xL * weight_yR * weight_zL +
                     f_xR_yR_zL * weight_xR * weight_yR * weight_zL +
                     f_xL_yL_zR * weight_xL * weight_yL * weight_zR +
                     f_xR_yL_zR * weight_xR * weight_yL * weight_zR +
                     f_xL_yR_zR * weight_xL * weight_yR * weight_zR +
                     f_xR_yR_zR * weight_xR * weight_yR * weight_zR;
/*
// Alternatively, call linear interpolation repeatedly
   const T f_xL_yL_z = Mis_LinearInterpolate( z, zL, zR, f_xL_yL_zL, f_xL_yL_zR );
   const T f_xR_yL_z = Mis_LinearInterpolate( z, zL, zR, f_xR_yL_zL, f_xR_yL_zR );
   const T f_xL_yR_z = Mis_LinearInterpolate( z, zL, zR, f_xL_yR_zL, f_xL_yR_zR );
   const T f_xR_yR_z = Mis_LinearInterpolate( z, zL, zR, f_xR_yR_zL, f_xR_yR_zR );

   const T f_z_y_x   = Mis_BilinearInterpolate( x, y,
                                                xL, xR, yL, yR,
                                                f_xL_yL_z, f_xR_yL_z, f_xL_yR_z, f_xR_yR_z );
*/

   return f_z_y_x;

} // FUNCTION : Mis_TrilinearInterpolate



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_MultilinearInterpolate
// Description :  Return the interpolated value of f=f(\vec{x}) for a given point
//                \vec{x} = (x_0, x_1, ..., x_{nDim-1}) on a n-dimensional rectilinear grid
//                with the coordinates and values at the 2^n corners provided
//
// Note        :  1. Assume xL[d] <= x[d] <= xR[d] for d from 0 to nDim-1; Otherwise, it is an extrapolation
//                2. Perform linear interpolation in each dimension recursively
//                3. Overloaded with different types
//                4. Explicit template instantiation is put at the end of this file
//
// Parameter   :  nDim     : Number of dimensions
//                x        : Target point \vec{x} at which to evaluate the interpolated value
//                           --> it is an array with size = nDim
//                xL       : Left corner coordinate in each dimension
//                           --> it is an array with size = nDim
//                xR       : Right corner coordinate in each dimension
//                           --> it is an array with size = nDim
//                fC       : Values of f( \vec{x}_corner ) at the corners
//                           --> it is an array with size = 2^(nDim)
//                           --> the fastest-changing index follows the major-axis order
//                               x_0 -> x_1 -> x_2 -> ... -> x_{nDim-1}
//
// Example     :  // For f(x,y,z) = 2 * x^3 + 3 * y^2 - z
//                double  x[3] = { 2.1, 6.2, 8.3 };
//                //                 x    y    z
//                double xL[3] = { 1.0, 4.0, 7.0 };
//                //                xL   yL   zL
//                double xR[3] = { 4.0, 7.0, 9.0 };
//                //                xR   yR   zR
//                double fC[8] = {      43.0,      169.0,      142.0,      268.0,       41.0,      167.0,       40.0,      266.0 };
//                //              f(xL,yL,zL) f(xR,yL,zL) f(xL,yR,zL) f(xR,yR,zL) f(xL,yL,zR) f(xR,yL,zR) f(xL,yR,zR) f(xR,yR,zR)
//                double f     = Mis_MultilinearInterpolate( 3, x, xL, xR, fC );
//                //           = 160.5
//
// Return      :  f_interpolated : Interpolated value f( \vec{x} )
//-------------------------------------------------------------------------------------------------------
template <typename T>
T Mis_MultilinearInterpolate( const int nDim, const T x[], const T xL[], const T xR[], const T fC[] )
{
// initial check
#  ifdef GAMER_DEBUG
   if ( nDim <= 0 )    Aux_Error( ERROR_INFO, "incorrect input parameter \"nDim (%d) <= 0\" !!\n", nDim );
   if ( xL == NULL )   Aux_Error( ERROR_INFO, "xL == NULL !!\n" );
   if ( xR == NULL )   Aux_Error( ERROR_INFO, "xR == NULL !!\n" );
   if ( fC == NULL )   Aux_Error( ERROR_INFO, "fC == NULL !!\n" );
   if ( x  == NULL )   Aux_Error( ERROR_INFO, "xT == NULL !!\n" );
   for (int d=0; d<nDim; d++)
      if ( xL[d] >= xR[d] )   Aux_Error( ERROR_INFO, "xL >= xR in d=%d !!\n", d );
#  endif

// end of the recursion
   if ( nDim == 1 )
      return Mis_LinearInterpolate( x[0], xL[0], xR[0], fC[0], fC[1] );

   if ( nDim == 2 )
      return Mis_BilinearInterpolate( x[0], x[1],
                                     xL[0], xR[0], xL[1], xR[1],
                                     fC[0], fC[1], fC[2], fC[3] );
   if ( nDim == 3 )
      return Mis_TrilinearInterpolate( x[0], x[1], x[2],
                                      xL[0], xR[0], xL[1], xR[1], xL[2], xR[2],
                                      fC[0], fC[1], fC[2], fC[3], fC[4], fC[5], fC[6], fC[7] );

// interpolation for the next lower dimension
   const int NCorner_LowerDim = (int)POW( 2, nDim-1 );
   T *fC_LowerDim = new T [NCorner_LowerDim];

   for (int i=0; i<NCorner_LowerDim; i++)
      fC_LowerDim[i] = Mis_LinearInterpolate( x[nDim-1], xL[nDim-1], xR[nDim-1], fC[i], fC[NCorner_LowerDim+i] );

   T f_interpolated = Mis_MultilinearInterpolate( nDim-1, x, xL, xR, fC_LowerDim );

   delete [] fC_LowerDim;

   return f_interpolated;

} // FUNCTION : Mis_MultilinearInterpolate



//-------------------------------------------------------------------------------------------------------
// Function    :  UnitTest_Mis_LinearInterpolate
// Description :  Test the functionality in the file
//
// Note        :
//
// Parameter   :
//
// Return      :  0 : All the tests finish
//-------------------------------------------------------------------------------------------------------
int UnitTest_Mis_LinearInterpolate()
{

// Test 1. Mis_LinearInterpolate
   double Result_1 = Mis_LinearInterpolate( 2.5, 2.0, 3.0, 4.0, -2.0 );
   double Answer_1 = 1.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_1, Answer_1, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 1 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 1 !!\n", __FUNCTION__ );


// Test 2. Mis_BilinearInterpolate
// f(x,y) = x^2 + y^2
   double Result_2 = Mis_BilinearInterpolate( 2.5, 1.5,
                                             -2.0, 4.0, 0.0, 5.0,
                                              4.0, 16.0, 29.0, 41.0 );
   double Answer_2 = 20.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_2, Answer_2, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 2 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 2 !!\n", __FUNCTION__ );


// Test 3. Mis_TrilinearInterpolate
// f(x,y,z) = 2 * x^3 + 3 * y^2 - z
   double Result_3 = Mis_TrilinearInterpolate( 2.1, 6.2, 8.3,
                                               1.0, 4.0, 4.0, 7.0, 7.0, 9.0,
                                               43.0, 169.0, 142.0, 268.0, 41.0, 167.0, 140.0, 266.0 );
   double Answer_3 = 160.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_3, Answer_3, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 3 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 3 !!\n", __FUNCTION__ );


// Test 4. Mis_MultilinearInterpolate for 1D
   float  x_4[1] = { 2.5 };
   float xL_4[1] = { 2.0 };
   float xR_4[1] = { 3.0 };
   float fC_4[2] = { 4.0, -2.0 };

   float Result_4 = Mis_MultilinearInterpolate( 1, x_4, xL_4, xR_4, fC_4 );
   float Answer_4 = 1.0;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_4, Answer_4, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 4 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 4 !!\n", __FUNCTION__ );


// Test 5. Mis_MultilinearInterpolate for 2D
// f(x,y) = x^2 + y^2
   float  x_5[2] = {  2.5, 1.5 };
   float xL_5[2] = { -2.0, 0.0 };
   float xR_5[2] = {  4.0, 5.0 };
   float fC_5[4] = {  4.0, 16.0, 29.0, 41.0 };

   float Result_5 = Mis_MultilinearInterpolate( 2, x_5, xL_5, xR_5, fC_5 );
   float Answer_5 = 20.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_5, Answer_5, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 5 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 5 !!\n", __FUNCTION__ );


// Test 6. Mis_MultilinearInterpolate for 3D
// f(x,y,z) = 2 * x^3 + 3 * y^2 - z
   double  x_6[3] = { 2.1, 6.2, 8.3 };
   double xL_6[3] = { 1.0, 4.0, 7.0 };
   double xR_6[3] = { 4.0, 7.0, 9.0 };
   double fC_6[8] = { 43.0, 169.0, 142.0, 268.0, 41.0, 167.0, 140.0, 266.0 };

   double Result_6 = Mis_MultilinearInterpolate( 3, x_6, xL_6, xR_6, fC_6 );
   double Answer_6 = 160.5;
   if ( MPI_Rank == 0 )
   if ( !Mis_CompareRealValue( Result_6, Answer_6, NULL, false ) )   Aux_Message( stdout, "Fail in %s Test 6 !!\n", __FUNCTION__ );
   else                                                              Aux_Message( stdout, "Pass in %s Test 6 !!\n", __FUNCTION__ );


   return 0;

} // FUNCTION : UnitTest_Mis_LinearInterpolate



// explicit template instantiation
template float  Mis_LinearInterpolate <float>  ( const float  x, const float  xL, const float  xR, const float  f_xL, const float  f_xR );
template double Mis_LinearInterpolate <double> ( const double x, const double xL, const double xR, const double f_xL, const double f_xR );

template float  Mis_BilinearInterpolate <float>  ( const float  x, const float  y,
                                                   const float  xL, const float  xR, const float  yL, const float  yR,
                                                   const float  f_xL_yL, const float  f_xR_yL, const float  f_xL_yR, const  float f_xR_yR );
template double Mis_BilinearInterpolate <double> ( const double x, const double y,
                                                   const double xL, const double xR, const double yL, const double yR,
                                                   const double f_xL_yL, const double f_xR_yL, const double f_xL_yR, const double f_xR_yR );

template float  Mis_TrilinearInterpolate <float>  ( const float  x, const float  y, const float  z,
                                                    const float  xL, const float  xR, const float  yL, const float  yR, const float  zL, const float  zR,
                                                    const float  f_xL_yL_zL, const float f_xR_yL_zL, const float  f_xL_yR_zL, const float  f_xR_yR_zL,
                                                    const float  f_xL_yL_zR, const float f_xR_yL_zR, const float  f_xL_yR_zR, const float  f_xR_yR_zR );
template double Mis_TrilinearInterpolate <double> ( const double x, const double y, const double z,
                                                    const double xL, const double xR, const double yL, const double yR, const double zL, const double zR,
                                                    const double f_xL_yL_zL, const double f_xR_yL_zL, const double f_xL_yR_zL, const double f_xR_yR_zL,
                                                    const double f_xL_yL_zR, const double f_xR_yL_zR, const double f_xL_yR_zR, const double f_xR_yR_zR );

template float  Mis_MultilinearInterpolate <float>  ( const int nDim, const float  x[], const float  xL[], const float  xR[], const float  fC[] );
template double Mis_MultilinearInterpolate <double> ( const int nDim, const double x[], const double xL[], const double xR[], const double fC[] );
