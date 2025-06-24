#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "../include/General.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Lambda
// Description :
// Note        :
// Parameter   :
// Return      :  lambdaAtTemp
//-------------------------------------------------------------------------------------------------------
real Lambda( real Temp, int numRow, real *tempTable, real *lambdaTable )
{
   real lambdaAtTemp;
   const int Idx = BinarySearch( tempTable, 0, numRow-1, Temp );

// The target value is not between the sorted array
   if ( Idx < 0  ||  Idx >= numRow-2 )   lambdaAtTemp = 0.0;
   else                                  lambdaAtTemp = LinearInterpolation( tempTable[Idx  ], lambdaTable[Idx  ],
                                                                             tempTable[Idx+1], lambdaTable[Idx+1], Temp );

//  gsl_interp_accel *acc    = gsl_interp_accel_alloc ();
//  gsl_spline       *spline = gsl_spline_alloc (gsl_interp_cspline, numRow);
//
//  gsl_spline_init (spline, tempTable, lambdaTable, numRow);
//
//  lambdaAtTemp = gsl_spline_eval (spline, Temp, acc);
//
//  gsl_spline_free (spline);
//  gsl_interp_accel_free (acc);

   return lambdaAtTemp;
} // FUNCTION : Lambda
