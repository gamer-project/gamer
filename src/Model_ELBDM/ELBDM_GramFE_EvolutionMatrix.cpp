
#include "CUFLU.h"
#include "GAMER.h"


#if ( MODEL == ELBDM && WAVE_SCHEME == WAVE_GRAMFE && GRAMFE_SCHEME == GRAMFE_MATMUL )

#include "GSL.h"

#include "GramFE_ExtensionTables.h"

#include <complex.h>


#ifdef GRAMFE_MM
namespace gramfe_mm_gsl = gsl_double_precision;
#define gramfe_mm_real double
#else // #ifdef GRAMFE_MM
namespace gramfe_mm_gsl = gsl_single_precision;
#define gramfe_mm_real float
#endif // #ifdef GRAMFE_MM ... # else


#ifdef GRAMFE_EV
namespace gramfe_ev_gsl = gsl_double_precision;
#define gramfe_ev_real double
#else // #ifdef GRAMFE_EV
namespace gramfe_ev_gsl = gsl_single_precision;
#define gramfe_ev_real float
#endif // #ifdef GRAMFE_EV ... # else


using gramfe_mm_complex_type = std::complex<gramfe_mm_real>;
using gramfe_ev_complex_type = std::complex<gramfe_ev_real>;

//-------------------------------------------------------------------------------------------------------
// Function    :  Factorial
// Description :  Compute the factorial of n.
// Parameter   :  n  : integer
// Return      :  n!
//-------------------------------------------------------------------------------------------------------
int Factorial(int n) {
     return (n==0) || (n==1) ? 1 : n* Factorial(n-1);
} // FUNCTION : Factorial

//-------------------------------------------------------------------------------------------------------
// Function    :  CosineTaylorExpansion
// Description :  Compute the taylor expansion of the cosine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
gramfe_mm_real CosineTaylorExpansion(gramfe_mm_real x, int Nterms) {
   gramfe_mm_real result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * (1 / ((gramfe_mm_real) Factorial(2 * i))) * pow(x, 2 * i   );
   }

   return result;
} // FUNCTION : CosineTaylorExpansion

//-------------------------------------------------------------------------------------------------------
// Function    :  SineTaylorExpansion
// Description :  Compute the taylor expansion of the sine function at the point x with Nterms terms
// Parameter   :  x       : Point at which to evaluate Taylor expansion
//                NTerms  : Number of terms to retain in expansion
// Return      :  Value of expansion
//-------------------------------------------------------------------------------------------------------
gramfe_mm_real SineTaylorExpansion(gramfe_mm_real x, int Nterms) {
   gramfe_mm_real result = 0;

   for (int i = 0; i  < Nterms; ++i) {
      result += pow(-1, i) * ( 1 / ((gramfe_mm_real) Factorial(2 * i + 1)) ) * pow(x, 2 * i + 1);
   }

   return result;
} // FUNCTION : SineTaylorExpansion


//-------------------------------------------------------------------------------------------------------
// Function    :  GramFE_ComputeTimeEvolutionMatrix
// Description :  Compute the time evolution matrix for the Schrödinger equation and store result in output
// Parameter   :  output  : Complex PS2 x 2 * FLU_NXT matrix (contiguous memory block of size 2 * FLU_NXT * PS2 * sizeof(real) bytes)
//                dt      : Time step
//                dh      : Grid spacing
//                Eta     : m/hbar
//-------------------------------------------------------------------------------------------------------
void ELBDM_GramFE_ComputeTimeEvolutionMatrix(real (*output)[2 * FLU_NXT], real dt, real dh, real Eta)
{
// set up time evolution operator and filter
   gramfe_mm_real K, Filter, Coeff;
   gramfe_mm_complex_type ExpCoeff;

   const gramfe_mm_real filterDecay  = (gramfe_mm_real) 32.0 * (gramfe_mm_real) 2.302585092994046; // decay of k-space filter ( 32 * log(10) )
   const gramfe_mm_real filterDegree = (gramfe_mm_real) 100;                                       // degree of k-space filter
   const gramfe_mm_real kmax         = (gramfe_mm_real) M_PI / dh;                                 // maximum value of k
   const gramfe_mm_real dk           = (gramfe_mm_real) + 2.0 * kmax / GRAMFE_FLU_NXT;             // k steps in k-space
   const gramfe_mm_real dT           = (gramfe_mm_real) - 0.5 * dt / Eta;                          // coefficient in time evolution operator

// naively "exp(1j * Coeff)" should give the exact time evolution of the free Schrödinger equation
// however, the time-evolution operator depends on higher-order derivatives
// since these are unavailable for a finite ghost boundary size, the series needs to be truncated
// the ideal Taylor expansion order using all available derivatives is FLU_GHOST_SIZE - 1
// for FLU_GHOST_SIZE == 8, 4 terms in the cosine series and 3 terms in the sine series are retained
   const int cosineNTerms = FLU_GHOST_SIZE / 2;
   const int sineNTerms   = cosineNTerms - (int) ((FLU_GHOST_SIZE % 2) == 0);

// set up momentum, filter and time evolution array
   gramfe_ev_complex_type (* Out)    [FLU_NXT]        = (gramfe_ev_complex_type (*)[       FLU_NXT]) output;
   gramfe_mm_complex_type (* FFT)    [GRAMFE_FLU_NXT] = (gramfe_mm_complex_type (*)[GRAMFE_FLU_NXT]) GramFE_FFT;

   gramfe_mm_complex_type GramFE_DFFT     [GRAMFE_FLU_NXT][GRAMFE_FLU_NXT];  //        exp(-i k^2 dt) * FFT
   gramfe_mm_complex_type GramFE_IFFTDFFT [           PS2][GRAMFE_FLU_NXT];  // IFFT * exp(-i k^2 dt) * FFT
   gramfe_mm_complex_type GramFE_Evolution[           PS2][       FLU_NXT];  // IFFT * exp(-i k^2 dt) * FFT * extension

   for (int i=0; i<GRAMFE_FLU_NXT; i++)
   {
      K           = ( i <= GRAMFE_FLU_NXT/2 ) ? dk*i : dk*(i-GRAMFE_FLU_NXT);
      Filter      = exp(-filterDecay * pow(fabs(K/kmax), 2*filterDegree));
      Coeff       = SQR(K)*dT;
      ExpCoeff    = gramfe_mm_complex_type(CosineTaylorExpansion(Coeff, cosineNTerms), SineTaylorExpansion(Coeff, sineNTerms)) * Filter;

//    multiply FFT matrix by diagonal time evolution operator in k-space
      for (int k=0; k<GRAMFE_FLU_NXT; k++)
      {
         GramFE_DFFT[i][k] = ExpCoeff * FFT[i][k];
      }
   } // for (int i=0; i<GRAMFE_FLU_NXT; i++)


   gramfe_mm_gsl::matrix_complex_const_view extend     = gramfe_mm_gsl::matrix_complex_const_view_array(GramFE_Extend                     , GRAMFE_FLU_NXT, FLU_NXT);
   gramfe_mm_gsl::matrix_complex_const_view ifft       = gramfe_mm_gsl::matrix_complex_const_view_array(GramFE_IFFT                       , PS2           , GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_const_view dfft       = gramfe_mm_gsl::matrix_complex_const_view_array((gramfe_mm_real*) GramFE_DFFT     , GRAMFE_FLU_NXT, GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_view       ifftdfft   = gramfe_mm_gsl::matrix_complex_view_array      ((gramfe_mm_real*) GramFE_IFFTDFFT , PS2,            GRAMFE_FLU_NXT);
   gramfe_mm_gsl::matrix_complex_view       evolution  = gramfe_mm_gsl::matrix_complex_view_array      ((gramfe_mm_real*) GramFE_Evolution, PS2           , FLU_NXT);

   gramfe_mm_gsl::blas_cgemm(CblasNoTrans, CblasNoTrans, {1.0,  0.0}, &ifft.matrix,     &dfft.matrix,   {0.0, 0.0}, &ifftdfft.matrix);
   gramfe_mm_gsl::blas_cgemm(CblasNoTrans, CblasNoTrans, {1.0,  0.0}, &ifftdfft.matrix, &extend.matrix, {0.0, 0.0}, &evolution.matrix);

   for (size_t i = 0; i < PS2; ++i) {
      for (size_t j = 0; j < FLU_NXT; ++j) {
         Out[i][j].real((real) GramFE_Evolution[i][j].real());
         Out[i][j].imag((real) GramFE_Evolution[i][j].imag());
      }
   }
}

#endif // #if ( MODEL == ELBDM && WAVE_SCHEME == WAVE_GRAMFE && GRAMFE_SCHEME == GRAMFE_MATMUL )