#include "GAMER.h"


#if ( defined(SUPPORT_FFTW) && defined(SUPPORT_GSL) )

// GSL includes
#include "GSL.h"

// FFTW includes
#include "FFTW.h"

// STL includes
#include <unordered_map>

//-------------------------------------------------------------------------------------------------------
// Structure   :  IPRContext
// Description :  Data structure of the IPR implementation
//
// Data Member :
//                changeOfBasisMatrix             : Complex  NxN array
//                interpolationPolynomials        : Complex 2NxN array
//                interpolationMatrix             : Complex 2NxN array
//                fftwPlan                        : FFTW plan
//                p                               : GSL permutation
//
// Method      :  IPRContext                      : Constructor
//               ~IPRContext                      : Destructor
//                GetPolynomial                   : Evaluate Gegenbauer polynomial
//                ComputeLUDecomposition          : Compute LU decomposition with pivoting
//                ComputeChangeOfBasisMatrix      : Change from k-space to Gegenbauer space
//                ComputeInterpolationPolynomials : Evaluate Gegenbauer polynomials at interpolation points
//                ComputeInterpolationMatrix      : Interpolate values given k-space values
//                GetInterpolationPolynomials     : Get interpolation values
//                GetChangeOfBasisMatrix          : Get change-of-basis matrix
//                GetInterpolationMatrix          : Get interpolation matrix
//                GetFFTWPlan                     : Get complex forward FFT plan of size N
//                GetPermutation                  : Get permutation of LU decomposition
//-------------------------------------------------------------------------------------------------------
class IPRContext
{
public:

   //-------------------------------------------------------------------------------------------------------
   // Function    :  IPRContext
   // Description :  Constructor of IPRContext, uses FFT and LU decomposition to compute a change-of-basis matrix from k-space to space of Gegenbauer polynomials
   //                Evaluates Gegenbauer polynomials of order up to N at 2 * N - 2 * ghostBoundary evaluation points and stores their values in N x ( 2 * N - 2 * ghostBoundary ) matrix
   //                Computes interpolation matrix of size N x ( 2 * N - 2 * ghostBoundary ) that takes k-space coefficients of (non-periodic) function as input and returns interpolated points
   //
   // Parameter   :  N             : Size of coarse grid function
   //                ghostBoundary : Number of coarse grid points on each side to be discarded during interpolation, must be greater than 1
   //                lambda        : Parameter of Gegenbauer polynomials
   //-------------------------------------------------------------------------------------------------------
    IPRContext(size_t N, size_t ghostBoundary, double lambda)
    {

//      complex N^2 array
        changeOfBasisMatrix      = new double[2 * N * N];

//      number of interpolation points
        size_t Nint              = 2 * (N - 2 * ghostBoundary);

//      real array with N polynomials of degree 0 - (N-1) evaluated at Nout points
        interpolationPolynomials = new double[2 * Nint * N];

//      interpolation matrix
        interpolationMatrix      = new double[2 * Nint * N];

//      allocate memory for permutation
        p                        = gsl_permutation_alloc(N);

//      plan FFT of function to be interpolated
        fftw_complex* planner_array = (fftw_complex*) fftw_malloc(N * sizeof(fftw_complex));
        fftwPlan                    = fftw_plan_dft_1d(N, planner_array, planner_array, FFTW_FORWARD, FFTW_MEASURE);
        fftw_free(planner_array);

//      compute change-of-basis matrix from polynomial to k-space
        ComputeChangeOfBasisMatrix(changeOfBasisMatrix, N, lambda);

//      compute LU decomposition of W
        int signum;
        ComputeLUDecomposition(changeOfBasisMatrix, N, p, &signum);

//      compute interpolation polynomials
        ComputeInterpolationPolynomials(interpolationPolynomials, N, lambda, Nint, ghostBoundary);

        ComputeInterpolationMatrix(interpolationMatrix, changeOfBasisMatrix, interpolationPolynomials, N, Nint);
    } // FUNCTION :  IPRContext


   //-------------------------------------------------------------------------------------------------------
   // Function    :  ~IPRContext
   // Description :  Destructor of IPRContext, frees memory
   //
   //-------------------------------------------------------------------------------------------------------
    ~IPRContext()
    {
        delete [] changeOfBasisMatrix;
        delete [] interpolationPolynomials;
        delete [] interpolationMatrix;
        gsl_permutation_free(p);
        fftw_destroy_plan(fftwPlan);
    } // FUNCTION : ~IPRContext

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetPolynomial
   // Description :  Return Gegenbauer polynomial of order n with parameter lambda evaluated at the point x
   //
   // Parameter   :  n      : Order of Polynomial
   //                lambda : Gegenbauer parameter (lambda = 0.5 -> Legendre polynomials, lambda = 1.0 -> Chebyshev polynomials)
   //                x      : Point in [-1, 1]
   //
   // Return      :  Value of Gegenbauer polynomial
   //-------------------------------------------------------------------------------------------------------
    double GetPolynomial(size_t n, double lambda, double x) const
    {
        return gsl_sf_gegenpoly_n(n, lambda, x);
    } // FUNCTION : GetPolynomial

   //-------------------------------------------------------------------------------------------------------
   // Function    :  ComputeLUDecomposition
   // Description :  Compute LU decomposition with pivoting of complex NxN matrix
   //
   // Note        :  Pointer to GSL permutation must be allocated using gsl_permutation_alloc(N) beforehand
   //
   // Parameter   :  input  : Pointer to complex double-precision array of size NxN
   //                N      : Matrix size
   //                p      : Pointer to allocated GSL permutation of size N that stores pivoting in LU decomposition
   //                signum : Pointer to integer storing the sign of the determinant of input matrix
   //
   //-------------------------------------------------------------------------------------------------------
    void ComputeLUDecomposition(double* input, size_t N, gsl_permutation* p, int* signum)
    {
        gsl_matrix_complex_view input_view = gsl_matrix_complex_view_array(input, N, N);

        gsl_linalg_complex_LU_decomp(&input_view.matrix, p, signum);
    } // FUNCTION : ComputeLUDecomposition

   //-------------------------------------------------------------------------------------------------------
   // Function    :  ComputeChangeOfBasisMatrix
   // Description :  Compute complex change-of-basis NxN matrix from k-space to Gegenbauer space with lambda
   //
   // Parameter   :  W      : Pointer to complex double-precision array of size NxN
   //                N      : Dimension of vector spaces
   //                lambda : Gegenbauer parameter (lambda = 0.5 -> Legendre polynomials, lambda = 1.0 -> Chebyshev polynomials)
   //
   //-------------------------------------------------------------------------------------------------------
    void ComputeChangeOfBasisMatrix(double* W, size_t N, double lambda) const
    {
//      allocate input array for FFT
        double*       input  = (double* )       fftw_malloc( N * 1 * sizeof( double )       );
        fftw_complex* output = (fftw_complex* ) W;

        for (size_t n = 0; n < N; ++n)
        {
//          evaluate polynomial in interval [-1, 1]
            for (size_t j = 0; j < N; ++j)
            {
                input[j] = GetPolynomial( n, lambda, -1 + j / ((double)N / 2.0) );
            }

//          create FFTW plan for double-to-complex FFT
            fftw_plan plan = fftw_plan_dft_r2c_1d(N, input, &output[n * N], FFTW_ESTIMATE);

//          compute forward FFT
            fftw_execute(plan);

//          destroy the FFTW plan
            fftw_destroy_plan(plan);

//          real-to-complex FFT maps from n to n/2 +1 because of symmetry of real FFT
//          fill up remaining n - n/2 - 1 values with complex conjugate to obtain square matrix
            for (size_t j = (N / 2 + 1); j < N; ++j)
            {
                output[ n * N + j ][0] =   output[ n * N + N - j ][0];
                output[ n * N + j ][1] = - output[ n * N + N - j ][1];
            }
        }

//      transpose output array
        gsl_matrix_complex_view output_view = gsl_matrix_complex_view_array((double* ) output, N, N);
        gsl_matrix_complex_transpose(&output_view.matrix);

//      free input array for FFT
        fftw_free(input);

    } // FUNCTION : ComputeChangeOfBasisMatrix

   //-------------------------------------------------------------------------------------------------------
   // Function    :  ComputeInterpolationPolynomials
   // Description :  Compute complex 2NxN matrix with values of Gegenbauer polynomials with parameter lambda
   //                of order [0 - N-1] evaluated at the 2 N points {-1 + 1/2N + ( j + 1 + 2 * ghostBoundary - 1)) * 1 / N for j = 0...2N}
   //                For N = 7 points and ghostBoundary = 2:
   //                i    0   1   2   3   4   5   6
   //                j           0 1 2 3 4 5
   //
   // Parameter   :  polynomialMatrix : Pointer to complex double-precision array of size 2NxN
   //                N                : Number of input points
   //                lambda           : Gegenbauer parameter (lambda = 0.5 -> Legendre polynomials, lambda = 1.0 -> Chebyshev polynomials)
   //                Nint             : Number of interpolation points
   //                ghostBoundary    : Size of ghost boundary
   //-------------------------------------------------------------------------------------------------------
    void ComputeInterpolationPolynomials(double *polynomialMatrix, size_t N, double lambda, size_t Nint, size_t ghostBoundary) const
    {
        const double dx    = 2.0 / N;
        const double x0    = -1.0 - dx / 4.0 + ghostBoundary * dx;
        double x;

//      iterate over cells in interval [-1, 1] and evaluate polynomial
        for (size_t cell = 0; cell < Nint; ++cell)
        {
//          iterate over polynomials
            for (size_t polyOrder = 0; polyOrder < N; ++polyOrder)
            {
                x = x0 + cell * dx / 2.0;
                polynomialMatrix[((cell * N) + polyOrder) * 2     ] = (real) GetPolynomial(polyOrder, lambda, x);
                polynomialMatrix[((cell * N) + polyOrder) * 2 + 1 ] = 0;
            }
        }
    } // FUNCTION : ComputeInterpolationPolynomials

    //-------------------------------------------------------------------------------------------------------
    // Function    :  ComputeInterpolationMatrix
    // Description :  Compute complex 2NxN matrix that converts the k-space value of an input function evaluated at the points {i * dx for i = 0...N}
    //                to the interpolated function values in real space at the points {i * dx/2 + dx/4 for i = 0...2N}
    //
    // Note        :  Compute inverse of upper triangular part of change-of-basis matrix and multiply with interpolation polynomials evaluated at 2N points
    // Parameter   :  transformationMatrix     : Pointer to complex NintxN output matrix
    //                changeOfBasisMatrix      : Pointer to complex    NxN change-of-basis matrix
    //                interpolationPolynomials : Pointer to complex NintxN matrix with interpolation polynomials
    //                N                        : Number of input points
    //                Nint                     : Number of interpolation points
    //-------------------------------------------------------------------------------------------------------
    void ComputeInterpolationMatrix(double* transformationMatrix, const double* changeOfBasisMatrix, const double* interpolationPolynomials, size_t N, size_t Nint)
    {
//      alllocate memory for matrix inversion and multiplication
        double* in         = new double[2 * N * N];
        double* out        = new double[2 * N * N];
//      allocate memory for LU inversion of U matrix
        gsl_permutation* p = gsl_permutation_alloc(N);
        int signum;
//      copy change of basis matrix since LU decomposition overwrites data
        memcpy(in, changeOfBasisMatrix, 2 * N * N * sizeof(double));

//      create GSL matrix views
        gsl_matrix_complex_view input_view       = gsl_matrix_complex_view_array       (in,  N, N);
        gsl_matrix_complex_view output_view      = gsl_matrix_complex_view_array       (out, N, N);
        gsl_matrix_complex_const_view poly_view  = gsl_matrix_complex_const_view_array (interpolationPolynomials, Nint, N);
        gsl_matrix_complex_view trans_view       = gsl_matrix_complex_view_array       (transformationMatrix,     Nint, N);

//      set lower triangular part of change-of-basis-matrix to zero
        for (size_t m = 1; m < N; ++m)
        {
            for (size_t n = 0; n < m; ++n)
            {
                gsl_matrix_complex_set(&input_view.matrix, m, n, {0, 0});
            }
        }

//      invert upper triangular part
        gsl_linalg_complex_LU_decomp(&input_view.matrix, p, &signum);
        gsl_linalg_complex_LU_invert(&input_view.matrix, p, &output_view.matrix);

//      multiply inverted U by polynomials matrix to obtain transformation matrix
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, {1.0, 0.0}, &poly_view.matrix, &output_view.matrix, {0.0, 0.0}, &trans_view.matrix);

//      free memory
        delete[] in;
        delete[] out;
        gsl_permutation_free(p);
    } // FUNCTION : ComputeInterpolationMatrix

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetInterpolationPolynomials
   // Description :  Return constant pointer to complex 2NxN matrix with values of Gegenbauer polynomials with parameter lambda
   //
   // Return      :  Pointer to double array of size 2 * N * N * 2
   //-------------------------------------------------------------------------------------------------------
    const double* GetInterpolationPolynomials() const
    {
        return interpolationPolynomials;
    } // FUNCTION : GetInterpolationPolynomials

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetChangeOfBasisMatrix
   // Description :  Return constant pointer to complex change-of-basis NxN matrix from k-space to Gegenbauer space with lambda
   //
   // Return      :  Pointer to double array of size 2 * N * N
   //-------------------------------------------------------------------------------------------------------
    const double* GetChangeOfBasisMatrix() const
    {
        return changeOfBasisMatrix;
    } // FUNCTION : GetChangeOfBasisMatrix

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetInterpolationMatrix
   // Description :  Return constant pointer to complex 2*NxN matrix that converts the k-space value of an input function evaluated at the points {i * dx for i = 0...N}
   //
   // Return      :  Pointer to double array of size 2 * N * N * 2
   //-------------------------------------------------------------------------------------------------------
    const double* GetInterpolationMatrix() const
    {
        return interpolationMatrix;
    } // FUNCTION : GetInterpolationMatrix

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetFFTWPlan
   // Description :  Return pointer to FFTW plan for complex, forward FFT of size N
   //
   // Return      :  Pointer to FFTW plan
   //-------------------------------------------------------------------------------------------------------
    fftw_plan GetFFTWPlan() const
    {
        return fftwPlan;
    } // FUNCTION : GetFFTWPlan

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GetPermutation
   // Description :  Return constant pointer to permutation of LU decomposition of change-of-basis matrix
   //
   // Return      :  Pointer to permutation
   //-------------------------------------------------------------------------------------------------------
    const gsl_permutation* GetPermutation() const
    {
        return p;
    } // FUNCTION : GetPermutation

private:
    double*                 changeOfBasisMatrix;
    double*                 interpolationPolynomials;
    double*                 interpolationMatrix;
    fftw_plan               fftwPlan;
    gsl_permutation*        p;
}; // CLASS : IPRContext

//-------------------------------------------------------------------------------------------------------
// Structure   :  IPR
// Description :  Interpolation using Inverse Polynomial Reconstruction (IPR)
//
// Data Member :  contexts                   : STL container for IPR contexts
//
// Method      :  InterpolateComplex         : Interpolate complex function
//                InterpolateReal            : Interpolate real function
//                GaussianElimWithTrunc      : Perform Gaussian elimination with truncation
//                Interpolate                : Interpolate FFT values using BLAS matrix multiplication
//-------------------------------------------------------------------------------------------------------
class IPR
{
public:
    //-------------------------------------------------------------------------------------------------------
    // Function    :  InterpolateComplex
    // Description :  Interpolate complex input vector given at the points {i * dx for i = 0...N-1} at the points  {i * dx/2 + dx/4 for i = 0...2N}
    //                using IPR
    // Note        :  Compute inverse of upper triangular part of change-of-basis matrix and multiply with interpolation polynomials evaluated at 2N points
    // Parameter   :  input         : Pointer to complex input vector of size N
    //                output        : Pointer to complex output vector of size 2 N
    //                N             : Size of vectors
    //                ghostBoundary : Size of ghost zone
    //-------------------------------------------------------------------------------------------------------
    void InterpolateComplex(fftw_complex *input, fftw_complex *output, size_t N, size_t ghostBoundary)
    {

        size_t Nint = 2 * (N - 2 * ghostBoundary);
        const double lambda = 20.0;

        if (contexts.find(N) == contexts.end()) {
            contexts.emplace(std::piecewise_construct, std::make_tuple(N), std::make_tuple(N, ghostBoundary, lambda));
        }

        IPRContext& c = contexts.at(N);

//      compute forward FFT
        fftw_execute_dft   (c.GetFFTWPlan(), (fftw_complex *) input, (fftw_complex *) input);

        size_t truncN = N;

        GaussianElimWithTrunc(c.GetChangeOfBasisMatrix(), (double*) input, N, truncN, 1e-9, c.GetPermutation());

        Interpolate(c.GetInterpolationMatrix(), input, output, N, truncN, Nint);
    } // FUNCTION : InterpolateComplex

    //-------------------------------------------------------------------------------------------------------
    // Function    :  InterpolateReal
    // Description :  Interpolate real input vector given at the points {i * dx for i = 0...N-1} at the points  {i * dx/2 + dx/4 for i = 0...2N}
    //                using IPR
    // Note        :  Compute inverse of upper triangular part of change-of-basis matrix and multiply with interpolation polynomials evaluated at 2N points
    // Parameter   :  input         : Pointer to real input vector of size N
    //                output        : Pointer to real output vector of size 2 N
    //                N             : Size of vectors
    //                ghostBoundary : Size of ghost zone
    //-------------------------------------------------------------------------------------------------------
    void InterpolateReal(double* re_input, double *re_output, size_t N, size_t ghostBoundary)
    {
        size_t Nint = 2 * (N - 2 * ghostBoundary);

//      compute forward FFT
        fftw_complex* input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_complex* output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nint);

        for (size_t i = 0; i < N; ++i) {
            c_re(input[i]) = re_input[i];
            c_im(input[i]) = 0.0;
        }

        InterpolateComplex(input, output, N, ghostBoundary);

        for (size_t i = 0; i < Nint; ++i) {
            re_output[i] = c_re(output[i]);
        }

        fftw_free(input);
        fftw_free(output);
    } // FUNCTION : InterpolateReal

//  interpolation coefficients
    const real R[ 1 + 2*1 ] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
    const real L[ 1 + 2*1 ] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };

    void InterpolateRealo(double* re_input, double *re_output, size_t N, size_t ghostBoundary)
    {

//      compute forward FFT
        size_t j = 0;
        for (size_t i = ghostBoundary; i < N - ghostBoundary; ++i) {
            re_output[2 * j    ] = re_input[i-1] * L[0] + re_input[i] * L[1] + re_input[i+1] * L[2];
            re_output[2 * j + 1] = re_input[i-1] * R[0] + re_input[i] * R[1] + re_input[i+1] * R[2];
            j++;
        }
    } // FUNCTION : InterpolateReal


private:

    //-------------------------------------------------------------------------------------------------------
    // Function    :  GaussianElimWithTrunc
    // Description :  Perform Gaussian elimination with truncation for stability of IPR
    //
    // Note        :  Perform forward substitution and truncate result
    // Parameter   :  LU            : Pointer to complex NxN matrix containing LU decomposition of change-of-basis matrix
    //                x             : Pointer to complex output vector of size N
    //                N             : Size of vector
    //                truncN        : Size of truncated vector
    //                truncThresh   : Threshold for truncation (should be around 10-100 * machine epsilon)
    //                p             : GSL permutation of LU decomposition
    //-------------------------------------------------------------------------------------------------------
    void GaussianElimWithTrunc(const double *LU, double *x, size_t N, size_t& truncN, double truncThresh, const gsl_permutation* p) const
    {
    //  create matrix views for input vector and matrix
        gsl_matrix_complex_const_view A_view = gsl_matrix_complex_const_view_array(LU, N, N);
        const gsl_matrix_complex          *a = &A_view.matrix;

        gsl_vector_complex_view B_view       = gsl_vector_complex_view_array(x, N);
        gsl_vector_complex          *b       = &B_view.vector;


    //  apply permutation p to b
        gsl_permute_vector_complex(p, b);

        // forward substitution to solve for Ly = B
        gsl_blas_ztrsv(CblasLower, CblasNoTrans, CblasUnit, a, b);


    //  truncation for IPR
    //  necessary for convergence for large N (roughly > 32)
    //  ( see Short Note: On the numerical convergence with the inverse polynomial reconstruction method for the resolution of the Gibbs phenomenon, Jung and Shizgal 2007)
        for (size_t m = 0; m < N; ++m)
        {
            if (gsl_complex_abs(gsl_vector_complex_get(b, m)) < truncThresh)
            {
                gsl_vector_complex_set(b, m, {0, 0});
            }
        }
    } // FUNCTION : GaussianElimWithTrunc

    //-------------------------------------------------------------------------------------------------------
    // Function    :  Interpolate
    // Description :  Interpolate by performing matrix multiplication
    //
    // Note        :  Perform forward substitution and truncate result
    // Parameter   :  interpolationMatrix : Pointer to complex NxN matrix containing LU decomposition of change-of-basis matrix
    //                g                   : Pointer to complex input vector of size N
    //                output              : Pointer to complex output vector of size N
    //                N                   : Size of vector
    //                truncN              : Size of truncated vector
    //-------------------------------------------------------------------------------------------------------
    void Interpolate(const double* interpolationMatrix, const fftw_complex* g, fftw_complex* output, size_t N, size_t truncN, size_t Nint) const
    {

        gsl_matrix_complex_const_view A_view          = gsl_matrix_complex_const_view_array((double*) interpolationMatrix, Nint, N);
        gsl_matrix_complex_const_view truncatedA_view = gsl_matrix_complex_const_submatrix(&A_view.matrix, 0, 0, Nint, truncN);

        const gsl_matrix_complex*     a      = &truncatedA_view.matrix;

        gsl_vector_complex_const_view B_view = gsl_vector_complex_const_view_array((double*) g, truncN);
        const gsl_vector_complex*     b      = &B_view.vector;

        gsl_vector_complex_view C_view       = gsl_vector_complex_view_array((double*) output, Nint);
        gsl_vector_complex*           c      = &C_view.vector;

        gsl_blas_zgemv(CblasNoTrans, {1.0, 0.0}, a, b, {0.0, 0.0}, c);

    } // FUNCTION : Interpolate

    std::unordered_map<size_t, IPRContext> contexts;
}; // CLASS : IPRContext

// Create IPR object that stores all IPR contexts during runtime of GAMER
IPR ipr;
int counter = 0;

//-------------------------------------------------------------------------------------------------------
// Function    :  Int_Spectral
// Description :  Perform spatial interpolation based on the Gram-Fourier extension method
//
// Note        :  1. The interpolation is spectrally accurate
//                2. The interpolation result is neither conservative nor monotonic
//                3. 3D interpolation is achieved by performing interpolation along x, y, and z directions
//                   in order
//
// Parameter   :  CData           : Input coarse-grid array
//                CSize           : Size of the CData array
//                CStart          : (x,y,z) starting indices to perform interpolation on the CData array
//                CRange          : Number of grids in each direction to perform interpolation
//                FData           : Output fine-grid array
//                FStart          : (x,y,z) starting indcies to store the interpolation results
//                NComp           : Number of components in the CData and FData array
//                UnwrapPhase     : Unwrap phase when OPT__INT_PHASE is on (for ELBDM only)
//                Monotonic       : Unused
//                MonoCoeff       : Unused
//                OppSign0thOrder : Unused
//-------------------------------------------------------------------------------------------------------
void Int_Spectral(  real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                    real FData[], const int FSize[3], const int FStart[3], const int NComp,
                    const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder )
{

// interpolation-scheme-dependent parameters
// ===============================================================================
// number of coarse-grid ghost zone
   const int CGhost    = 2;
// ===============================================================================

   const int maxSize   = MAX(MAX(CSize[0], CSize[1]), CSize[2]) +  2 * CGhost;

   double* Input, *Output;
   Input  = (double*) fftw_malloc( 1 * maxSize * sizeof(double) * 2);
   Output = (double*) fftw_malloc( 2 * maxSize * sizeof(double) * 2);

// index stride of the coarse-grid input array
   const int Cdx    = 1;
   const int Cdy    = Cdx*CSize[0];
   const int Cdz    = Cdy*CSize[1];

// index stride of the temporary arrays storing the data after x and y interpolations
   const int Tdx    = 1;
   const int Tdy    = Tdx* CRange[0]*2;
   const int TdzX   = Tdy*(CRange[1]+2*CGhost);    // array after x interpolation
   const int TdzY   = Tdy* CRange[1]*2;            // array after y interpolation

// index stride of the fine-grid output array
   const int Fdx    = 1;
   const int Fdy    = Fdx*FSize[0];
   const int Fdz    = Fdy*FSize[1];

// index stride of different components
   const int CDisp  = CSize[0]*CSize[1]*CSize[2];
   const int FDisp  = FSize[0]*FSize[1]*FSize[2];


   real *CPtr   = CData;
   real *FPtr   = FData;
   real *TDataX = new real [ (CRange[2]+2*CGhost)*TdzX ];   // temporary array after x interpolation
   real *TDataY = new real [ (CRange[2]+2*CGhost)*TdzY ];   // temporary array after y interpolation

   int Idx_InL, Idx_InC, Idx_InR, Idx_Out;

   for (int v=0; v<NComp; v++)
   {
    //    interpolation along x direction
      for (int In_z=CStart[2]-CGhost, Out_z=0;  In_z<CStart[2]+CRange[2]+CGhost;  In_z++, Out_z++)
      for (int In_y=CStart[1]-CGhost, Out_y=0;  In_y<CStart[1]+CRange[1]+CGhost;  In_y++, Out_y++)
      {
//       fill FFT array with one column of input data (including ghost zones) in x-direction
         for (int In_x=CStart[0]-CGhost, Out_x=0;  In_x<CStart[0]+CRange[0]+CGhost;  In_x++, Out_x++)
         {
            Idx_InC                 = In_z*Cdz  +  In_y*Cdy +  In_x*Cdx;
            Input[Out_x]            = CPtr[Idx_InC];
         } // i

//       interpolate data using IPR
         ipr.InterpolateReal(Input, Output, CRange[0] + 2 * CGhost, CGhost);

//       write result of Fourier interpolation (excluding ghost zones) to temporary array
         for (int In_x=0, Out_x=0; In_x<CRange[0]; In_x++, Out_x+=2)
         {
            Idx_Out                 = Out_z*TdzX + Out_y*Tdy + Out_x*Tdx;
            TDataX[ Idx_Out       ] = Output[2 * In_x    ];
            TDataX[ Idx_Out + Tdx ] = Output[2 * In_x + 1];
         } // i
      } // k,j

//    interpolation along y direction
      for (int InOut_z=0;             InOut_z<CRange[2]+2*CGhost;  InOut_z++)
      for (int InOut_x=0;             InOut_x<2*CRange[0];         InOut_x++)
      {

//       fill FFT array with one column of input data (including ghost boundary) in y-direction
         for (int In_y=0, Out_y=0;    In_y < CRange[1]+2*CGhost;   In_y++, Out_y++) {
            Idx_InC                 = InOut_z*TdzX +  In_y*Tdy + InOut_x*Tdx;
            Input[Out_y]            = TDataX[Idx_InC];
         } // j

//       interpolate data using IPR
         ipr.InterpolateReal(Input, Output, CRange[1] + 2 * CGhost, CGhost);


//       write result of Fourier interpolation to temporary array
         for (int In_y=0, Out_y=0;  In_y<CRange[1];   In_y++, Out_y+=2) {
            Idx_Out = InOut_z*TdzY + Out_y*Tdy + InOut_x*Tdx;
            TDataY[ Idx_Out       ] = Output[2 * In_y    ];
            TDataY[ Idx_Out + Tdy ] = Output[2 * In_y + 1];
         } // j
      } // k,i

//    interpolation along z direction
      for (int In_y=0,      Out_y=FStart[1];  In_y<2*CRange[1];       In_y++, Out_y++)
      for (int In_x=0,      Out_x=FStart[0];  In_x<2*CRange[0];       In_x++, Out_x++)
      {

//       fill FFT array with one column of input data in z-direction
         for (int In_z=0, Out_z=0;  In_z<CRange[2]+2*CGhost;  In_z++, Out_z++)
         {
            Idx_InC                 = In_z*TdzY +  In_y*Tdy +  In_x*Tdx;
            Input[Out_z]            = TDataY[Idx_InC];
         }

//       interpolate data using IPR
         ipr.InterpolateReal(Input, Output, CRange[2] + 2 * CGhost, CGhost);

         for (int In_z=0, Out_z=FStart[2];  In_z<CRange[2];  In_z++, Out_z+=2)
         {
            Idx_Out                 = Out_z*Fdz  + Out_y*Fdy + Out_x*Fdx;
            FPtr[ Idx_Out       ]   = Output[ 2 * In_z     ];
            FPtr[ Idx_Out + Fdz ]   = Output[ 2 * In_z + 1 ];
         }
      } // k,j,i

      CPtr += CDisp;
      FPtr += FDisp;
   } // for (int v=0; v<NComp; v++)

   delete [] TDataX;
   delete [] TDataY;

   fftw_free( Input  );
   fftw_free( Output );

} // FUNCTION : Int_Spectral

#endif // #if ( defined(SUPPORT_FFTW) && defined(SUPPORT_GSL) )