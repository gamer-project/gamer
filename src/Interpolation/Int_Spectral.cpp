#include "GAMER.h"


#if ( defined(SUPPORT_FFTW) && defined(SUPPORT_GSL) )

// GSL includes
#include "GSL.h"

// FFTW includes
#include "FFTW.h"

// STL includes
#include <unordered_map>
#include <array>
#include <memory>


#define GFEI_FFTW_FLOAT8
#define GFEI_GSL_FLOAT8

// Accuracy for FFT in Gram-FE extension interpolation (GFEI)
#if ( SUPPORT_FFTW == FFTW3 )
#ifdef GFEI_FFTW_FLOAT8
namespace gfei_fftw = fftw3_double_precision;
#else // #ifdef GRAMFE_FLOAT8
namespace gfei_fftw = fftw3_single_precision;
#endif // #ifdef GRAMFE_FLOAT8 ... # else
#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace gfei_fftw = fftw2;
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... #else


// Accuracy for matrix multiplication in Gram-FE extension interpolation (GFEI)
#ifdef GFEI_GSL_FLOAT8
namespace gfei_gsl  = gsl_double_precision;
#else
namespace gfei_gsl  = gsl_single_precision;
#endif

// Accuracy for matrix multiplication in precomputed interpolation (PI)
#ifdef FLOAT8
namespace pi_gsl  = gsl_double_precision;
#else
namespace pi_gsl  = gsl_single_precision;
#endif

// Base class for interpolation context
struct InterpolationContext {
    virtual void InterpolateReal(const real* input, real *output, size_t N, size_t ghostBoundary) = 0;
};

//-------------------------------------------------------------------------------------------------------
// Structure   :  GramFEInterpolationContext
// Description :  Data structure of the Gram-Fourier extension interpolation implementation
//
// Data Member :
//             :  nInput                          : Size of input array
//             :  nGhostBoundary                  : Size of ghost boundary
//             :  nInterpolated                   : Size of interpolated input data
//             :  nExtension                      : Size of extension region
//             :  nExtended                       : Size of domain after extension
//             :  nExtendedPadded                 : Complex size of domain after extension with padding for real FFT
//             :  nDelta                          : Size of boundary domain
//             :  inputExtended                   : Array to store input array after extension
//             :  outputL                         : Array to store left interpolated values
//             :  outputR                         : Array to store right interpolated values
//             :  translationCoeffL               : Array to store left translation coefficients
//             :  translationCoeffR               : Array to store right translation coefficients
//             :  extensionMatrix                 : Array to store Gram-Fourier table
//             :  r2cPlan, c2rPlan                : Real-To-Complex and Complex-To-Real FFT plans
//
// Method      :  GramFEInterpolationContext      : Constructor
//               ~GramFEInterpolationContext      : Destructor
//                InterpolateReal                 : Interpolate real function of
//                ComputeLUDecomposition          : Compute LU decomposition with pivoting
//                ComputeChangeOfBasisMatrix      : Change from k-space to Gegenbauer space
//                ComputeInterpolationPolynomials : Evaluate Gegenbauer polynomials at interpolation points
//                ComputeInterpolationMatrix      : Interpolate values given k-space values
//                GetChangeOfBasisMatrix          : Get change-of-basis matrix
//                GetInterpolationMatrix          : Get interpolation matrix
//                GetC2CPlan                      : Get complex forward FFT plan of size nInput
//                GetR2CPlan                      : Get real-to-complex FFT plan of size nInput
//                GetPermutation                  : Get permutation of LU decomposition
//-------------------------------------------------------------------------------------------------------
class GramFEInterpolationContext  : public InterpolationContext {
//  data members
public:
    const size_t nInput, nGhostBoundary, nInterpolated, nExtension, nExtended, nExtendedPadded, nDelta;
    gfei_fftw::fft_real* inputExtended;
    gfei_fftw::fft_real* outputL;
    gfei_fftw::fft_real* outputR;
    std::array<gfei_fftw::fft_real> boundary;

private:
    std::array<gfei_fftw::std_complex> translationCoeffL, translationCoeffR;
    std::array<gfei_fftw::fft_real>    extensionMatrix;
    gfei_fftw::real_plan                r2cPlan, c2rPlan;

public:

   //-------------------------------------------------------------------------------------------------------
   // Function    :  GramFEInterpolationContext
   // Description :  Constructor of GramFEInterpolationContext, sets up FFTs and translation arrays in k-space for fast interpolation, load Gram-Fourier extension table from disk
   //
   // Parameter   :  nInput        : Size of input function to be interpolated
   //                nGhostBoundary : Size of ghost boundary
   //                nExtension    : Size of extension domain
   //                nDelta        : Size of boundary domain used for Gram-Fourier extension
   //
   //-------------------------------------------------------------------------------------------------------
    GramFEInterpolationContext(size_t nInput, size_t nGhostBoundary, size_t nExtension, size_t nDelta) :
        nInput              (nInput),
        nGhostBoundary      (nGhostBoundary),
        nInterpolated       (2 * (nInput - 2 * nGhostBoundary)),
        nExtension          (nExtension),
        nExtended           (nInput + nExtension),
        nExtendedPadded     (nExtended / 2 + 1),
        nDelta              (nDelta),
        inputExtended       (NULL),
        outputL             (NULL),
        outputR             (NULL),
        boundary            (2*nDelta),
        translationCoeffL   (nExtendedPadded),
        translationCoeffR   (nExtendedPadded),
        extensionMatrix     (nExtension * 2 * nDelta)
        {

//      allocate memory for planning FFTs
        gfei_fftw::fft_complex* planner_array = (gfei_fftw::fft_complex*) gfei_fftw::fft_malloc(nExtended * sizeof(gfei_fftw::fft_complex));

//      plan FFT of function to be interpolated
        r2cPlan  = gfei_fftw_create_1d_r2c_plan(nExtended, planner_array);
        c2rPlan  = gfei_fftw_create_1d_c2r_plan(nExtended, planner_array);

//      free memory
        gfei_fftw::fft_free(planner_array);

//      set up translation arrays
        const gfei_fftw::fft_real Norm = 1.0 / ( (gfei_fftw::fft_real) nExtended );
        gfei_fftw::fft_real  K;
        for (int i = 0; i < nExtendedPadded; i++)
        {
            K                       = 2.0 * M_PI / nExtended * i;
            translationCoeffL.at(i) = gfei_fftw::std_complex({cos(-K/4.0), sin(-K/4.0)}) * Norm;
            translationCoeffR.at(i) = gfei_fftw::std_complex({cos(+K/4.0), sin(+K/4.0)}) * Norm;
        }

//      load Gram-Fourier extension tables from files
        char filename[100];
        sprintf(filename, "gram/%ld_%ld_%d_%d.bin", nExtension, nDelta, 150, 63);
        std::array<double> doubleExtensionMatrix(nExtension * 2 * nDelta);
        ReadBinaryFile(filename, doubleExtensionMatrix.data(), 2 * nDelta * nExtension);

//      convert extension table to precision of DFT
        for (size_t i = 0; i < nExtension; ++i) {
            for (size_t j = 0; j < 2 * nDelta; ++j) {
                extensionMatrix.at(2 * nDelta * i + j) = (gfei_fftw::fft_real) doubleExtensionMatrix.at(2 * nDelta * i + j);
            }
        }

//      allocate memory for FFT arrays
        inputExtended = (gfei_fftw::fft_real*) gfei_fftw::fft_malloc(2 * sizeof(gfei_fftw::fft_real) * nExtendedPadded); // extension region
        outputL       = (gfei_fftw::fft_real*) gfei_fftw::fft_malloc(2 * sizeof(gfei_fftw::fft_real) * nExtendedPadded); // left output
        outputR       = (gfei_fftw::fft_real*) gfei_fftw::fft_malloc(2 * sizeof(gfei_fftw::fft_real) * nExtendedPadded); // right output

    } // CONSTRUCTOR : GramFEInterpolationContext


   //-------------------------------------------------------------------------------------------------------
   // Function    :  ~GramFEInterpolationContext
   // Description :  Destructor of GramFEInterpolationContext, frees memory
   //
   //-------------------------------------------------------------------------------------------------------
    ~GramFEInterpolationContext() {
        gfei_fftw::destroy_real_plan(r2cPlan);
        gfei_fftw::destroy_real_plan(c2rPlan);
        gfei_fftw::fft_free(inputExtended);
        gfei_fftw::fft_free(outputL);
        gfei_fftw::fft_free(outputR);
    } // DESTRUCTOR : GramFEInterpolationContext


   //-------------------------------------------------------------------------------------------------------
   // Function    :  InterpolateReal
   // Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
   //
   // Parameter   :  input          : Real input  array of size nInput
   //                output         : Real output array of size 2 * (nInput - nGhostBoundary)
   //                nInput         : Size of input array
   //                nGhostBoundary : Size of ghost boundary
   //
   //-------------------------------------------------------------------------------------------------------
    void InterpolateReal(const real* input, real *output)
    {
    //  convert input to FFT precision
        for (size_t i = 0; i < nInput; ++i) {
            inputExtended[i]    = (gfei_fftw::fft_real) input[ i];
        }

    //  write boundaries of input array to boundary array
        for (size_t i = 0; i < nDelta; ++i) {
            boundary[i         ] = (gfei_fftw::fft_real) input[                  i];
            boundary[i + nDelta] = (gfei_fftw::fft_real) input[nInput - nDelta + i];
        }

    //  compute extension using GSL BLAS real matrix-vector multiplication
        gfei_gsl::matrix_const_view A_view = gfei_gsl::matrix_const_view_array(extensionMatrix.data(), nExtension, 2 * nDelta);
        gfei_gsl::vector_const_view B_view = gfei_gsl::vector_const_view_array(boundary.data(),                    2 * nDelta);
        gfei_gsl::vector_view C_view       = gfei_gsl::vector_view_array      (inputExtended + nInput, nExtension);

        gfei_gsl::blas_sgemv(CblasNoTrans, 1.0, &A_view.matrix, &B_view.vector, 0.0, &C_view.vector);

    //  compute forward FFT
        gfei_fftw_r2c   (r2cPlan, inputExtended);


        gfei_fftw::std_complex* InputK   = (gfei_fftw::std_complex* ) inputExtended;
        gfei_fftw::std_complex* OutputLK = (gfei_fftw::std_complex* ) outputL;
        gfei_fftw::std_complex* OutputRK = (gfei_fftw::std_complex* ) outputR;

    //  interpolate by shifting samples to the left by 0.25 dh
        for (size_t i = 0; i < nExtendedPadded; i++)
        {
            OutputLK[i] = InputK[i] * translationCoeffL[i];
            OutputRK[i] = InputK[i] * translationCoeffR[i];
        }

    //  compute backward FFT
        gfei_fftw_c2r( c2rPlan, OutputLK );

    //  compute backward FFT
        gfei_fftw_c2r( c2rPlan, OutputRK );

    //  store output data
        for (size_t i = nGhostBoundary; i < nInput - nGhostBoundary; i++)
        {
            output[2 * (i - nGhostBoundary)    ] = outputL[i];
            output[2 * (i - nGhostBoundary) + 1] = outputR[i];
        }
    } // FUNCTION : InterpolateReal

    void ReadBinaryFile(const char* filename, double* array, int size)
    {
        FILE* file = fopen(filename, "rb");
        if (file == NULL) {
            fprintf(stderr, "Failed to open the file %s.\n", filename);
            return;
        }
        // Read the array from the binary file
        fread(array, sizeof(double), size, file);

        // Close the file
        fclose(file);
    }
}; // CLASS : GramFEInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Structure   :  PrecomputedInterpolationContext
// Description :  Data structure of the precomputed Gram-Fourier extension interpolation implementation
//
// Data Member :
//             :  nInput                          : Size of input array
//             :  nGhostBoundary                  : Size of ghost boundary
//             :  nInterpolated                   : Size of interpolated input data
//             :  interpolationMatrix             : GSL matrix of size nInterpolated x nInput
//
// Method      :  PrecomputedInterpolationContext      : Constructor
//               ~PrecomputedInterpolationContext      : Destructor
//                InterpolateReal                      : Interpolate real function of
//-------------------------------------------------------------------------------------------------------
struct PrecomputedInterpolationContext : public InterpolationClass
{
private:
    const size_t nInput, nGhostBoundary, nInterpolated;
    pi_gsl::matrix* interpolationMatrix;

public:
    PrecomputedInterpolationContext(size_t nInput, size_t nGhostBoundary) :
    nInput          (nInput),
    nGhostBoundary  (nGhostBoundary),
    nInterpolated   (2 * (nInput - 2 * nGhostBoundary))
    {
//      load interpolation table from binary file
        char filename[100];
        sprintf(filename, "gfei/N=%ld.bin", nInput);

        std::array<double> doubleInterpolationMatrix(2 * (nInput - 2) * nInput);
        ReadBinaryFile(filename, doubleInterpolationMatrix.data(), 2 * (nInput - 2) * nInput);

        interpolationMatrix = pi_gsl::matrix_alloc(nInterpolated, nInput);

        for (size_t i = 0; i < nInterpolated; ++i) {
            for (size_t j = 0; j < nInput; ++j) {
                pi_gsl::matrix_set(interpolationMatrix, i, j, (pi_gsl::gsl_real) doubleInterpolationMatrix[(i + 2 * (nGhostBoundary - 1)) * nInput + j]);
            }
        }
    } // CONSTRUCTOR : PrecomputedInterpolationContext

   //-------------------------------------------------------------------------------------------------------
   // Function    :  ~PrecomputedInterpolationContext
   // Description :  Destructor of GramFEInterpolationContext, frees memory
   //
   //-------------------------------------------------------------------------------------------------------
    ~PrecomputedInterpolationContext()
    {
        pi_gsl::matrix_free(interpolationMatrix);
    } // DESTRUCTOR : PrecomputedInterpolationContext

   //-------------------------------------------------------------------------------------------------------
   // Function    :  InterpolateReal
   // Description :  Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
   //
   // Parameter   :  input          : Real input  array of size nInput
   //                output         : Real output array of size 2 * (nInput - nGhostBoundary)
   //                nInput         : Size of input array
   //                nGhostBoundary : Size of ghost boundary
   //
   //-------------------------------------------------------------------------------------------------------
    void InterpolateReal(const real *input, real *output, size_t nInput, size_t nGhostBoundary)
    {
        pi_gsl::vector_const_view in = pi_gsl::vector_const_view_array( input, nInput       );
        pi_gsl::vector_view out      = pi_gsl::vector_view_array      (output, nInterpolated);

        pi_gsl::blas_sgemv(CblasNoTrans, 1.0, interpolationMatrix, &in.vector, 0.0, &out.vector);
    } // FUNCTION : InterpolateReal

private:

    void ReadBinaryFile(const char *filename, double *array, int size)
    {
        FILE *file = fopen(filename, "rb");
        if (file == NULL)
        {
            printf("Failed to open the file.\n");
            return;
        }

        // Read the array from the binary file
        fread(array, sizeof(double), size, file);

        // Close the file
        fclose(file);
    } // FUNCTION : ReadBinaryFile
}; // CLASS : PrecomputedInterpolationContext


//-------------------------------------------------------------------------------------------------------
// Structure   :  SpectralInterpolation
// Description :  Data structure to switch between depending spectral interpolation schemes depending on the input size
//
// Data Member :  contexts                        : Unordered map that holds interpolation contexts for different interpolations sizes
//
// Method      :  GramFEInterpolationContext      : Constructor
//               ~GramFEInterpolationContext      : Destructor
//                GetPolynomial                   : Evaluate Gegenbauer polynomial
//                ComputeLUDecomposition          : Compute LU decomposition with pivoting
//                ComputeChangeOfBasisMatrix      : Change from k-space to Gegenbauer space
//                ComputeInterpolationPolynomials : Evaluate Gegenbauer polynomials at interpolation points
//                ComputeInterpolationMatrix      : Interpolate values given k-space values
//                GetChangeOfBasisMatrix          : Get change-of-basis matrix
//                GetInterpolationMatrix          : Get interpolation matrix
//                GetC2CPlan                      : Get complex forward FFT plan of size nInput
//                GetR2CPlan                      : Get real-to-complex FFT plan of size nInput
//                GetPermutation                  : Get permutation of LU decomposition
//-------------------------------------------------------------------------------------------------------
class InterpolationHandler {
// data members
private:
    std::unordered_map<size_t, std::shared_ptr<InterpolationContext>> contexts;
    std::array<size_t> nWithSmallPrimeFactorisations = {16, 18, 24, 32, 36, 40, 48, 54, 60, 64, 72, 80, 90, 96, 100, 108, 112, 120, 128, 135, 256, 288, 512, 540};


public:
    void InterpolateReal(const real* input, real* output, size_t nInput, size_t nGhostBoundary)  {
        if (contexts.find(nInput) == contexts.end()) {

    //      for small N <= 32 pick precomputed interpolation
            if ( nInput <= 32 ) {
                contexts.insert({nInput, new PrecomputedInterpolationContext(nInput, nGhostBoundary)});
            } else {
        //      for large N >  32 use Gram-Fourier extension scheme
                size_t nExtension = 32, nDelta = 14;

                const size_t minimumExtensionSize = 18;
                const size_t maximumExtensionSize = 38;
                for (auto smallPrimeFactors : nWithSmallPrimeFactorisations) {
                    if (nInput + minimumExtensionSize < smallPrimeFactors && nInput + maximumExtensionSize >= smallPrimeFactors) {
                        nExtension = smallPrimeFactors - nInput;
                        break;
                    }
                }
                contexts.insert({nInput, new GramFEInterpolationContext(nInput, nGhostBoundary, nExtension, nDelta)})
            }

            contexts.emplace(std::piecewise_construct, std::make_tuple(nInput), std::make_tuple(nInput, nGhostBoundary, nExtension, nDelta));
        }

        contexts.at(nInput)->InterpolateReal(input, output);
    }
};


// Create InterpolationHandler object that stores all interpolationHandler contexts during runtime of GAMER
InterpolationHandler interpolationHandler;

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
         interpolationHandler.InterpolateReal(Input, Output, CRange[0] + 2 * CGhost, CGhost);

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

//       interpolate data using interpolationHandler
         interpolationHandler.InterpolateReal(Input, Output, CRange[1] + 2 * CGhost, CGhost);


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

//       interpolate data using interpolationHandler
         interpolationHandler.InterpolateReal(Input, Output, CRange[2] + 2 * CGhost, CGhost);

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