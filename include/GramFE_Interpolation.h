#ifndef __SPECTRAL_INT_H__
#define __SPECTRAL_INT_H__

#ifdef SUPPORT_SPECTRAL_INT



// ****************************************************************************
// ** This header declares the classes used for spectral interpolation.      **
// ** They are defined in Interpolation/Int_Spectral.cpp                     **
// ****************************************************************************

#include "GSL.h"
#include "FFTW.h"
#include <unordered_map>
#include <memory>


#define GFEI_GSL_FLOAT8
// accuracy for matrix multiplication in Gram-FE extension interpolation (GFEI)
#ifdef GFEI_GSL_FLOAT8
namespace gfei_gsl = gsl_double_precision;
#else
namespace gfei_gsl = gsl_single_precision;
#endif

// accuracy for matrix multiplication in precomputed interpolation (PI)
#ifdef FLOAT8
namespace pi_gsl   = gsl_double_precision;
#else
namespace pi_gsl   = gsl_single_precision;
#endif



//-------------------------------------------------------------------------------------------------------
// Class       :  InterpolationContext
// Description :  Base class for interpolation contexts
//
// Data Member :  nInput               : Size of input array
//                nGhostBoundary       : Size of ghost boundary
//                nInterpolated        : Size of interpolated input data
//
// Method      :  InterpolationContext : Constructor
//               ~InterpolationContext : Destructor
//                Preprocess           : Prepare input array by unwrapping phase etc.
//                Postprocess          : Postprocess interpolation results by enforcing monotonicity etc.
//                ReadBinaryFile       : Read "size" doubles from binary file "filename" into double array "array"
//                GetWorkspaceSize     : Purely virtual; Return size of interpolation workspace in bytes
//                InterpolateReal      : Purely virtual; Interpolate input array of size nInput and store interpolation results
//                                       of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
class InterpolationContext {
public:
    InterpolationContext( size_t nInput, size_t nGhostBoundary );
    virtual void   Preprocess( real* input, const bool UnwrapPhase ) const;
    virtual void   Postprocess( const real* input, real* output, const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder ) const;
    void           ReadBinaryFile( const char* filename, double* array, const int size ) const;

    virtual size_t GetWorkspaceSize() const = 0;
    virtual void   InterpolateReal( const real* input, real *output, char* workspace ) const = 0;

    const size_t   nInput;
    const size_t   nGhostBoundary;
    const size_t   nInterpolated;
}; // CLASS : InterpolationContext



//-------------------------------------------------------------------------------------------------------
// Class       :  GramFEInterpolationContext
// Description :  Data structure of the Gram-Fourier extension interpolation implementation derived from InterpolationContext
//
// Data Member :  nInput                     : Size of input array
//                nGhostBoundary             : Size of ghost boundary
//                nInterpolated              : Size of interpolated input data
//                nExtension                 : Size of extension region
//                nExtended                  : Size of domain after extension
//                nExtendedPadded            : Complex size of domain after extension with padding for real FFT
//                nDelta                     : Size of boundary domain
//                r2cPlan                    : Real-To-Complex FFT plan
//                c2rPlan                    : Complex-To-Real FFT plan
//                translationCoeffL          : Array of nExtendedPadded complex coefficients for translation to left by dx/4 in k-space
//                translationCoeffR          : Array of nExtendedPadded complex coefficients for translation to right by dx/4 in k-space
//                extensionMatrix            : GSL matrix of size 2*nDelta*nExtension to store Gram-Fourier table
//
// Method      :  GramFEInterpolationContext : Constructor
//               ~GramFEInterpolationContext : Destructor
//                GetWorkspaceSize           : Return size of interpolation workspace in bytes
//                InterpolateReal            : Interpolate input array of size nInput and store interpolation results
//                                             of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
class GramFEInterpolationContext : public InterpolationContext
{
public:
    GramFEInterpolationContext( size_t nInput, size_t nGhostBoundary, size_t nExtension, size_t nDelta );
    ~GramFEInterpolationContext();
    size_t GetWorkspaceSize() const;
    void   InterpolateReal( const real* input, real *output, char* workspace ) const;

private:
    const size_t nExtension;
    const size_t nExtended;
    const size_t nExtendedPadded;
    const size_t nDelta;
    gfei_fftw::real_plan_1d r2cPlan;
    gfei_fftw::real_plan_1d c2rPlan;
    gfei_fftw::fft_complex* translationCoeffL;
    gfei_fftw::fft_complex* translationCoeffR;
    gfei_gsl::matrix*       extensionMatrix;
}; // CLASS : GramFEInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Structure   :  PrecomputedInterpolationContext
// Description :  Data structure of the precomputed Gram-Fourier extension interpolation implementation derived from InterpolationContext
//
// Data Member :  interpolationMatrix             : GSL matrix of size nInterpolated x nInput
//
// Method      :  PrecomputedInterpolationContext : Constructor
//               ~PrecomputedInterpolationContext : Destructor
//                GetWorkspaceSize                : Return size of interpolation workspace in bytes
//                InterpolateReal                 : Interpolate input array of size nInput and store interpolation results
//                                                  of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
struct PrecomputedInterpolationContext : public InterpolationContext
{
public:
    PrecomputedInterpolationContext( size_t nInput, size_t nGhostBoundary );
    ~PrecomputedInterpolationContext();
    size_t GetWorkspaceSize() const;
    void   InterpolateReal( const real *input, real *output, char* workspace ) const;

private:
    pi_gsl::matrix* interpolationMatrix;
}; // CLASS : PrecomputedInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Structure   :  QuarticInterpolationContext
// Description :  Data structure of quartic interpolator derived from InterpolationContext
//
// Data Member :  QuarticR         : Quartic right interpolation coefficients
//                QuarticL         : Quartic left interpolation coefficients
//
// Method      :  GetWorkspaceSize : Return size of interpolation workspace in bytes
//                InterpolateReal  : Interpolate input array of size nInput and store interpolation results
//                                   of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
struct QuarticInterpolationContext : public InterpolationContext
{
public:
    static const real QuarticR[5];
    static const real QuarticL[5];

    QuarticInterpolationContext( size_t nInput, size_t nGhostBoundary );
    size_t GetWorkspaceSize() const;
    void   InterpolateReal( const real *input, real *output, char* workspace ) const;
}; // CLASS : QuarticInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Structure   :  CQuarticInterpolationContext
// Description :  Data structure of quartic interpolator derived from InterpolationContext
//
// Data Member :  QuarticR         : Quartic right interpolation coefficients
//                QuarticL         : Quartic left interpolation coefficients
//
// Method      :  GetWorkspaceSize : Return size of interpolation workspace in bytes
//                InterpolateReal  : Interpolate input array of size nInput and store interpolation results
//                                   of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
struct CQuarticInterpolationContext : public InterpolationContext
{
public:
    static const real CQuartic[5];

    CQuarticInterpolationContext( size_t nInput, size_t nGhostBoundary );
    size_t GetWorkspaceSize() const;
    void   InterpolateReal( const real *input, real *output, char* workspace ) const;
}; // CLASS : CQuarticInterpolationContext



//-------------------------------------------------------------------------------------------------------
// Class       :  InterpolationHandler
// Description :  Data structure to switch between spectral interpolation schemes depending on the input size
//
// Data Member :  contexts                   : Unordered map that holds interpolation contexts for different interpolations sizes
//
// Method      :  GramFEInterpolationContext : Constructor
//               ~GramFEInterpolationContext : Destructor
//                AddInterpolationContext    : Add new interpolation context of size nInput with nGhostBoundary to contexts map
//                GetWorkspaceSize           : Return size of interpolation workspace in bytes
//                InterpolateReal            : Interpolate input array of size nInput and store interpolation results
//                                             of size 2*(nInput - nGhostBoundary) in output array
//-------------------------------------------------------------------------------------------------------
class InterpolationHandler {
public:
    void   AddInterpolationContext( size_t nInput, size_t nGhostBoundary );
    size_t GetWorkspaceSize( size_t nInput, size_t nGhostBoundary ) const;
    void   InterpolateReal( real* input, real* output, size_t nInput, size_t nGhostBoundary, char* workspace,
                            const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder ) const;

private:
    std::unordered_map< size_t, std::shared_ptr<InterpolationContext> > contexts;
}; // CLASS : InterpolationHandler



#endif // #ifdef SUPPORT_SPECTRAL_INT

#endif // #ifdef __SPECTRAL_INT_H__
