#ifndef __SPECTRAL_INTERPOLATION_H__
#define __SPECTRAL_INTERPOLATION_H__


// ****************************************************************************
// ** This header declares the classes used for spectral interpolation.      **
// ** They are defined in Interpolation/Int_Spectral.cpp                     **
// ****************************************************************************

#ifdef SUPPORT_SPECTRAL_INTERPOLATION

// GSL includes
#include "GSL.h"

// FFTW includes
#include "FFTW.h"

#include <unordered_map>
#include <memory>


#ifdef SUPPORT_SPECTRAL_INTERPOLATION
#define GFEI_GSL_FLOAT8
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
#endif // #ifdef SUPPORT_SPECTRAL_INTERPOLATION

class InterpolationContext {
public:
    InterpolationContext(size_t nInput, size_t nGhostBoundary);
    virtual void   Preprocess       (real* input, const bool UnwrapPhase) const;
    virtual void   Postprocess      (const real* input, real* output, const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder) const;
    void           ReadBinaryFile   (const char* filename, double* array, int size) const;

    virtual size_t GetWorkspaceSize () const = 0;
    virtual void   InterpolateReal  (const real* input, real *output, char* workspace) const = 0;

    const size_t nInput;
    const size_t nGhostBoundary;
    const size_t nInterpolated;
};

//-------------------------------------------------------------------------------------------------------
// Structure   :  GramFEInterpolationContext
// Description :  Data structure of the Gram-Fourier extension interpolation implementation
//
// Data Member :  nInput                          : Size of input array
//                nGhostBoundary                  : Size of ghost boundary
//                nInterpolated                   : Size of interpolated input data
//                nExtension                      : Size of extension region
//                nExtended                       : Size of domain after extension
//                nExtendedPadded                 : Complex size of domain after extension with padding for real FFT
//                nDelta                          : Size of boundary domain
//                extensionMatrix                 : GSL matrix of size 2 * nDelta x nExtension to store Gram-Fourier table
//                r2cPlan, c2rPlan                : Real-To-Complex and Complex-To-Real FFT plans
//
// Method      :  GramFEInterpolationContext      : Constructor
//               ~GramFEInterpolationContext      : Destructor
//                GetWorkspaceSize                : Return size of interpolation workspace in bytes
//                InterpolateReal                 : Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
//-------------------------------------------------------------------------------------------------------
class GramFEInterpolationContext  : public InterpolationContext
{
public:
    GramFEInterpolationContext(size_t nInput, size_t nGhostBoundary, size_t nExtension, size_t nDelta);
    ~GramFEInterpolationContext();
    size_t GetWorkspaceSize () const;
    void   InterpolateReal  (const real* input, real *output, char* workspace) const;
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
// Description :  Data structure of the precomputed Gram-Fourier extension interpolation implementation
//
// Data Member :  nInput                          : Size of input array
//             :  nGhostBoundary                  : Size of ghost boundary
//             :  nInterpolated                   : Size of interpolated input data
//             :  interpolationMatrix             : GSL matrix of size nInterpolated x nInput
//
// Method      :  PrecomputedInterpolationContext : Constructor
//               ~PrecomputedInterpolationContext : Destructor
//                GetWorkspaceSize                : Return size of interpolation workspace in bytes
//                InterpolateReal                 : Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
//-------------------------------------------------------------------------------------------------------
struct PrecomputedInterpolationContext : public InterpolationContext
{
public:
    PrecomputedInterpolationContext(size_t nInput, size_t nGhostBoundary);
    ~PrecomputedInterpolationContext();
    size_t GetWorkspaceSize () const;
    void   InterpolateReal  (const real *input, real *output, char* workspace) const;
private:
    pi_gsl::matrix* interpolationMatrix;
}; // CLASS : PrecomputedInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Structure   :  QuarticInterpolationContext
// Description :  Data structure of quadratic interpolator
//
// Method      :  GetWorkspaceSize                : Return size of interpolation workspace in bytes
//                InterpolateReal                 : Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
//-------------------------------------------------------------------------------------------------------
struct QuarticInterpolationContext : public InterpolationContext
{
public:
    static const real QuarticR[ 1 + 2*2 ];
    static const real QuarticL[ 1 + 2*2 ];

    QuarticInterpolationContext(size_t nInput, size_t nGhostBoundary);
    size_t GetWorkspaceSize () const;
    void   InterpolateReal  (const real *input, real *output, char* workspace) const;
}; // CLASS : QuarticInterpolationContext

//-------------------------------------------------------------------------------------------------------
// Class       :  InterpolationHandler
// Description :  Data structure to switch between depending spectral interpolation schemes depending on the input size
//
// Data Member :  contexts                        : Unordered map that holds interpolation contexts for different interpolations sizes
//
// Method      :  GramFEInterpolationContext      : Constructor
//               ~GramFEInterpolationContext      : Destructor
//                AddInterpolationContext         : Add new interpolation context of size nInput with nGhostBoundary to contexts map
//                GetWorkspaceSize                : Return size of interpolation workspace in bytes
//                InterpolateReal                 : Interpolate input array of size nInput and store interpolation results of size 2 * (nInput - nGhostBoundary) in output array
//
//-------------------------------------------------------------------------------------------------------
class InterpolationHandler {
public:
    void   AddInterpolationContext  (size_t nInput, size_t nGhostBoundary);
    size_t GetWorkspaceSize         (size_t nInput, size_t nGhostBoundary) const;
    void   InterpolateReal          (real* input, real* output, size_t nInput, size_t nGhostBoundary,char* workspace, const bool UnwrapPhase, const bool Monotonic, const real MonoCoeff, const bool OppSign0thOrder) const;
private:
    std::unordered_map<size_t, std::shared_ptr<InterpolationContext>> contexts;
}; // CLASS : InterpolationHandler

#endif // #ifdef SUPPORT_SPECTRAL_INTERPOLATION

#endif // #ifdef __SPECTRAL_INTERPOLATION_H__