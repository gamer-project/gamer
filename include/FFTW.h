#ifndef __FFTW_H__
#define __FFTW_H__



// ****************************************************************************
// ** This header defines wrappers for functions in the FFTW 2 and 3         **
// ** libraries in single and double precision.                              **
// ****************************************************************************


#include "GAMER.h"


#if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )

// forward declaration of std::complex
namespace std {
    template<typename T> class complex;
}

// wrappers for fftw3 single and double precision routines
#if ( SUPPORT_FFTW == FFTW3 )

namespace fftw3_single_precision {
using      fft_real                     = float;
using      fft_complex                  = fftwf_complex;
using      std_complex                  = std::complex<float>;
using      plan                         = fftwf_plan;
using      real_plan                    = fftwf_plan;
using      real_nd_plan                 = fftwf_plan;
using      complex_plan                 = fftwf_plan;
using      complex_nd_plan              = fftwf_plan;
const auto fft_malloc                   = fftwf_malloc;
const auto fft_free                     = fftwf_free;
const auto destroy_plan                 = fftwf_destroy_plan;
const auto destroy_real_plan            = fftwf_destroy_plan;
const auto destroy_real_nd_plan         = fftwf_destroy_plan;
const auto destroy_complex_plan         = fftwf_destroy_plan;
const auto destroy_complex_nd_plan      = fftwf_destroy_plan;
const auto execute_dft_r2c              = fftwf_execute_dft_r2c;
const auto execute_dft_c2r              = fftwf_execute_dft_c2r;
const auto execute_dft_c2c              = fftwf_execute_dft;
const auto execute_nd_dft_r2c           = fftwf_execute_dft_r2c;
const auto execute_nd_dft_c2r           = fftwf_execute_dft_c2r;
const auto execute_nd_dft_c2c           = fftwf_execute_dft;
const auto plan_dft_r2c_3d              = fftwf_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d              = fftwf_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d              = fftwf_plan_dft_3d;
const auto plan_dft_c2c_1d              = fftwf_plan_dft_1d;
const auto plan_dft_c2r_1d              = fftwf_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d              = fftwf_plan_dft_r2c_1d;
const auto cleanup                      = fftwf_cleanup;
#ifndef SERIAL
using      real_mpi_nd_plan             = fftwf_plan;
using      complex_mpi_nd_plan          = fftwf_plan;
const auto destroy_real_mpi_nd_plan     = fftwf_destroy_plan;
const auto destroy_complex_mpi_nd_plan  = fftwf_destroy_plan;
const auto mpi_execute_nd_dft_r2c       = fftwf_mpi_execute_dft_r2c;
const auto mpi_execute_nd_dft_c2r       = fftwf_mpi_execute_dft_c2r;
const auto mpi_execute_nd_dft_c2c       = fftwf_mpi_execute_dft;
const auto mpi_plan_dft_r2c_3d          = fftwf_mpi_plan_dft_r2c_3d;
const auto mpi_plan_dft_c2r_3d          = fftwf_mpi_plan_dft_c2r_3d;
const auto mpi_plan_dft_c2c_3d          = fftwf_mpi_plan_dft_3d;
const auto mpi_cleanup                  = fftwf_mpi_cleanup;
#endif // #ifndef SERIAL
};

namespace fftw3_double_precision {
using      fft_real                     = double;
using      fft_complex                  = fftw_complex;
using      std_complex                  = std::complex<double>;
using      real_plan                    = fftw_plan;
using      real_nd_plan                 = fftw_plan;
using      complex_plan                 = fftw_plan;
using      complex_nd_plan              = fftw_plan;
const auto fft_malloc                   = fftw_malloc;
const auto fft_free                     = fftw_free;
const auto destroy_plan                 = fftw_destroy_plan;
const auto destroy_real_plan            = fftw_destroy_plan;
const auto destroy_real_nd_plan         = fftw_destroy_plan;
const auto destroy_complex_plan         = fftw_destroy_plan;
const auto destroy_complex_nd_plan      = fftw_destroy_plan;
const auto execute_dft_r2c              = fftw_execute_dft_r2c;
const auto execute_dft_c2r              = fftw_execute_dft_c2r;
const auto execute_dft_c2c              = fftw_execute_dft;
const auto execute_nd_dft_r2c           = fftw_execute_dft_r2c;
const auto execute_nd_dft_c2r           = fftw_execute_dft_c2r;
const auto execute_nd_dft_c2c           = fftw_execute_dft;
const auto plan_dft_r2c_3d              = fftw_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d              = fftw_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d              = fftw_plan_dft_3d;
const auto plan_dft_c2c_1d              = fftw_plan_dft_1d;
const auto plan_dft_c2r_1d              = fftw_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d              = fftw_plan_dft_r2c_1d;
const auto cleanup                      = fftw_cleanup;
#ifndef SERIAL
using      real_mpi_nd_plan             = fftw_plan;
using      complex_mpi_nd_plan          = fftw_plan;
const auto destroy_real_mpi_nd_plan     = fftw_destroy_plan;
const auto destroy_complex_mpi_nd_plan  = fftw_destroy_plan;
const auto mpi_execute_nd_dft_r2c       = fftw_mpi_execute_dft_r2c;
const auto mpi_execute_nd_dft_c2r       = fftw_mpi_execute_dft_c2r;
const auto mpi_execute_nd_dft_c2c       = fftw_mpi_execute_dft;
const auto mpi_plan_dft_r2c_3d          = fftw_mpi_plan_dft_r2c_3d;
const auto mpi_plan_dft_c2r_3d          = fftw_mpi_plan_dft_c2r_3d;
const auto mpi_plan_dft_c2c_3d          = fftw_mpi_plan_dft_3d;
const auto mpi_cleanup                  = fftw_mpi_cleanup;
#endif // #ifndef SERIAL
};

#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace fftw2 {
using      fft_real                     = real;
using      fft_complex                  = fftw_complex;
using      std_complex                  = std::complex<real>;
using      real_plan                    =   rfftw_plan;
using      real_nd_plan                 = rfftwnd_plan;
using      complex_plan                 =    fftw_plan;
using      complex_nd_plan              =  fftwnd_plan;
const auto fft_malloc                   = malloc;
const auto fft_free                     = free;
const auto destroy_real_plan            =   rfftw_destroy_plan;
const auto destroy_real_nd_plan         = rfftwnd_destroy_plan;
const auto destroy_complex_plan         =    fftw_destroy_plan;
const auto destroy_complex_nd_plan      =  fftwnd_destroy_plan;
const auto execute_dft_r2c              =   rfftw_one;
const auto execute_nd_dft_r2c           = rfftwnd_one_real_to_complex;
const auto execute_dft_c2r              =   rfftw_one;
const auto execute_nd_dft_c2r           = rfftwnd_one_complex_to_real;
const auto execute_dft_c2c              =    fftw_one;
const auto execute_nd_dft_c2c           =  fftwnd_one;
const auto plan_dft_r2c_3d              = rfftw3d_create_plan;
const auto plan_dft_c2r_3d              = rfftw3d_create_plan;
const auto plan_dft_c2c_3d              =  fftw3d_create_plan;
const auto plan_dft_c2c_1d              =  fftw_create_plan;
const auto plan_dft_c2r_1d              = rfftw_create_plan;
const auto plan_dft_r2c_1d              = rfftw_create_plan;
#ifndef SERIAL
using      real_mpi_nd_plan             = rfftwnd_mpi_plan;
using      complex_mpi_nd_plan          =  fftwnd_mpi_plan;
const auto destroy_real_mpi_nd_plan     = rfftwnd_mpi_destroy_plan;
const auto destroy_complex_mpi_nd_plan  =  fftwnd_mpi_destroy_plan;
const auto mpi_execute_nd_dft_r2c       = rfftwnd_mpi;
const auto mpi_execute_nd_dft_c2r       = rfftwnd_mpi;
const auto mpi_execute_nd_dft_c2c       =  fftwnd_mpi;
const auto mpi_plan_dft_r2c_3d          = rfftw3d_mpi_create_plan;
const auto mpi_plan_dft_c2r_3d          = rfftw3d_mpi_create_plan;
const auto mpi_plan_dft_c2c_3d          =  fftw3d_mpi_create_plan;
#endif // #ifndef SERIAL
};
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... # else

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef FLOAT8
namespace gamer_fftw = fftw3_double_precision;
#else // #ifdef FLOAT8
namespace gamer_fftw = fftw3_single_precision;
#endif // #ifdef FLOAT8 ... # else
#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace gamer_fftw = fftw2;
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... #else


//real and imaginary part should only be accessed through c_re and c_im since this macro is also defined in FFTW2
# if ( SUPPORT_FFTW == FFTW3 )
# define c_re(c) ((c)[0])
# define c_im(c) ((c)[1])
# endif

// define index types for mpi_local_size function that uses int in FFTW2 and long int in FFTW3
# if ( SUPPORT_FFTW == FFTW3 )
# define mpi_index_int long int
# else
# define mpi_index_int int
# endif


//wrappers for fftw plans and real-to-complex as well as complex to real n-dimensional transforms on the root level of the AMR hierarchy
//used for Poisson solver and for computing power spectra
namespace root_fftw {
const auto malloc               = gamer_fftw::fft_malloc;
const auto free                 = gamer_fftw::fft_free;
#ifdef SERIAL
using      real_plan            = gamer_fftw::real_nd_plan;
using      complex_plan         = gamer_fftw::complex_nd_plan;
const auto destroy_real_plan    = gamer_fftw::destroy_real_nd_plan;
const auto destroy_complex_plan = gamer_fftw::destroy_complex_nd_plan;
#else // #ifdef SERIAL
using      real_plan            = gamer_fftw::real_mpi_nd_plan;
using      complex_plan         = gamer_fftw::complex_mpi_nd_plan;
const auto destroy_real_plan    = gamer_fftw::destroy_real_mpi_nd_plan;
const auto destroy_complex_plan = gamer_fftw::destroy_complex_mpi_nd_plan;
#endif // #endif
};

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr,                startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr,                startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_FORWARD , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_BACKWARD, startup )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_nd_dft_r2c     ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_nd_dft_c2r     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_nd_dft_c2c     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#else  // #ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_IN  )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_FORWARD , startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_BACKWARD, startup | FFTW_MPI_TRANSPOSED_IN  )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_r2c ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_c2r ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_c2c ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#endif // #ifdef SERIAL ... # else
#else // # if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d        ( size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d        ( size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], FFTW_FORWARD        , startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], FFTW_BACKWARD       , startup | FFTW_IN_PLACE )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_nd_dft_r2c     ( plan,                      (gamer_fftw::fft_real*)    arr, NULL )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_nd_dft_c2r     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_nd_dft_c2c     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#else  // #ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_FORWARD        , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_BACKWARD       , startup )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_r2c ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_c2r ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_nd_dft_c2c ( plan, 1,                   (gamer_fftw::fft_complex*) arr, NULL, FFTW_TRANSPOSED_ORDER )
#endif // #ifdef SERIAL ... # else
#endif // # if ( SUPPORT_FFTW == FFTW3 )

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef FFTW_FLOAT8
namespace gfei_fftw = fftw3_double_precision;
#else // #ifdef GRAMFE_FLOAT8
namespace gfei_fftw = fftw3_single_precision;
#endif // #ifdef GRAMFE_FLOAT8 ... # else
#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace gfei_fftw = fftw2;
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... #else

//wrappers for fftw plans and complex 1D-transform used in Gram-Fourier extension interpolation algorithm
#if ( SUPPORT_FFTW == FFTW3 )
#define gfei_fftw_c2c(plan, arr)                                    gfei_fftw::execute_dft_c2c ( plan, (gfei_fftw::fft_complex*) arr, (gfei_fftw::fft_complex*) arr )
#define gfei_fftw_r2c(plan, arr)                                    gfei_fftw::execute_dft_r2c ( plan, (gfei_fftw::fft_real*)    arr, (gfei_fftw::fft_complex*) arr )
#define gfei_fftw_c2r(plan, arr)                                    gfei_fftw::execute_dft_c2r ( plan, (gfei_fftw::fft_complex*) arr, (gfei_fftw::fft_real*)    arr )
#define gfei_fftw_create_1d_forward_c2c_plan(size, arr)             gfei_fftw::plan_dft_c2c_1d ( size, (gfei_fftw::fft_complex*) arr, (gfei_fftw::fft_complex*) arr, FFTW_FORWARD , FFTW_MEASURE )
#define gfei_fftw_create_1d_backward_c2c_plan(size, arr)            gfei_fftw::plan_dft_c2c_1d ( size, (gfei_fftw::fft_complex*) arr, (gfei_fftw::fft_complex*) arr, FFTW_BACKWARD, FFTW_MEASURE )
#define gfei_fftw_create_1d_r2c_plan(size, arr)                     gfei_fftw::plan_dft_r2c_1d ( size, (gfei_fftw::fft_real*)    arr, (gfei_fftw::fft_complex*) arr, FFTW_MEASURE )
#define gfei_fftw_create_1d_c2r_plan(size, arr)                     gfei_fftw::plan_dft_c2r_1d ( size, (gfei_fftw::fft_complex*) arr, (gfei_fftw::fft_real*)    arr, FFTW_MEASURE )
#else // #if ( SUPPORT_FFTW == FFTW3 )
#define gfei_fftw_c2c(plan, arr)                                    gfei_fftw::execute_dft_c2c ( plan, (gfei_fftw::fft_complex*) arr, NULL )
#define gfei_fftw_r2c(plan, arr)                                    gfei_fftw::execute_dft_r2c ( plan, (gfei_fftw::fft_real*)    arr, NULL )
#define gfei_fftw_c2r(plan, arr)                                    gfei_fftw::execute_dft_c2r ( plan, (gfei_fftw::fft_real*)    arr, NULL )
#define gfei_fftw_create_1d_forward_c2c_plan(size, arr)             gfei_fftw::plan_dft_c2c_1d ( size, FFTW_FORWARD , FFTW_MEASURE | FFTW_IN_PLACE )
#define gfei_fftw_create_1d_backward_c2c_plan(size, arr)            gfei_fftw::plan_dft_c2c_1d ( size, FFTW_BACKWARD, FFTW_MEASURE | FFTW_IN_PLACE )
#define gfei_fftw_create_1d_r2c_plan(size, arr)                     gfei_fftw::plan_dft_r2c_1d ( size, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE | FFTW_IN_PLACE )
#define gfei_fftw_create_1d_c2r_plan(size, arr)                     gfei_fftw::plan_dft_c2r_1d ( size, FFTW_COMPLEX_TO_REAL, FFTW_MEASURE | FFTW_IN_PLACE )
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... # else

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef GRAMFE_FLOAT8
namespace gramfe_fftw = fftw3_double_precision;
#else // #ifdef GRAMFE_FLOAT8
namespace gramfe_fftw = fftw3_single_precision;
#endif // #ifdef GRAMFE_FLOAT8 ... # else
#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace gramfe_fftw = fftw2;
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... #else

//wrappers for fftw plans and complex 1D-transform used in Gram-Fourier extension algorithm
#if ( SUPPORT_FFTW == FFTW3 )
#define gramfe_fftw_c2c(plan, arr)                          gramfe_fftw::execute_1d_dft_c2c ( plan, (gramfe_fftw::fft_complex*) arr, (gramfe_fftw::fft_complex*) arr )
#define gramfe_fftw_r2c(plan, arr)                          gramfe_fftw::execute_1d_dft_r2c ( plan, (gramfe_fftw::fft_real*)    arr, (gramfe_fftw::fft_complex*) arr )
#define gramfe_fftw_c2r(plan, in, out)                      gramfe_fftw::execute_1d_dft_c2r ( plan, (gramfe_fftw::fft_complex*) in,  (gramfe_fftw::fft_real*)    out )
#define gramfe_fftw_create_1d_forward_c2c_plan(size, arr)   gramfe_fftw::plan_dft_c2c_1d ( size, (gramfe_fftw::fft_complex*) arr, (gramfe_fftw::fft_complex*) arr, FFTW_FORWARD , FFTW_PATIENT )
#define gramfe_fftw_create_1d_backward_c2c_plan(size, arr)  gramfe_fftw::plan_dft_c2c_1d ( size, (gramfe_fftw::fft_complex*) arr, (gramfe_fftw::fft_complex*) arr, FFTW_BACKWARD, FFTW_PATIENT )
#else // #if ( SUPPORT_FFTW == FFTW3 )
#define gramfe_fftw_c2c(plan, arr)                          gramfe_fftw::execute_1d_dft_c2c ( plan, (gramfe_fftw::fft_complex*) arr, NULL )
#define gramfe_fftw_r2c(plan, arr)                          gramfe_fftw::execute_1d_dft_r2c ( plan, (gramfe_fftw::fft_real*)    arr, NULL )
#define gramfe_fftw_c2r(plan, in, out)                      gramfe_fftw::execute_1d_dft_c2r ( plan, (gramfe_fftw::fft_complex*) arr, NULL )
#define gramfe_fftw_create_1d_forward_c2c_plan(size, arr)   gramfe_fftw::plan_dft_c2c_1d ( size, FFTW_FORWARD , FFTW_PATIENT | FFTW_IN_PLACE )
#define gramfe_fftw_create_1d_backward_c2c_plan(size, arr)  gramfe_fftw::plan_dft_c2c_1d ( size, FFTW_BACKWARD, FFTW_PATIENT | FFTW_IN_PLACE )
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... # else

#endif  // #if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )

#endif  // #ifndef __FFTW_H__
