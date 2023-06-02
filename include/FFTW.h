#ifndef __FFTW_H__
#define __FFTW_H__



// ****************************************************************************
// ** This header defines wrappers for functions in the FFTW 2 and 3         **
// ** libraries in single and double precision.                              **
// ****************************************************************************


#include "GAMER.h"


#if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )
// wrappers for fftw3 single and double precision routines
#if ( SUPPORT_FFTW == FFTW3 )

namespace fftw3_single_precision {
using      fft_real                 = float;
using      fft_complex              = fftwf_complex;
using      std_complex              = std::complex<float>;
using      plan                     = fftwf_plan;
using      real_plan                = fftwf_plan;
using      complex_plan             = fftwf_plan;
const auto fft_malloc               = fftwf_malloc;
const auto fft_free                 = fftwf_free;
const auto destroy_plan             = fftwf_destroy_plan;
const auto destroy_real_plan        = fftwf_destroy_plan;
const auto destroy_complex_plan     = fftwf_destroy_plan;
const auto execute_dft_r2c          = fftwf_execute_dft_r2c;
const auto execute_dft_c2r          = fftwf_execute_dft_c2r;
const auto execute_dft_c2c          = fftwf_execute_dft;
const auto plan_dft_r2c_3d          = fftwf_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d          = fftwf_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d          = fftwf_plan_dft_3d;
const auto plan_dft_c2c_1d          = fftwf_plan_dft_1d;
const auto plan_dft_c2r_1d          = fftwf_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d          = fftwf_plan_dft_r2c_1d;
const auto cleanup                  = fftwf_cleanup;
#ifndef SERIAL
using      real_mpi_plan            = fftwf_plan;
using      complex_mpi_plan         = fftwf_plan;
const auto destroy_real_mpi_plan    = fftwf_destroy_plan;
const auto destroy_complex_mpi_plan = fftwf_destroy_plan;
const auto mpi_execute_dft_r2c      = fftwf_mpi_execute_dft_r2c;
const auto mpi_execute_dft_c2r      = fftwf_mpi_execute_dft_c2r;
const auto mpi_execute_dft_c2c      = fftwf_mpi_execute_dft;
const auto mpi_plan_dft_r2c_3d      = fftwf_mpi_plan_dft_r2c_3d;
const auto mpi_plan_dft_c2r_3d      = fftwf_mpi_plan_dft_c2r_3d;
const auto mpi_plan_dft_c2c_3d      = fftwf_mpi_plan_dft_3d;
const auto mpi_cleanup              = fftwf_mpi_cleanup;
#endif // #ifndef SERIAL
};

namespace fftw3_double_precision {
using      fft_real                 = double;
using      fft_complex              = fftw_complex;
using      std_complex              = std::complex<double>;
using      real_plan                = fftw_plan;
using      complex_plan             = fftw_plan;
const auto fft_malloc               = fftw_malloc;
const auto fft_free                 = fftw_free;
const auto destroy_plan             = fftw_destroy_plan;
const auto destroy_real_plan        = fftw_destroy_plan;
const auto destroy_complex_plan     = fftw_destroy_plan;
const auto execute_dft_r2c          = fftw_execute_dft_r2c;
const auto execute_dft_c2r          = fftw_execute_dft_c2r;
const auto execute_dft_c2c          = fftw_execute_dft;
const auto plan_dft_r2c_3d          = fftw_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d          = fftw_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d          = fftw_plan_dft_3d;
const auto plan_dft_c2c_1d          = fftw_plan_dft_1d;
const auto plan_dft_c2r_1d          = fftw_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d          = fftw_plan_dft_r2c_1d;
const auto cleanup                  = fftw_cleanup;
#ifndef SERIAL
using      real_mpi_plan            = fftw_plan;
using      complex_mpi_plan         = fftw_plan;
const auto destroy_real_mpi_plan    = fftw_destroy_plan;
const auto destroy_complex_mpi_plan = fftw_destroy_plan;
const auto mpi_execute_dft_r2c      = fftw_mpi_execute_dft_r2c;
const auto mpi_execute_dft_c2r      = fftw_mpi_execute_dft_c2r;
const auto mpi_execute_dft_c2c      = fftw_mpi_execute_dft;
const auto mpi_plan_dft_r2c_3d      = fftw_mpi_plan_dft_r2c_3d;
const auto mpi_plan_dft_c2r_3d      = fftw_mpi_plan_dft_c2r_3d;
const auto mpi_plan_dft_c2c_3d      = fftw_mpi_plan_dft_3d;
const auto mpi_cleanup              = fftw_mpi_cleanup;
#endif // #ifndef SERIAL
};

#else // #if ( SUPPORT_FFTW == FFTW3 )
namespace fftw2 {
using      fft_real                 = real;
using      fft_complex              = fftw_complex;
using      std_complex              = std::complex<real>;
using      real_plan                = rfftwnd_plan;
using      complex_plan             =  fftwnd_plan;
const auto fft_malloc               = malloc;
const auto fft_free                 = free;
const auto destroy_real_plan        = rfftwnd_destroy_plan;
const auto destroy_complex_plan     =  fftwnd_destroy_plan;
const auto execute_dft_r2c          = rfftwnd_one_real_to_complex;
const auto execute_dft_c2r          = rfftwnd_one_complex_to_real;
const auto execute_dft_c2c          =  fftwnd_one;
const auto plan_dft_r2c_3d          = rfftw3d_create_plan;
const auto plan_dft_c2r_3d          = rfftw3d_create_plan;
const auto plan_dft_c2c_3d          =  fftw3d_create_plan;
const auto plan_dft_c2c_1d          =  fftw_create_plan;
const auto plan_dft_c2r_1d          = rfftw_create_plan;
const auto plan_dft_r2c_1d          = rfftw_create_plan;
#ifndef SERIAL
using      real_mpi_plan            = rfftwnd_mpi_plan;
using      complex_mpi_plan         =  fftwnd_mpi_plan;
const auto destroy_real_mpi_plan    = rfftwnd_mpi_destroy_plan;
const auto destroy_complex_mpi_plan =  fftwnd_mpi_destroy_plan;
const auto mpi_execute_dft_r2c      = rfftwnd_mpi;
const auto mpi_execute_dft_c2r      = rfftwnd_mpi;
const auto mpi_execute_dft_c2c      =  fftwnd_mpi;
const auto mpi_plan_dft_r2c_3d      = rfftw3d_mpi_create_plan;
const auto mpi_plan_dft_c2r_3d      = rfftw3d_mpi_create_plan;
const auto mpi_plan_dft_c2c_3d      =  fftw3d_mpi_create_plan;
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
using      real_plan            = gamer_fftw::real_plan;
using      complex_plan         = gamer_fftw::fft_complex_plan;
const auto destroy_real_plan    = gamer_fftw::destroy_real_plan;
const auto destroy_complex_plan = gamer_fftw::destroy_complex_plan;
#else // #ifdef SERIAL
using      real_plan            = gamer_fftw::real_mpi_plan;
using      complex_plan         = gamer_fftw::fft_complex_mpi_plan;
const auto destroy_real_plan    = gamer_fftw::destroy_real_mpi_plan;
const auto destroy_complex_plan = gamer_fftw::destroy_complex_mpi_plan;
#endif // #endif
};

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_dft_r2c     ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_dft_c2r     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_dft_c2c     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d     ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr,                startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d     ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr,                startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d     ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_FORWARD , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d     ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_BACKWARD, startup )
#else  // #ifdef SERIAL
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_r2c ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2r ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2c ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_IN  )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_FORWARD , startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_BACKWARD, startup | FFTW_MPI_TRANSPOSED_IN  )
#endif // #ifdef SERIAL ... # else
#else // # if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d     ( size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d     ( size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d     ( size[2], size[1], size[0], FFTW_FORWARD        , startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d     ( size[2], size[1], size[0], FFTW_BACKWARD       , startup | FFTW_IN_PLACE )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_dft_r2c     ( plan,                      (gamer_fftw::fft_real*)    arr, NULL )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_dft_c2r     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_dft_c2c     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#else  // #ifdef SERIAL
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_r2c ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2r ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2c ( plan, 1,                   (gamer_fftw::fft_complex*) arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_FORWARD        , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_BACKWARD       , startup )
#endif // #ifdef SERIAL ... # else
#endif // # if ( SUPPORT_FFTW == FFTW3 )  ... # else

#endif  // #if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )

#endif  // #ifndef __FFTW_H__
