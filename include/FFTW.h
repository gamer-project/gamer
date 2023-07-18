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
using      fft_real                     = float;
using      fft_complex                  = fftwf_complex;
using      plan                         = fftwf_plan;
using      real_plan_1d                 = fftwf_plan;
using      real_plan_nd                 = fftwf_plan;
using      complex_plan_1d              = fftwf_plan;
using      complex_plan_nd              = fftwf_plan;
const auto fft_malloc                   = fftwf_malloc;
const auto fft_free                     = fftwf_free;
const auto destroy_plan                 = fftwf_destroy_plan;
const auto destroy_real_plan_1d         = fftwf_destroy_plan;
const auto destroy_real_plan_nd         = fftwf_destroy_plan;
const auto destroy_complex_plan_1d      = fftwf_destroy_plan;
const auto destroy_complex_plan_nd      = fftwf_destroy_plan;
const auto execute_dft_r2c_1d           = fftwf_execute_dft_r2c;
const auto execute_dft_c2r_1d           = fftwf_execute_dft_c2r;
const auto execute_dft_c2c_1d           = fftwf_execute_dft;
const auto execute_dft_r2c_nd           = fftwf_execute_dft_r2c;
const auto execute_dft_c2r_nd           = fftwf_execute_dft_c2r;
const auto execute_dft_c2c_nd           = fftwf_execute_dft;
const auto plan_dft_r2c_3d              = fftwf_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d              = fftwf_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d              = fftwf_plan_dft_3d;
const auto plan_dft_c2c_1d              = fftwf_plan_dft_1d;
const auto plan_dft_c2r_1d              = fftwf_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d              = fftwf_plan_dft_r2c_1d;
const auto cleanup                      = fftwf_cleanup;
#ifndef SERIAL
using      real_mpi_plan_nd             = fftwf_plan;
using      complex_mpi_plan_nd          = fftwf_plan;
const auto destroy_real_mpi_plan_nd     = fftwf_destroy_plan;
const auto destroy_complex_mpi_plan_nd  = fftwf_destroy_plan;
const auto mpi_execute_dft_r2c_nd       = fftwf_mpi_execute_dft_r2c;
const auto mpi_execute_dft_c2r_nd       = fftwf_mpi_execute_dft_c2r;
const auto mpi_execute_dft_c2c_nd       = fftwf_mpi_execute_dft;
const auto mpi_plan_dft_r2c_3d          = fftwf_mpi_plan_dft_r2c_3d;
const auto mpi_plan_dft_c2r_3d          = fftwf_mpi_plan_dft_c2r_3d;
const auto mpi_plan_dft_c2c_3d          = fftwf_mpi_plan_dft_3d;
const auto mpi_cleanup                  = fftwf_mpi_cleanup;
#endif // #ifndef SERIAL
};

namespace fftw3_double_precision {
using      fft_real                     = double;
using      fft_complex                  = fftw_complex;
using      plan                         = fftw_plan;
using      real_plan_1d                 = fftw_plan;
using      real_plan_nd                 = fftw_plan;
using      complex_plan_1d              = fftw_plan;
using      complex_plan_nd              = fftw_plan;
const auto fft_malloc                   = fftw_malloc;
const auto fft_free                     = fftw_free;
const auto destroy_plan_1d              = fftw_destroy_plan;
const auto destroy_real_plan_1d         = fftw_destroy_plan;
const auto destroy_real_plan_nd         = fftw_destroy_plan;
const auto destroy_complex_plan_1d      = fftw_destroy_plan;
const auto destroy_complex_plan_nd      = fftw_destroy_plan;
const auto execute_dft_r2c_1d           = fftw_execute_dft_r2c;
const auto execute_dft_c2r_1d           = fftw_execute_dft_c2r;
const auto execute_dft_c2c_1d           = fftw_execute_dft;
const auto execute_dft_r2c_nd           = fftw_execute_dft_r2c;
const auto execute_dft_c2r_nd           = fftw_execute_dft_c2r;
const auto execute_dft_c2c_nd           = fftw_execute_dft;
const auto plan_dft_r2c_3d              = fftw_plan_dft_r2c_3d;
const auto plan_dft_c2r_3d              = fftw_plan_dft_c2r_3d;
const auto plan_dft_c2c_3d              = fftw_plan_dft_3d;
const auto plan_dft_c2c_1d              = fftw_plan_dft_1d;
const auto plan_dft_c2r_1d              = fftw_plan_dft_c2r_1d;
const auto plan_dft_r2c_1d              = fftw_plan_dft_r2c_1d;
const auto cleanup                      = fftw_cleanup;
#ifndef SERIAL
using      real_mpi_plan_nd             = fftw_plan;
using      complex_mpi_plan_nd          = fftw_plan;
const auto destroy_real_mpi_plan_nd     = fftw_destroy_plan;
const auto destroy_complex_mpi_plan_nd  = fftw_destroy_plan;
const auto mpi_execute_dft_r2c_nd       = fftw_mpi_execute_dft_r2c;
const auto mpi_execute_dft_c2r_nd       = fftw_mpi_execute_dft_c2r;
const auto mpi_execute_dft_c2c_nd       = fftw_mpi_execute_dft;
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
using      real_plan_1d                 =   rfftw_plan;
using      real_plan_nd                 = rfftwnd_plan;
using      complex_plan_1d              =    fftw_plan;
using      complex_plan_nd              =  fftwnd_plan;
const auto fft_malloc                   = malloc;
const auto fft_free                     = free;
const auto destroy_real_plan_1d         =   rfftw_destroy_plan;
const auto destroy_real_plan_nd         = rfftwnd_destroy_plan;
const auto destroy_complex_plan_1d      =    fftw_destroy_plan;
const auto destroy_complex_plan_nd      =  fftwnd_destroy_plan;
const auto execute_dft_r2c_1d           =   rfftw_one;
const auto execute_dft_c2r_1d           =   rfftw_one;
const auto execute_dft_c2c_1d           =    fftw_one;
const auto execute_dft_r2c_nd           = rfftwnd_one_real_to_complex;
const auto execute_dft_c2r_nd           = rfftwnd_one_complex_to_real;
const auto execute_dft_c2c_nd           =  fftwnd_one;
const auto plan_dft_r2c_3d              = rfftw3d_create_plan;
const auto plan_dft_c2r_3d              = rfftw3d_create_plan;
const auto plan_dft_c2c_3d              =  fftw3d_create_plan;
const auto plan_dft_c2c_1d              =  fftw_create_plan;
const auto plan_dft_c2r_1d              = rfftw_create_plan;
const auto plan_dft_r2c_1d              = rfftw_create_plan;
#ifndef SERIAL
using      real_mpi_plan_nd             = rfftwnd_mpi_plan;
using      complex_mpi_plan_nd          =  fftwnd_mpi_plan;
const auto destroy_real_mpi_plan_nd     = rfftwnd_mpi_destroy_plan;
const auto destroy_complex_mpi_plan_nd  =  fftwnd_mpi_destroy_plan;
const auto mpi_execute_dft_r2c_nd       = rfftwnd_mpi;
const auto mpi_execute_dft_c2r_nd       = rfftwnd_mpi;
const auto mpi_execute_dft_c2c_nd       =  fftwnd_mpi;
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
const auto fft_malloc              = gamer_fftw::fft_malloc;
const auto fft_free                = gamer_fftw::fft_free;
#ifdef SERIAL
using      real_plan_nd            = gamer_fftw::real_plan_nd;
using      complex_plan_nd         = gamer_fftw::complex_plan_nd;
const auto destroy_real_plan_nd    = gamer_fftw::destroy_real_plan_nd;
const auto destroy_complex_plan_nd = gamer_fftw::destroy_complex_plan_nd;
#else // #ifdef SERIAL
using      real_plan_nd            = gamer_fftw::real_mpi_plan_nd;
using      complex_plan_nd         = gamer_fftw::complex_mpi_plan_nd;
const auto destroy_real_plan_nd    = gamer_fftw::destroy_real_mpi_plan_nd;
const auto destroy_complex_plan_nd = gamer_fftw::destroy_complex_mpi_plan_nd;
#endif // #endif
};

#if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr,                startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr,                startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_FORWARD , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, FFTW_BACKWARD, startup )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_dft_r2c_nd     ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_dft_c2r_nd     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_dft_c2c_nd     ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#else  // #ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_IN  )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_FORWARD , startup | FFTW_MPI_TRANSPOSED_OUT )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d    ( size[2], size[1], size[0], (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr, MPI_COMM_WORLD, FFTW_BACKWARD, startup | FFTW_MPI_TRANSPOSED_IN  )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_r2c_nd ( plan,                      (gamer_fftw::fft_real*)    arr, (gamer_fftw::fft_complex*) arr )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2r_nd ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_real*)    arr )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2c_nd ( plan,                      (gamer_fftw::fft_complex*) arr, (gamer_fftw::fft_complex*) arr )
#endif // #ifdef SERIAL ... # else
#else // # if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::plan_dft_r2c_3d        ( size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::plan_dft_c2r_3d        ( size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], FFTW_FORWARD        , startup | FFTW_IN_PLACE )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::plan_dft_c2c_3d        ( size[2], size[1], size[0], FFTW_BACKWARD       , startup | FFTW_IN_PLACE )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::execute_dft_r2c_nd     ( plan,                      (gamer_fftw::fft_real*)    arr, NULL )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::execute_dft_c2r_nd     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::execute_dft_c2c_nd     ( plan,                      (gamer_fftw::fft_complex*) arr, NULL )
#else  // #ifdef SERIAL
#define root_fftw_create_3d_r2c_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_r2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup )
#define root_fftw_create_3d_c2r_plan(size, arr, startup)            gamer_fftw::mpi_plan_dft_c2r_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup )
#define root_fftw_create_3d_forward_c2c_plan(size, arr, startup)    gamer_fftw::mpi_plan_dft_c2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_FORWARD        , startup )
#define root_fftw_create_3d_backward_c2c_plan(size, arr, startup)   gamer_fftw::mpi_plan_dft_c2c_3d    ( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_BACKWARD       , startup )
#define root_fftw_r2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_r2c_nd ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2r(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2r_nd ( plan, 1,                   (gamer_fftw::fft_real*)    arr, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2c(plan, arr)                                    gamer_fftw::mpi_execute_dft_c2c_nd ( plan, 1,                   (gamer_fftw::fft_complex*) arr, NULL, FFTW_TRANSPOSED_ORDER )
#endif // #ifdef SERIAL ... # else
<<<<<<< HEAD
#endif // # if ( SUPPORT_FFTW == FFTW3 )

#if ( WAVE_SCHEME == WAVE_GRAMFE && SUPPORT_FFTW == FFTW3 )
#ifdef GRAMFE_FLOAT8
#define gramfe_float_fftw3_malloc              fftw_malloc
#define gramfe_float_fftw3_free                fftw_free
#define gramfe_float_fftw3_plan                fftw_plan
#define gramfe_float_fftw3_destroy_plan        fftw_destroy_plan
#define gramfe_float_fftw3_complex             fftw_complex
#define gramfe_float_fftw3_execute_dft_c2c     fftw_execute_dft
#define gramfe_float_fftw3_execute_dft_r2c     fftw_execute_dft_r2c
#define gramfe_float_fftw3_execute_dft_c2r     fftw_execute_dft_c2r
#define gramfe_float_fftw3_plan_dft_c2c_1d     fftw_plan_dft_1d
#define gramfe_float_fftw3_plan_dft_c2r_1d     fftw_plan_dft_c2r_1d
#define gramfe_float_fftw3_plan_dft_r2c_1d     fftw_plan_dft_r2c_1d
#else // # ifdef GRAMFE_FLOAT8
#define gramfe_float_fftw3_malloc              fftwf_malloc
#define gramfe_float_fftw3_free                fftwf_free
#define gramfe_float_fftw3_plan                fftwf_plan
#define gramfe_float_fftw3_destroy_plan        fftwf_destroy_plan
#define gramfe_float_fftw3_complex             fftwf_complex
#define gramfe_float_fftw3_execute_dft_c2c     fftwf_execute_dft
#define gramfe_float_fftw3_execute_dft_r2c     fftwf_execute_dft_r2c
#define gramfe_float_fftw3_execute_dft_c2r     fftwf_execute_dft_c2r
#define gramfe_float_fftw3_plan_dft_c2c_1d     fftwf_plan_dft_1d
#define gramfe_float_fftw3_plan_dft_c2r_1d     fftwf_plan_dft_c2r_1d
#define gramfe_float_fftw3_plan_dft_r2c_1d     fftwf_plan_dft_r2c_1d
#endif // #ifdef GRAMFE_FLOAT8 ... # else
#endif // # if ( WAVE_SCHEME == WAVE_GRAMFE && SUPPORT_FFTW == FFTW3 )

//wrappers for fftw plans and complex 1D-transform used in Gram-Fourier extension algorithm
#if ( SUPPORT_FFTW == FFTW3 )
#define gramfe_float_complex                                        gramfe_float_fftw3_complex
#define gramfe_real_fftw_plan                                       gramfe_float_fftw3_plan
#define gramfe_complex_fftw_plan                                    gramfe_float_fftw3_plan
#define gramfe_fftw_malloc                                          gramfe_float_fftw3_malloc
#define gramfe_fftw_free                                            gramfe_float_fftw3_free
#define gramfe_fftw_c2c(plan, array)                                gramfe_float_fftw3_execute_dft_c2c    ( plan, (gramfe_float_complex*) array, (gramfe_float_complex*)  array )
#define gramfe_fftw_r2c(plan, array)                                gramfe_float_fftw3_execute_dft_r2c    ( plan, (gramfe_float*)         array, (gramfe_float_complex*)  array )
#define gramfe_fftw_c2r(plan, in, out)                              gramfe_float_fftw3_execute_dft_c2r    ( plan, (gramfe_float_complex*) in,    (gramfe_float*)          out   )
#define gramfe_create_fftw_1d_forward_c2c_plan(size, arr, startup)  gramfe_float_fftw3_plan_dft_c2c_1d(     size, (gramfe_float_fftw3_complex*) arr, (gramfe_float_fftw3_complex*) arr, FFTW_FORWARD , startup )
#define gramfe_create_fftw_1d_backward_c2c_plan(size, arr, startup) gramfe_float_fftw3_plan_dft_c2c_1d(     size, (gramfe_float_fftw3_complex*) arr, (gramfe_float_fftw3_complex*) arr, FFTW_BACKWARD, startup )
#define gramfe_destroy_complex_fftw_plan                            gramfe_float_fftw3_destroy_plan
#else // #if ( SUPPORT_FFTW == FFTW3 )
#define gramfe_float_complex                                        fftw_complex
#define gramfe_fftw_malloc                                          malloc
#define gramfe_fftw_free                                            free
#define gramfe_complex_fftw_plan                                    fftw_plan
#define gramfe_fftw_c2c(plan, array)                                fftw_one                   ( plan, (gramfe_float_complex*) array, NULL )
#define gramfe_fftw_r2c(plan, array)                                rfftw_one_real_to_complex  ( plan, (gramfe_float*)         array, NULL )
#define gramfe_fftw_c2r(plan, in, out)                              rfftw_one_complex_to_real  ( plan, (gamer_float_complex*)  array, NULL )
#define gramfe_create_fftw_1d_forward_c2c_plan(size, arr, startup)  fftw_create_plan( size, FFTW_FORWARD , startup | FFTW_IN_PLACE )
#define gramfe_create_fftw_1d_backward_c2c_plan(size, arr, startup) fftw_create_plan( size, FFTW_BACKWARD, startup | FFTW_IN_PLACE )
#define gramfe_destroy_complex_fftw_plan                            fftw_destroy_plan
#endif // #if ( SUPPORT_FFTW == FFTW3 ) ... # else



#endif  // # if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )
#endif  // #ifndef __FFTW_H__
=======
#endif // # if ( SUPPORT_FFTW == FFTW3 )  ... # else

#endif  // #if ( SUPPORT_FFTW == FFTW2 || SUPPORT_FFTW == FFTW3 )

#endif  // #ifndef __FFTW_H__
>>>>>>> origin/main
