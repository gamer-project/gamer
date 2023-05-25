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
#ifdef FLOAT8
#define gamer_float_fftw3_malloc              fftw_malloc
#define gamer_float_fftw3_free                fftw_free
#define gamer_float_fftw3_plan                fftw_plan
#define gamer_float_fftw3_destroy_plan        fftw_destroy_plan
#define gamer_float_fftw3_complex             fftw_complex
#define gamer_float_fftw3_execute_dft_r2c     fftw_execute_dft_r2c
#define gamer_float_fftw3_execute_dft_c2r     fftw_execute_dft_c2r
#define gamer_float_fftw3_execute_dft_c2c     fftw_execute_dft
#define gamer_float_fftw3_plan_dft_r2c_3d     fftw_plan_dft_r2c_3d
#define gamer_float_fftw3_plan_dft_c2r_3d     fftw_plan_dft_c2r_3d
#define gamer_float_fftw3_plan_dft_c2c_3d     fftw_plan_dft_3d
#define gamer_float_fftw3_cleanup             fftw_cleanup
#ifndef SERIAL
#define gamer_float_fftw3_mpi_execute_dft_r2c fftw_mpi_execute_dft_r2c
#define gamer_float_fftw3_mpi_execute_dft_c2r fftw_mpi_execute_dft_c2r
#define gamer_float_fftw3_mpi_execute_dft_c2c fftw_mpi_execute_dft
#define gamer_float_fftw3_mpi_plan_dft_r2c_3d fftw_mpi_plan_dft_r2c_3d
#define gamer_float_fftw3_mpi_plan_dft_c2r_3d fftw_mpi_plan_dft_c2r_3d
#define gamer_float_fftw3_mpi_plan_dft_c2c_3d fftw_mpi_plan_dft_3d
#define gamer_float_fftw3_mpi_cleanup         fftw_mpi_cleanup
#endif // #ifndef SERIAL
#else // #ifdef FLOAT8
#define gamer_float_fftw3_malloc              fftwf_malloc
#define gamer_float_fftw3_free                fftwf_free
#define gamer_float_fftw3_plan                fftwf_plan
#define gamer_float_fftw3_destroy_plan        fftwf_destroy_plan
#define gamer_float_fftw3_complex             fftwf_complex
#define gamer_float_fftw3_execute_dft_r2c     fftwf_execute_dft_r2c
#define gamer_float_fftw3_execute_dft_c2r     fftwf_execute_dft_c2r
#define gamer_float_fftw3_execute_dft_c2c     fftwf_execute_dft
#define gamer_float_fftw3_plan_dft_r2c_3d     fftwf_plan_dft_r2c_3d
#define gamer_float_fftw3_plan_dft_c2r_3d     fftwf_plan_dft_c2r_3d
#define gamer_float_fftw3_plan_dft_c2c_3d     fftwf_plan_dft_3d
#define gamer_float_fftw3_cleanup             fftwf_cleanup
#ifndef SERIAL
#define gamer_float_fftw3_mpi_execute_dft_r2c fftwf_mpi_execute_dft_r2c
#define gamer_float_fftw3_mpi_execute_dft_c2r fftwf_mpi_execute_dft_c2r
#define gamer_float_fftw3_mpi_execute_dft_c2c fftwf_mpi_execute_dft
#define gamer_float_fftw3_mpi_plan_dft_r2c_3d fftwf_mpi_plan_dft_r2c_3d
#define gamer_float_fftw3_mpi_plan_dft_c2r_3d fftwf_mpi_plan_dft_c2r_3d
#define gamer_float_fftw3_mpi_plan_dft_c2c_3d fftwf_mpi_plan_dft_3d
#define gamer_float_fftw3_mpi_cleanup         fftwf_mpi_cleanup
#endif // #ifndef SERIAL
#endif // #ifdef FLOAT8 ... else
#endif // #if ( SUPPORT_FFTW == FFTW3 )


//wrapper for fftw complex type
//note that real and imaginary part should only be accessed through c_re and c_im since this macro is also defined in FFTW2
#if ( SUPPORT_FFTW == FFTW3 )
#define gamer_float_complex gamer_float_fftw3_complex
#define c_re(c) ((c)[0])
#define c_im(c) ((c)[1])
#else // # if ( SUPPORT_FFTW == FFTW3 )
#define gamer_float_complex fftw_complex
#endif // # if ( SUPPORT_FFTW == FFTW3 )  ... # else

// define index types for mpi_local_size function that uses int in FFTW2 and long int in FFTW3
# if ( SUPPORT_FFTW == FFTW3 )
# define mpi_index_int long int
# else
# define mpi_index_int int
# endif

//wrappers for fftw plans and real-to-complex as well as complex to real n-dimensional transforms on the root level of the AMR hierarchy
//used for Poisson solver and for computing power spectra
#if ( SUPPORT_FFTW == FFTW3 )
#define root_real_fftw_plan         gamer_float_fftw3_plan
#define root_complex_fftw_plan      gamer_float_fftw3_plan
#define root_fftw_malloc            gamer_float_fftw3_malloc
#define root_fftw_free              gamer_float_fftw3_free
#ifdef SERIAL
#define root_fftw_r2c(plan, array)  gamer_float_fftw3_execute_dft_r2c    ( plan, (real*)                array, (gamer_float_complex*)  array )
#define root_fftw_c2r(plan, array)  gamer_float_fftw3_execute_dft_c2r    ( plan, (gamer_float_complex*) array, (real*)                 array )
#define root_fftw_c2c(plan, array)  gamer_float_fftw3_execute_dft_c2c    ( plan, (gamer_float_complex*) array, (gamer_float_complex*)  array )
#else  // #ifdef SERIAL
#define root_fftw_r2c(plan, array)  gamer_float_fftw3_mpi_execute_dft_r2c( plan, (real*)                array, (gamer_float_complex*)  array )
#define root_fftw_c2r(plan, array)  gamer_float_fftw3_mpi_execute_dft_c2r( plan, (gamer_float_complex*) array, (real*)                 array )
#define root_fftw_c2c(plan, array)  gamer_float_fftw3_mpi_execute_dft_c2c( plan, (gamer_float_complex*) array, (gamer_float_complex*)  array )
#endif // #ifdef SERIAL ... # else
#else // # if ( SUPPORT_FFTW == FFTW3 )
#define root_fftw_malloc            malloc
#define root_fftw_free              free
#ifdef SERIAL
#define root_real_fftw_plan         rfftwnd_plan
#define root_complex_fftw_plan       fftwnd_plan
#define root_fftw_r2c(plan, array)  rfftwnd_one_real_to_complex( plan, (real*)                array, NULL )
#define root_fftw_c2r(plan, array)  rfftwnd_one_complex_to_real( plan, (gamer_float_complex*) array, NULL )
#define root_fftw_c2c(plan, array)   fftwnd_one                ( plan, (gamer_float_complex*) array, NULL )
#else  // #ifdef SERIAL
#define root_real_fftw_plan         rfftwnd_mpi_plan
#define root_complex_fftw_plan       fftwnd_mpi_plan
#define root_fftw_r2c(plan, array)  rfftwnd_mpi( plan, 1, (real*)                array, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2r(plan, array)  rfftwnd_mpi( plan, 1, (real*)                array, NULL, FFTW_TRANSPOSED_ORDER )
#define root_fftw_c2c(plan, array)   fftwnd_mpi( plan, 1, (gamer_float_complex*) array, NULL, FFTW_TRANSPOSED_ORDER )
#endif // #ifdef SERIAL ... # else
#endif // # if ( SUPPORT_FFTW == FFTW3 )  ... # else

//wrappers for fftw create and destroy plan functions used in Init_FFTW
#if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define create_fftw_3d_r2c_plan(size, arr, startup)          gamer_float_fftw3_plan_dft_r2c_3d(     size[2], size[1], size[0], (real*)                      arr, (gamer_float_fftw3_complex*) arr,                startup )
#define create_fftw_3d_c2r_plan(size, arr, startup)          gamer_float_fftw3_plan_dft_c2r_3d(     size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (real*)                      arr,                startup )
#define create_fftw_3d_forward_c2c_plan(size, arr, startup)  gamer_float_fftw3_plan_dft_c2c_3d(     size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (gamer_float_fftw3_complex*) arr, FFTW_FORWARD , startup )
#define create_fftw_3d_backward_c2c_plan(size, arr, startup) gamer_float_fftw3_plan_dft_c2c_3d(     size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (gamer_float_fftw3_complex*) arr, FFTW_BACKWARD, startup )
#define destroy_real_fftw_plan                               gamer_float_fftw3_destroy_plan
#define destroy_complex_fftw_plan                            gamer_float_fftw3_destroy_plan
#else  // #ifdef SERIAL
#define create_fftw_3d_r2c_plan(size, arr, startup)          gamer_float_fftw3_mpi_plan_dft_r2c_3d( size[2], size[1], size[0], (real*)                      arr, (gamer_float_fftw3_complex*) arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_OUT )
#define create_fftw_3d_c2r_plan(size, arr, startup)          gamer_float_fftw3_mpi_plan_dft_c2r_3d( size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (real*)                      arr, MPI_COMM_WORLD,                startup | FFTW_MPI_TRANSPOSED_IN  )
#define create_fftw_3d_forward_c2c_plan(size, arr, startup)  gamer_float_fftw3_mpi_plan_dft_c2c_3d( size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (gamer_float_fftw3_complex*) arr, MPI_COMM_WORLD, FFTW_FORWARD , startup | FFTW_MPI_TRANSPOSED_OUT )
#define create_fftw_3d_backward_c2c_plan(size, arr, startup) gamer_float_fftw3_mpi_plan_dft_c2c_3d( size[2], size[1], size[0], (gamer_float_fftw3_complex*) arr, (gamer_float_fftw3_complex*) arr, MPI_COMM_WORLD, FFTW_BACKWARD, startup | FFTW_MPI_TRANSPOSED_IN  )
#define destroy_real_fftw_plan                               gamer_float_fftw3_destroy_plan
#define destroy_complex_fftw_plan                            gamer_float_fftw3_destroy_plan
#endif // #ifdef SERIAL ... # else
#else // # if ( SUPPORT_FFTW == FFTW3 )
#ifdef SERIAL
#define create_fftw_3d_r2c_plan(size, arr, startup)          rfftw3d_create_plan(                     size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup | FFTW_IN_PLACE )
#define create_fftw_3d_c2r_plan(size, arr, startup)          rfftw3d_create_plan(                     size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup | FFTW_IN_PLACE )
#define create_fftw_3d_forward_c2c_plan(size, arr, startup)   fftw3d_create_plan(                     size[2], size[1], size[0], FFTW_FORWARD        , startup | FFTW_IN_PLACE )
#define create_fftw_3d_backward_c2c_plan(size, arr, startup)  fftw3d_create_plan(                     size[2], size[1], size[0], FFTW_BACKWARD       , startup | FFTW_IN_PLACE )
#define destroy_real_fftw_plan                               rfftwnd_destroy_plan
#define destroy_complex_fftw_plan                             fftwnd_destroy_plan
#else  // #ifdef SERIAL
#define create_fftw_3d_r2c_plan(size, arr, startup)          rfftw3d_mpi_create_plan( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_REAL_TO_COMPLEX, startup )
#define create_fftw_3d_c2r_plan(size, arr, startup)          rfftw3d_mpi_create_plan( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_COMPLEX_TO_REAL, startup )
#define create_fftw_3d_forward_c2c_plan(size, arr, startup)   fftw3d_mpi_create_plan( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_FORWARD        , startup )
#define create_fftw_3d_backward_c2c_plan(size, arr, startup)  fftw3d_mpi_create_plan( MPI_COMM_WORLD, size[2], size[1], size[0], FFTW_BACKWARD       , startup )
#define destroy_real_fftw_plan                               rfftwnd_mpi_destroy_plan
#define destroy_complex_fftw_plan                             fftwnd_mpi_destroy_plan
#endif // #ifdef SERIAL ... # else
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
