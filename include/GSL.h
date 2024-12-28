#ifndef __GAMER_GSL_H__
#define __GAMER_GSL_H__

#ifdef SUPPORT_GSL



// ****************************************************************************
// ** This header defines wrappers for functions in the GSL library          **
// ** in single and double precision.                                        **
// ****************************************************************************


#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_permute_matrix.h>
#include <gsl/gsl_linalg.h>


namespace gsl_double_precision {
    using      gsl_real                        = double;
    using      complex                         = gsl_complex;
    using      vector                          = gsl_vector;
    using      vector_view                     = gsl_vector_view;
    using      vector_const_view               = gsl_vector_const_view;
    using      vector_complex                  = gsl_vector_complex;
    using      vector_complex_view             = gsl_vector_complex_view;
    using      vector_complex_const_view       = gsl_vector_complex_const_view;
    using      matrix                          = gsl_matrix;
    using      matrix_view                     = gsl_matrix_view;
    using      matrix_const_view               = gsl_matrix_const_view;
    using      matrix_complex                  = gsl_matrix_complex;
    using      matrix_complex_view             = gsl_matrix_complex_view;
    using      matrix_complex_const_view       = gsl_matrix_complex_const_view;
    using      permutation                     = gsl_permutation;
    const auto sf_gegenpoly_n                  = gsl_sf_gegenpoly_n;
    const auto vector_alloc                    = gsl_vector_alloc;
    const auto vector_free                     = gsl_vector_free;
    const auto vector_view_array               = gsl_vector_view_array;
    const auto vector_const_view_array         = gsl_vector_const_view_array;
    const auto vector_const_subvector          = gsl_vector_const_subvector;
    const auto vector_complex_alloc            = gsl_vector_complex_alloc;
    const auto vector_complex_free             = gsl_vector_complex_free;
    const auto vector_complex_view_array       = gsl_vector_complex_view_array;
    const auto vector_complex_const_view_array = gsl_vector_complex_const_view_array;
    const auto vector_complex_const_subvector  = gsl_vector_complex_const_subvector;
    const auto matrix_alloc                    = gsl_matrix_alloc;
    const auto matrix_free                     = gsl_matrix_free;
    const auto matrix_view_array               = gsl_matrix_view_array;
    const auto matrix_const_view_array         = gsl_matrix_const_view_array;
    const auto matrix_const_submatrix          = gsl_matrix_const_submatrix;
    const auto matrix_transpose                = gsl_matrix_transpose;
    const auto matrix_complex_alloc            = gsl_matrix_complex_alloc;
    const auto matrix_complex_free             = gsl_matrix_complex_free;
    const auto matrix_complex_view_array       = gsl_matrix_complex_view_array;
    const auto matrix_complex_const_view_array = gsl_matrix_complex_const_view_array;
    const auto matrix_complex_const_submatrix  = gsl_matrix_complex_const_submatrix;
    const auto matrix_complex_transpose        = gsl_matrix_complex_transpose;
    const auto permutation_alloc               = gsl_permutation_alloc;
    const auto permutation_free                = gsl_permutation_free;
    const auto permute_vector                  = gsl_permute_vector;
    const auto permute_vector_complex          = gsl_permute_vector_complex;
    const auto linalg_complex_LU_decomp        = gsl_linalg_complex_LU_decomp;
    const auto linalg_complex_LU_svx           = gsl_linalg_complex_LU_svx;
    const auto vector_set                      = gsl_vector_set;
    const auto vector_get                      = gsl_vector_get;
    const auto vector_complex_set              = gsl_vector_complex_set;
    const auto vector_complex_get              = gsl_vector_complex_get;
    const auto matrix_set                      = gsl_matrix_set;
    const auto matrix_get                      = gsl_matrix_get;
    const auto matrix_complex_set              = gsl_matrix_complex_set;
    const auto matrix_complex_get              = gsl_matrix_complex_get;
    const auto blas_cgemv                      = gsl_blas_zgemv;
    const auto blas_ctrsv                      = gsl_blas_ztrsv;
    const auto blas_sgemv                      = gsl_blas_dgemv;
    const auto blas_sgemm                      = gsl_blas_dgemm;
    const auto blas_cgemm                      = gsl_blas_zgemm;
};

namespace gsl_single_precision {
    using      gsl_real                        = float;
    using      complex                         = gsl_complex_float;
    using      vector                          = gsl_vector_float;
    using      vector_view                     = gsl_vector_float_view;
    using      vector_const_view               = gsl_vector_float_const_view;
    using      vector_complex                  = gsl_vector_complex_float;
    using      vector_complex_view             = gsl_vector_complex_float_view;
    using      vector_complex_const_view       = gsl_vector_complex_float_const_view;
    using      matrix                          = gsl_matrix_float;
    using      matrix_view                     = gsl_matrix_float_view;
    using      matrix_const_view               = gsl_matrix_float_const_view;
    using      matrix_complex                  = gsl_matrix_complex_float;
    using      matrix_complex_view             = gsl_matrix_complex_float_view;
    using      matrix_complex_const_view       = gsl_matrix_complex_float_const_view;
    using      permutation                     = gsl_permutation;
    const auto sf_gegenpoly_n                  = gsl_sf_gegenpoly_n;
    const auto vector_alloc                    = gsl_vector_float_alloc;
    const auto vector_free                     = gsl_vector_float_free;
    const auto vector_view_array               = gsl_vector_float_view_array;
    const auto vector_const_view_array         = gsl_vector_float_const_view_array;
    const auto vector_const_subvector          = gsl_vector_float_const_subvector;
    const auto vector_complex_alloc            = gsl_vector_complex_float_alloc;
    const auto vector_complex_free             = gsl_vector_complex_float_free;
    const auto vector_complex_view_array       = gsl_vector_complex_float_view_array;
    const auto vector_complex_const_view_array = gsl_vector_complex_float_const_view_array;
    const auto vector_complex_const_subvector  = gsl_vector_complex_float_const_subvector;
    const auto matrix_alloc                    = gsl_matrix_float_alloc;
    const auto matrix_free                     = gsl_matrix_float_free;
    const auto matrix_view_array               = gsl_matrix_float_view_array;
    const auto matrix_const_view_array         = gsl_matrix_float_const_view_array;
    const auto matrix_const_submatrix          = gsl_matrix_float_const_submatrix;
    const auto matrix_transpose                = gsl_matrix_float_transpose;
    const auto matrix_complex_alloc            = gsl_matrix_complex_float_alloc;
    const auto matrix_complex_free             = gsl_matrix_complex_float_free;
    const auto matrix_complex_view_array       = gsl_matrix_complex_float_view_array;
    const auto matrix_complex_const_view_array = gsl_matrix_complex_float_const_view_array;
    const auto matrix_complex_const_submatrix  = gsl_matrix_complex_float_const_submatrix;
    const auto matrix_complex_transpose        = gsl_matrix_complex_float_transpose;
    const auto permutation_alloc               = gsl_permutation_alloc;
    const auto permutation_free                = gsl_permutation_free;
    const auto permute_vector                  = gsl_permute_vector_float;
    const auto permute_vector_complex          = gsl_permute_vector_complex_float;
    const auto linalg_complex_LU_decomp        = gsl_linalg_complex_LU_decomp;
    const auto linalg_complex_LU_svx           = gsl_linalg_complex_LU_svx;
    const auto vector_set                      = gsl_vector_float_set;
    const auto vector_get                      = gsl_vector_float_get;
    const auto vector_complex_set              = gsl_vector_complex_float_set;
    const auto vector_complex_get              = gsl_vector_complex_float_get;
    const auto matrix_set                      = gsl_matrix_float_set;
    const auto matrix_get                      = gsl_matrix_float_get;
    const auto matrix_complex_set              = gsl_matrix_complex_float_set;
    const auto matrix_complex_get              = gsl_matrix_complex_float_get;
    const auto blas_cgemv                      = gsl_blas_cgemv;
    const auto blas_ctrsv                      = gsl_blas_ctrsv;
    const auto blas_sgemv                      = gsl_blas_sgemv;
    const auto blas_sgemm                      = gsl_blas_sgemm;
    const auto blas_cgemm                      = gsl_blas_cgemm;
};



#endif // #ifdef SUPPORT_GSL
#endif // #ifndef __GAMER_GSL_H__
