#ifndef __GAMER_GSL_H__
#define __GAMER_GSL_H__

#ifdef SUPPORT_GSL

#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_permute_matrix.h>
#include <gsl/gsl_linalg.h>

namespace gsl_double_precision {
    using       complex                          =  gsl_complex;
    using       vector                           =  gsl_vector;
    using       vector_complex                   =  gsl_vector_complex;
    using       matrix                           =  gsl_matrix;
    using       matrix_complex                   =  gsl_matrix_complex;
    using       matrix_complex_view              =  gsl_matrix_complex_view;
    using       matrix_complex_const_view        =  gsl_matrix_complex_const_view;
    using       vector_complex_view              =  gsl_vector_complex_view;
    using       vector_complex_const_view        =  gsl_vector_complex_const_view;
    using       permutation                      =  gsl_permutation;
    const auto  sf_gegenpoly_n                   =  gsl_sf_gegenpoly_n;
    const auto  matrix_complex_view_array        =  gsl_matrix_complex_view_array;
    const auto  matrix_complex_const_view_array  =  gsl_matrix_complex_const_view_array;
    const auto  matrix_complex_const_submatrix   =  gsl_matrix_complex_const_submatrix;
    const auto  vector_complex_view_array        =  gsl_vector_complex_view_array;
    const auto  vector_complex_const_view_array  =  gsl_vector_complex_const_view_array;
    const auto  matrix_complex_transpose         =  gsl_matrix_complex_transpose;
    const auto  permutation_alloc                =  gsl_permutation_alloc;
    const auto  permutation_free                 =  gsl_permutation_free;
    const auto  permute_vector                   =  gsl_permute_vector;
    const auto  permute_vector_complex           =  gsl_permute_vector_complex;
    const auto  linalg_complex_LU_decomp         =  gsl_linalg_complex_LU_decomp;
    const auto  linalg_complex_LU_svx            =  gsl_linalg_complex_LU_svx;
    const auto  vector_complex_set               =  gsl_vector_complex_set;
    const auto  vector_complex_get               =  gsl_vector_complex_get;
    const auto  matrix_complex_set               =  gsl_matrix_complex_set;
    const auto  matrix_complex_get               =  gsl_matrix_complex_get;
    const auto  complex_abs                      =  gsl_complex_abs;
    const auto  complex_mul                      =  gsl_complex_mul;
    const auto  complex_div                      =  gsl_complex_div;
    const auto  complex_add                      =  gsl_complex_add;
    const auto  complex_sub                      =  gsl_complex_sub;
    const auto  blas_cgemv                       =  gsl_blas_zgemv;
};

namespace gsl_single_precision {
    using      complex                          =  gsl_complex_float;
    using      vector                           =  gsl_vector_float;
    using      vector_complex                   =  gsl_vector_complex_float;
    using      matrix                           =  gsl_matrix_float;
    using      matrix_complex                   =  gsl_matrix_complex_float;
    using      matrix_complex_view              =  gsl_matrix_complex_float_view;
    using      matrix_complex_const_view        =  gsl_matrix_complex_float_const_view;
    using      vector_complex_view              =  gsl_vector_complex_float_view;
    using      vector_complex_const_view        =  gsl_vector_complex_float_const_view;
    using      permutation                      =  gsl_permutation;
    const auto matrix_complex_view_array        =  gsl_matrix_complex_float_view_array;
    const auto matrix_complex_const_view_array  =  gsl_matrix_complex_float_const_view_array;
    const auto matrix_complex_const_submatrix   =  gsl_matrix_complex_float_const_submatrix;
    const auto vector_complex_view_array        =  gsl_vector_complex_float_view_array;
    const auto vector_complex_const_view_array  =  gsl_vector_complex_float_const_view_array;
    const auto matrix_complex_transpose         =  gsl_matrix_complex_float_transpose;
    const auto permutation_alloc                =  gsl_permutation_alloc;
    const auto permutation_free                 =  gsl_permutation_free;
    const auto permute_vector                   =  gsl_permute_vector_float;
    const auto permute_vector_complex           =  gsl_permute_vector_complex_float;
    const auto linalg_complex_LU_decomp         =  gsl_linalg_complex_LU_decomp;
    const auto linalg_complex_LU_svx            =  gsl_linalg_complex_LU_svx;
    const auto vector_complex_set               =  gsl_vector_complex_float_set;
    const auto vector_complex_get               =  gsl_vector_complex_float_get;
    const auto matrix_complex_set               =  gsl_matrix_complex_float_set;
    const auto matrix_complex_get               =  gsl_matrix_complex_float_get;
    const auto complex_abs                      =  gsl_complex_float_abs;
    const auto complex_mul                      =  gsl_complex_float_mul;
    const auto complex_div                      =  gsl_complex_float_div;
    const auto complex_add                      =  gsl_complex_float_add;
    const auto complex_sub                      =  gsl_complex_float_sub;
    const auto sf_gegenpoly_n                   =  gsl_sf_gegenpoly_n;
    const auto blas_cgemv                       =  gsl_blas_cgemv;
};

#endif // #ifdef SUPPORT_GSL
#endif // #ifndef __GAMER_GSL_H__