#ifndef __CUPOT_H__
#define __CUPOT_H__



#ifdef GRAVITY



// *****************************************************************
// ** This header will be included by all CPU/GPU gravity solvers **
// *****************************************************************


// include "Macro.h" and "Typedef.h" here since the header "GAMER.h" is NOT included in GPU solvers
#ifdef __CUDACC__
# include "Macro.h"
# include "Typedef.h"
#else
# include "GAMER.h"
#endif


// allow GPU to output messages in the debug mode
#ifdef GAMER_DEBUG
#  include "stdio.h"
#endif


// experiments show that the following macros give lower performance even in Fermi GPUs
/*
// faster integer multiplication in GPU
#if ( defined __CUDACC__  &&  __CUDA_ARCH__ >= 200 )
   #define __umul24( a, b )   ( (a)*(b) )
   #define  __mul24( a, b )   ( (a)*(b) )
#endif
*/



// ####################
// ## macros for SOR ##
// ####################
#if   ( POT_SCHEME == SOR )

// determine the version of the GPU Poisson solver (10to14cube vs. 16to18cube)
#define USE_PSOLVER_10TO14


// blockDim.z for the GPU Poisson solver
#if   ( POT_GHOST_SIZE == 1 )
#        define POT_BLOCK_SIZE_Z      4
#elif ( POT_GHOST_SIZE == 2 )
#        define POT_BLOCK_SIZE_Z      2
#elif ( POT_GHOST_SIZE == 3 )
#        define POT_BLOCK_SIZE_Z      2
#elif ( POT_GHOST_SIZE == 4 )
#  ifdef USE_PSOLVER_10TO14
#        define POT_BLOCK_SIZE_Z      5
#  else
#        define POT_BLOCK_SIZE_Z      1
#  endif
#elif ( POT_GHOST_SIZE == 5 )
#  ifdef USE_PSOLVER_10TO14
#     if   ( GPU_ARCH == FERMI )
#        ifdef FLOAT8
#        define POT_BLOCK_SIZE_Z      2
#        else
#        define POT_BLOCK_SIZE_Z      4
#        endif
#     elif ( GPU_ARCH == KEPLER )
#        ifdef FLOAT8
#        define POT_BLOCK_SIZE_Z      2
#        else
#        define POT_BLOCK_SIZE_Z      8
#        endif
#     elif ( GPU_ARCH == MAXWELL )
#        ifdef FLOAT8
#        define POT_BLOCK_SIZE_Z      2      // not optimized yet
#        else
#        define POT_BLOCK_SIZE_Z      4      // not optimized yet
#        endif
#     elif ( GPU_ARCH == PASCAL )
#        ifdef FLOAT8
#        define POT_BLOCK_SIZE_Z      2      // not optimized yet
#        else
#        define POT_BLOCK_SIZE_Z      4      // not optimized yet
#        endif
#     elif ( GPU_ARCH == VOLTA )
#        ifdef FLOAT8
#        define POT_BLOCK_SIZE_Z      2      // not optimized yet
#        else
#        define POT_BLOCK_SIZE_Z      4      // not optimized yet
#        endif
#     else
#        define POT_BLOCK_SIZE_Z      NULL_INT
#        ifdef GPU
#        error : UNKNOWN GPU_ARCH !!
#        endif
#     endif // GPU_ARCH
#  else
#        define POT_BLOCK_SIZE_Z      1
#  endif // #ifdef USE_PSOLVER_10TO14 ... else ...
#endif // POT_GHOST_SIZE


// optimization options for CUPOT_PoissonSolver_SOR_10to14cube.cu
#ifdef USE_PSOLVER_10TO14

// load density into shared memory for higher performance
#  ifndef FLOAT8
#     define SOR_RHO_SHARED
#  endif

#  if ( POT_GHOST_SIZE == 5 )
// use shuffle reduction
// --> only work for POT_GHOST_SIZE == 5 since # threads must be a multiple of warpSize
// --> although strickly speaking the shuffle functions do NOT work for double precision, but experiments
//     show that residual_sum += (float)residual, where residual_sum is double, gives acceptable accuracy
#  if ( GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA )
#     define SOR_USE_SHUFFLE
#  endif

// use padding to reduce shared memory bank conflict (optimized for POT_GHOST_SIZE == 5 only)
// --> does NOT work for FLOAT8 due to the lack of shared memory
// --> does NOT work with FERMI GPUs because SOR_USE_PADDING requires POT_BLOCK_SIZE_Z == 8 but FERMI does NOT support that
#  ifndef FLOAT8
#     define SOR_USE_PADDING
#  endif
#  endif // #if ( POT_GHOST_SIZE == 5 )

// frequency of reduction
#  define SOR_MOD_REDUCTION 2

// load coarse-grid potential into shared memory for higher performance
#  if ( !defined FLOAT8  &&  !defined SOR_USE_PADDING  &&  GPU_ARCH != KEPLER )
#     define SOR_CPOT_SHARED
#  endif

#endif // #ifdef USE_PSOLVER_10TO14



// ###################
// ## macros for MG ##
// ###################
#elif ( POT_SCHEME == MG  )

// blockDim.x for the GPU Poisson solver
#if   ( POT_GHOST_SIZE == 1 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 2 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 3 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 4 )
      #define POT_BLOCK_SIZE_X      256
#elif ( POT_GHOST_SIZE == 5 )
      #define POT_BLOCK_SIZE_X      256
#endif

#endif // POT_SCHEME


// blockDim.z for the GPU Gravity solver
#define GRA_BLOCK_SIZE              256


// dt solver for gravity
#define DT_GRA_BLOCK_SIZE           256

// use shuffle reduction in the KEPLER and later GPUs
#if ( GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA )
#   define DT_GRA_USE_SHUFFLE
#endif


// maximum size of the arrays ExtPot_AuxArray and ExtAcc_AuxArray
#define EXT_POT_NAUX_MAX            10
#define EXT_ACC_NAUX_MAX            10


// warp size (which must be the same as the CUDA predefined constant "warpSize")
// --> please refer to https://en.wikipedia.org/wiki/CUDA#Version_features_and_specifications
//     for information on warp size
#ifdef __CUDACC__
#if ( GPU_ARCH == FERMI  ||  GPU_ARCH == KEPLER  ||  GPU_ARCH == MAXWELL  ||  GPU_ARCH == PASCAL  ||  GPU_ARCH == VOLTA )
// CUFLU.h will define WARP_SIZE as well
#  ifndef WARP_SIZE
#  define WARP_SIZE 32
#  endif
#elif defined GPU
#  error : UNKNOWN GPU_ARCH !!
#endif
#endif // #ifdef __CUDACC__



// #########################
// ## CPU/GPU integration ##
// #########################

// GPU device function specifier
#ifdef __CUDACC__
# define GPU_DEVICE __forceinline__ __device__
#else
# define GPU_DEVICE
#endif

// unified CPU/GPU loop
#ifdef __CUDACC__
# define CGPU_LOOP( var, niter )    for (int (var)=threadIdx.x; (var)<(niter); (var)+=blockDim.x)
#else
# define CGPU_LOOP( var, niter )    for (int (var)=0;           (var)<(niter); (var)++          )
#endif



#endif // #ifdef GRAVITY



#endif // #ifndef __CUPOT_H__
