#ifndef __CUPOT_H__
#define __CUPOT_H__



#ifdef GRAVITY



// *****************************************************************
// ** This header will be included by all CPU/GPU gravity solvers **
// *****************************************************************


// include "Typedef" here since the header "GAMER.h" is NOT included in GPU solvers
#include "Typedef.h"


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
#     else
#        define POT_BLOCK_SIZE_Z      NULL_INT
#        ifdef GPU
#        error : UNKNOWN GPU_ARCH !!
#        endif
#     endif
#  else
#        define POT_BLOCK_SIZE_Z      1
#  endif // #ifdef USE_PSOLVER_10TO14 ... else ...
#endif // POT_GHOST_SIZE


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
#define GRA_BLOCK_SIZE_Z            4


// maximum size of the arrays ExtPot_AuxArray and ExtAcc_AuxArray
#define EXT_POT_NAUX_MAX            10
#define EXT_ACC_NAUX_MAX            10



#endif // #ifdef GRAVITY



#endif // #ifndef __CUPOT_H__
