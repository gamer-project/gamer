#ifndef __CUDA_CHECK_ERROR_H__
#define __CUDA_CHECK_ERROR_H__



#include "Macro.h"


void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );


// CUDA error check
#define CUDA_CHECK_ERROR( Call )   CUDA_Check_Error( Call, __FILE__, __LINE__, __FUNCTION__ )

inline void CUDA_Check_Error( cudaError Return, const char *File, const int Line, const char *Func )
{
   if ( Return != cudaSuccess )
      Aux_Error( File, Line, Func, "CUDA ERROR : %s !!\n", cudaGetErrorString(Return) );
}



// in CUDA_CHECK_MALLOC(), we must use "Call; cudaError_t Return = cudaGetLastError();" instead of "cudaError_t Return = Call;"
// since cudaGetLastError() will reset the last error to cudaSuccess
// --> otherwise CUDA_CHECK_ERROR( cudaGetLastError() ) in, for example, CUAPI_Asyn_FluidSolver(),
//     will fail since the last error has not been reset!
#define CUDA_CHECK_MALLOC( Call )                                                                        \
{                                                                                                        \
   Call;                                                                                                 \
   const cudaError_t Return = cudaGetLastError();                                                        \
   if      ( Return == cudaErrorMemoryAllocation )                                                       \
      return GAMER_FAILED;                                                                               \
   else if ( Return != cudaSuccess )                                                                     \
      Aux_Error( ERROR_INFO, "CUDA ERROR in memory allocation : %s !!\n", cudaGetErrorString(Return) );  \
}



#endif // #ifndef __CUDA_CHECK_ERROR_H__
