#ifndef __CUDA_CHECK_ERROR_H__
#define __CUDA_CHECK_ERROR_H__



#include "Macro.h"


void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );


// CUDA error check
#define CUDA_CHECK_ERROR( Call )   CUDA_Check_Error( Call, __FILE__, __LINE__, __FUNCTION__ )

inline void CUDA_Check_Error( cudaError Return, const char *File, const int Line, const char *Func )
{
   if ( Return != cudaSuccess )
      Aux_Error( ERROR_INFO, "CUDA ERROR : %s !!\n", cudaGetErrorString(Return) );
}



#endif // #ifndef __CUDA_CHECK_ERROR_H__
