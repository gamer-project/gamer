#ifndef __CUDA_CONSTMEMORY_H__
#define __CUDA_CONSTMEMORY_H__



#include "Macro.h"
#include "Typedef.h"
#include "CUPOT.h"



#if ( NCOMP_PASSIVE > 0 )
SET_GLOBAL( __constant__ int  c_NormIdx[NCOMP_PASSIVE] );
#else
SET_GLOBAL( __constant__ int *c_NormIdx, NULL );
#endif

#ifdef GRAVITY
SET_GLOBAL( __constant__ double c_ExtAcc_AuxArray[EXT_ACC_NAUX_MAX] );
SET_GLOBAL( __constant__ double c_ExtPot_AuxArray[EXT_POT_NAUX_MAX] );

SET_GLOBAL( __constant__ real c_Mp[3] );
SET_GLOBAL( __constant__ real c_Mm[3] );
#endif



#endif // #ifndef __CUDA_CONSTMEMORY_H__
