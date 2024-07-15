#ifndef __CUDA_CONSTMEMORY_H__
#define __CUDA_CONSTMEMORY_H__



#include "Macro.h"
#include "Typedef.h"



#if ( MODEL == HYDRO )
SET_GLOBAL( __constant__ double c_EoS_AuxArray_Flt[EOS_NAUX_MAX  ] );
SET_GLOBAL( __constant__ int    c_EoS_AuxArray_Int[EOS_NAUX_MAX  ] );
SET_GLOBAL( __constant__ real*  c_EoS_Table       [EOS_NTABLE_MAX] );
#endif

#if ( NCOMP_PASSIVE > 0 )
SET_GLOBAL( __constant__ int  c_NormIdx[NCOMP_PASSIVE] );
SET_GLOBAL( __constant__ int  c_FracIdx[NCOMP_PASSIVE] );
#else
SET_GLOBAL( __constant__ int *c_NormIdx, NULL );
SET_GLOBAL( __constant__ int *c_FracIdx, NULL );
#endif

#ifdef GRAVITY
SET_GLOBAL( __constant__ double c_ExtAcc_AuxArray    [EXT_ACC_NAUX_MAX] );
SET_GLOBAL( __constant__ double c_ExtPot_AuxArray_Flt[EXT_POT_NAUX_MAX] );
SET_GLOBAL( __constant__ int    c_ExtPot_AuxArray_Int[EXT_POT_NAUX_MAX] );

SET_GLOBAL( __constant__ real c_Mp[3] );
SET_GLOBAL( __constant__ real c_Mm[3] );
#endif

#if ( MODEL == HYDRO )
SET_GLOBAL( __constant__ double c_Src_Dlep_AuxArray_Flt[SRC_NAUX_DLEP] );
SET_GLOBAL( __constant__ int    c_Src_Dlep_AuxArray_Int[SRC_NAUX_DLEP] );
SET_GLOBAL( __constant__ double c_Src_EC_AuxArray_Flt[SRC_NAUX_EC] );
SET_GLOBAL( __constant__ int    c_Src_EC_AuxArray_Int[SRC_NAUX_EC] );
#endif
SET_GLOBAL( __constant__ double c_Src_User_AuxArray_Flt[SRC_NAUX_USER] );
SET_GLOBAL( __constant__ int    c_Src_User_AuxArray_Int[SRC_NAUX_USER] );



#endif // #ifndef __CUDA_CONSTMEMORY_H__
