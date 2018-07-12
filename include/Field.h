#ifndef __FIELD_H__
#define __FIELD_H__



#include "Macro.h"

typedef int FieldIdx_t;
const FieldIdx_t Idx_Undefined = -1;   // must be outside the range [0, NCOMP_TOTAL-1]


// ***********************************************************
// ** Add new field indices below with the format:          **
// **                                                       **
// ** SET_GLOBAL( FieldIdx_t Idx_NewField, Idx_Undefined ); **
// **                                                       **
// ***********************************************************
#if   ( MODEL == HYDRO  ||  MODEL == MHD )
SET_GLOBAL( FieldIdx_t Idx_Dens,    Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomX,    Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomY,    Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomZ,    Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Engy,    Idx_Undefined );


#elif ( MODEL == ELBDM )


#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// field labels
SET_GLOBAL( const char *FieldLabel[NCOMP_TOTAL] );
#ifdef GRAVITY
SET_GLOBAL( const char *PotLabel, "Pote" );  // potential label is currently fixed
#endif



#endif  // #ifndef __FIELD_H__
