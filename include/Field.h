#ifndef __FIELD_H__
#define __FIELD_H__



#include "Macro.h"
#include "Typedef.h"

const FieldIdx_t Idx_Undefined = -1;   // must be outside the range [0, NCOMP_TOTAL-1]


// ***********************************************************
// ** Add new field indices below with the format:          **
// **                                                       **
// ** SET_GLOBAL( FieldIdx_t Idx_NewField, Idx_Undefined ); **
// **                                                       **
// ***********************************************************

#if   ( MODEL == HYDRO )
// main fields
SET_GLOBAL( FieldIdx_t Idx_Dens,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomX,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomY,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_MomZ,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Engy,          Idx_Undefined );
#ifdef COSMIC_RAY
SET_GLOBAL( FieldIdx_t Idx_CRay,          Idx_Undefined );
#endif
#ifdef DUAL_ENERGY
SET_GLOBAL( FieldIdx_t Idx_Dual,          Idx_Undefined );
#endif

// Grackle fields
SET_GLOBAL( FieldIdx_t Idx_e,             Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HI,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HII,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeI,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeII,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HeIII,         Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HM,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_H2I,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_H2II,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_DI,            Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_DII,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_HDI,           Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Metal,         Idx_Undefined );


#elif ( MODEL == ELBDM )

#if   ( ELBDM_SCHEME == ELBDM_WAVE )
SET_GLOBAL( FieldIdx_t Idx_Dens,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Real,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Imag,          Idx_Undefined );
#elif ( ELBDM_SCHEME == ELBDM_HYBRID )
SET_GLOBAL( FieldIdx_t Idx_Dens,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Real,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Imag,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Phas,          Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_Stub,          Idx_Undefined );
#else
#  error : ERROR : unsupported ELBDM_SCHEME !!
#endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL



// ***************************************************************
// ** Add new particle attribute indices below with the format: **
// **                                                           **
// ** SET_GLOBAL( FieldIdx_t Idx_NewAttribute, Idx_Undefined ); **
// **                                                           **
// ***************************************************************

#ifdef PARTICLE
SET_GLOBAL( FieldIdx_t Idx_ParMass,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParPosZ,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParVelZ,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParTime,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParType,       Idx_Undefined );
#ifdef STORE_PAR_ACC
SET_GLOBAL( FieldIdx_t Idx_ParAccX,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParAccY,       Idx_Undefined );
SET_GLOBAL( FieldIdx_t Idx_ParAccZ,       Idx_Undefined );
#endif
#ifdef STAR_FORMATION
SET_GLOBAL( FieldIdx_t Idx_ParCreTime,    Idx_Undefined );
#endif
SET_GLOBAL( FieldIdx_t Idx_ParMetalFrac,  Idx_Undefined );
#endif // #ifdef PARTICLE



// field and particle attribute labels
SET_GLOBAL( char FieldLabel[NCOMP_TOTAL][MAX_STRING] );
#ifdef GRAVITY
SET_GLOBAL( const char *PotLabel, "Pote" );  // potential label is currently fixed
#endif
#ifdef MHD
SET_GLOBAL( char MagLabel[NCOMP_MAG][MAX_STRING] );
#endif
#ifdef PARTICLE
SET_GLOBAL( char ParAttFltLabel[PAR_NATT_FLT_TOTAL][MAX_STRING] );
SET_GLOBAL( char ParAttIntLabel[PAR_NATT_INT_TOTAL][MAX_STRING] );
#endif



#endif  // #ifndef __FIELD_H__
