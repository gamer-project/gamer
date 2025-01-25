#ifndef __TESTPROB_H__
#define __TESTPROB_H__



// common headers
#include "ReadPara.h"


// common function prototypes
static void Validate();
static void SetParameter();
static void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                       const int lv, double AuxArray[] );
#ifdef MHD
static void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                         const int lv, double AuxArray[] );
#endif


// function pointers of various user-specified routines
extern void (*Init_Function_User_Ptr)( real fluid[], const double x, const double y, const double z, const double Time,
                                       const int lv, double AuxArray[] );
#ifdef MHD
extern void (*Init_Function_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                              const int lv, double AuxArray[] );
extern double (*Init_BField_ByVecPot_User_Ptr)( const double x, const double y, const double z, const double Time,
                                                const int lv, const char Component, double AuxArray[] );
#endif
extern void (*Init_ByFile_User_Ptr)( real fluid_out[], const real fluid_in[], const int nvar_in,
                                     const double x, const double y, const double z, const double Time,
                                     const int lv, double AuxArray[] );
extern void (*Init_Field_User_Ptr)();
extern void (*Init_User_Ptr)();
extern void (*Init_User_AfterPoisson_Ptr)();
extern void (*Output_User_Ptr)();
extern void (*Output_UserWorkBeforeOutput_Ptr)();
extern bool (*Flag_Region_Ptr)( const int i, const int j, const int k, const int lv, const int PID );
extern bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
extern void (*Flag_UserWorkBeforeFlag_Ptr)( const double Time, const int lv );
extern double (*Mis_GetTimeStep_User_Ptr)( const int lv, const double dTime_dt );
extern void (*Mis_UserWorkBeforeNextLevel_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt );
extern void (*Mis_UserWorkBeforeNextSubstep_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt );
extern void (*Aux_Record_User_Ptr)();
extern void (*BC_User_Ptr)( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                            const int GhostSize, const int idx[], const double pos[], const double Time,
                            const int lv, const int TFluVarIdxList[], double AuxArray[] );
#ifdef MHD
extern void (*BC_BField_User_Ptr)( real magnetic[], const double x, const double y, const double z, const double Time,
                                   const int lv, double AuxArray[] );
#endif
extern int  (*Flu_ResetByUser_Func_Ptr)( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                         const double dt, const int lv, double AuxArray[] );
extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );
#ifdef MHD
extern double (*MHD_ResetByUser_VecPot_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[] );
extern double (*MHD_ResetByUser_BField_Ptr)( const double x, const double y, const double z, const double Time,
                                             const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                             const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                             const int i, const int j, const int k );
#endif
extern void (*End_User_Ptr)();
#ifdef GRAVITY
extern real (*Poi_AddExtraMassForGravity_Ptr)( const double x, const double y, const double z, const double Time,
                                               const int lv, double AuxArray[] );
extern void (*Poi_UserWorkBeforePoisson_Ptr)( const double Time, const int lv );
extern void (*Init_ExtAcc_Ptr)();
extern void (*Init_ExtPot_Ptr)();
extern void (*End_ExtAcc_Ptr)();
extern void (*End_ExtPot_Ptr)();
#endif
#ifdef PARTICLE
extern void (*Par_Init_ByFunction_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                        real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                        real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                        long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                        long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
extern void (*Par_Init_ByFile_User_Ptr)();
extern void (*Par_Init_Attribute_User_Ptr)();
#endif
#if ( MODEL == HYDRO )
extern void (*EoS_Init_Ptr)();
extern void (*EoS_End_Ptr)();
#endif
extern void (*Src_Init_User_Ptr)();
extern void (*Src_End_User_Ptr)();
extern void (*Src_WorkBeforeMajorFunc_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                                 double AuxArray_Flt[], int AuxArray_Int[] );
extern void (*Init_DerivedField_User_Ptr)();
#ifdef FEEDBACK
extern void (*FB_Init_User_Ptr)();
extern void (*FB_End_User_Ptr)();
extern int  (*FB_User_Ptr)( const int lv, const double TimeNew, const double TimeOld, const double dt,
                            const int NPar, const long *ParSortID, real_par *ParAttFlt[PAR_NATT_FLT_TOTAL], long_par *ParAttInt[PAR_NATT_INT_TOTAL],
                            real (*Fluid)[FB_NXT][FB_NXT][FB_NXT], const double EdgeL[], const double dh, bool CoarseFine[],
                            const int TID, RandomNumber_t *RNG );
#endif



#endif // #ifndef __TESTPROB_H__
