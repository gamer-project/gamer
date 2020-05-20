#ifndef __GLOBAL_H__
#define __GLOBAL_H__



extern int        CanMin[3][3], CanMax[3][3], CanBuf;
extern AMR_t      amr;
extern ParaVar_t  ParaVar;
extern double     GAMMA, MU, H0;
extern int        NX0_TOT[3], NX0[3], MyRank, MyRank_X[3], SibRank[26], NGPU, NGPU_X[3], TargetLevel, DumpID;
extern int        *BaseP, *BounP_IDMap[NLEVEL][26], SendP_NList[NLEVEL][26], *SendP_IDList[NLEVEL][26];
extern int        RecvP_NList[NLEVEL][26], *RecvP_IDList[NLEVEL][26], NPatchComma[NLEVEL][28];
extern int        TargetX[3], X[3], NLoad, NOut, ShiftScale[3], ExtBC;
extern double     Time[NLEVEL];
extern long int   Step;
extern bool       WithUnit, Comoving;
extern double     Unit_L, Unit_M, Unit_T, Unit_V, Unit_D, Unit_E, Unit_P, Convert2Temp;
extern bool       OutputPot, Shift2Center, OutputTemp;
extern int        BufSize, OutputParDens;
extern tree_t    *tree;
extern bool       UseTree;
extern char      *FileName_In;

#if ( MODEL == ELBDM )
extern double     ELBDM_ETA;
extern double     ELBDM_MASS;
extern bool       ELBDM_IntPhase;
#endif



#endif // #ifndef __GLOBAL_H__
