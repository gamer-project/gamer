#ifndef __GLOBAL_H__
#define __GLOBAL_H__



extern AMR_t       amr;
extern char       *FileName_In;
extern tree_t     *tree;
extern int        *BaseP;
extern int         NShell, NIn, NOut, NX0_TOT[3], DumpID, OutputParDens;
extern long        Step;
extern bool        OutputPot, Periodic, InputScale;
extern double      ShellWidth, GetNShell, LogBin, Center[3], Center_Map[3], MaxRadius, UseMaxRhoPos_R, GAMMA, INT_MONO_COEFF;
extern double      Time[NLEVEL];
extern bool        UseTree;
extern bool        NeedGhost;
extern bool        GetAvePot;
extern double      NewtonG;

#if ( MODEL == ELBDM )
extern bool        ELBDM_GetVir;
extern double      ELBDM_ETA;
#endif



#endif // #ifndef __GLOBAL_H__
