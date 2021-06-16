#ifndef __GLOBAL_H__
#define __GLOBAL_H__



extern AMR_t    patch1, patch2;
extern char    *FileName_In1, *FileName_In2, *FileName_Out;
extern bool     WithPot1, WithPot2, WithMagCC1, WithMagCC2, WithMagFC1, WithMagFC2;
extern bool     UseCorner, WithPar1, WithPar2;
extern int      WithParDens1, WithParDens2, NParVarOut1, NParVarOut2;
extern long     NPar1, NPar2;
extern double   TolErr;
extern real   **ParData1, **ParData2;



#endif // #ifndef __GLOBAL_H__
