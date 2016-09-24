#ifndef __GLOBAL_H__
#define __GLOBAL_H__



extern AMR_t    amr;
extern char    *FileName_In;
extern double   Time[NLEVEL];
extern int      DumpID;

#if   ( MODEL == HYDRO )
extern double   GAMMA;

#elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
extern double   ELBDM_ETA;

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL

extern int  SizeScale  [3];
extern int  CornerScale[3];
extern bool OutputPot;



#endif // #ifndef __GLOBAL_H__
