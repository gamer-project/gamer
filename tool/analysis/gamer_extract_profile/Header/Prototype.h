#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// ExtractProfile functions
void TakeNote( int argc, char **argv );
void ReadOption( int argc, char **argv );
void CheckParameter();
void SetMaxRhoPos( const int AveN );
void End();
void Init_ShellAve();
void Output_ShellAve();
void ShellAverage();
void GetRMS();
void GetMaxRho();
double GetMinShellWidth( const double TCen[], const double TCen_Map[] );
void GetR( const int n, double &R, double &dR );

// GAMER functions
void LoadData();
void FindFather( const int lv );
void Init_RecordBasePatch();
int  TABLE_01( const int SibIndex, const char dim, const int w0, const int w1, const int w2 );
int  TABLE_02( const int LocalID, const char dim, const int w0, const int w1 );
int  TABLE_03( const int SibID, const int Count );
int  TABLE_04( const int SibID );
void SiblingSearch( const int lv );
void SiblingSearch_Base();
void Flu_Restrict( const int lv, const bool GetAvePot, const bool GetAveParDens );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void Aux_Message( FILE *Type, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
void Prepare_PatchData( const int lv, real *h_Input_Array, const int GhostSize, const int NPG, const int *PID0_List, 
                        const long TVar, const IntScheme_t IntScheme, const NSide_t NSide, const bool IntPhase );
void Int_Table( const IntScheme_t IntScheme, int &NSide, int &NGhost );
void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3], 
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic );
void Output_Patch( const int lv, const int PID, const char *comment );
#if ( MODEL == ELBDM )
real ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped );
#endif



#endif // __PROTOTYPE_H__
