#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif



// GAMER_GetSlice functions
void ReadOption( int argc, char **argv );
void Init_MPI( int *argc, char ***argv );
void CheckParameter();
void TakeNote( int argc, char *argv[] );
void Refine2TargetLevel();
void PreparePatch( const int lv, const int PID, const int Buffer, real FData[], real CData[], const int BC_Face[] );
void AllocateOutputArray();
void StoreData( const int lv, const int PID, real FData[], const int Buffer, real *Out );
void Output();
#ifndef SERIAL
void SumOverRanks();
#endif
void Init_TargetDomain();
void Init_Convert2Temp();
void GetCandidateBox();
bool WithinCandidateBox( const int *Corner, const int Size, const int Buf );
void LoadTree();
void GetDivVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh );
void GetCurlVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh );
void GetPres( real FData[], const int FSize, const int Buffer, const int NextIdx );
void GetTemp( real FData[], const int FSize, const int Buffer, const int NextIdx, const real Con2T );
#if ( MODEL == ELBDM )
void GetELBDM_Vel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh,
                   const bool Com2Phy );
#endif
#ifdef SUPPORT_HDF5
void SetHDF5Info( hid_t &H5_FileID );
#endif


// GAMER functions
void LoadData();
void Init_MemAllocate();
void Init_RecordBasePatch();
void End_MemFree();
#ifndef SERIAL
void Buf_AllocateBufferPatch( const int lv );
void Buf_AllocateBufferPatch_Base();
void Buf_RecordBoundaryPatch( const int lv );
void Buf_RecordBoundaryPatch_Base();
void Buf_GetBufferData( const int lv, const int Package, const int ParaBuffer );
void Buf_RecordExchangeDataPatchID( const int lv );
void MPI_Exit();
void MPI_ExchangeBufferPosition( int NSend[26], int NRecv[26], int *Send_PosList[26], int *Recv_PosList[26] );
void MPI_ExchangeInfo( const int TargetRank[2], const int SendSize[2], const int RecvSize[2],
                       real *SendBuffer[2], real *RecvBuffer[2] );
#endif
void Buf_SortBoundaryPatch( const int NPatch, int *IDList, int *PosList );
void FindFather( const int lv );
void SiblingSearch( const int lv );
void SiblingSearch_Base();
void Flu_Restrict( const int lv, const bool GetAvePot, const bool GetAveParDens );
int  TABLE_01( const int SibIndex, const char dim, const int w0, const int w1, const int w2 );
int  TABLE_02( const int LocalID, const char dim, const int w0, const int w1 );
int  TABLE_03( const int SibID, const int Count );
int  TABLE_04( const int SibID );
int  TABLE_05( const int SibID );
int  TABLE_07( const int SibID, const int Count );
void Aux_Message( FILE *Type, const char *Format, ... );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
void Hydro_BoundaryCondition_Outflow( real *Array, const int BC_Face, const int NVar_Flu, const int GhostSize, 
                                      const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ, 
                                      const int Idx_Start[], const int Idx_End[] );
void Int_Table( const IntScheme_t IntScheme, int &NSide, int &NGhost );
void Interpolate( real CData [], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData [], const int FSize[3], const int FStart[3], 
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase, const bool Monotonic,
                  const real MonoCoeff );
#if ( MODEL == ELBDM )
real ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped );
#endif



#endif // __PROTOTYPE_H__
