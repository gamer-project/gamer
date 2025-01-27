#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



#include "Macro.h"
#include "Typedef.h"
#include "AMR.h"
#include "GatherTree.h"

// Auxiliary
void Aux_Check_MemFree( const double MinMemFree_Total, const char *comment );
void Aux_Check_Conservation( const char *comment );
void Aux_Check_NormalizePassive( const int lv, const char *comment );
void Aux_Check();
void Aux_Check_Finite( const int lv, const char *comment );
void Aux_Check_FluxAllocate( const int lv, const char *comment );
void Aux_Check_Parameter();
void Aux_Check_PatchAllocate( const int lv, const char *comment );
void Aux_Check_ProperNesting( const int lv, const char *comment );
void Aux_Check_Refinement( const int lv, const char *comment );
void Aux_Check_Restrict( const int lv, const char *comment );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
bool Aux_CheckFileExist( const char *FileName );
bool Aux_CheckFolderExist( const char *FolderName );
bool Aux_CheckPermission( const char *FileName, const int perms );
void Aux_GetCPUInfo( const char *FileName );
void Aux_GetMemInfo();
void Aux_Message( FILE *Type, const char *Format, ... );
void Aux_TakeNote();
void Aux_CreateTimer();
void Aux_DeleteTimer();
void Aux_ResetTimer();
void Aux_AccumulatedTiming( const double TotalT, double InitT, double OtherT );
void Aux_Record_Timing();
void Aux_Record_PatchCount();
void Aux_Record_Performance( const double ElapsedTime );
void Aux_Record_CorrUnphy();
void Aux_Record_Center();
int  Aux_CountRow( const char *FileName );
void Aux_ComputeProfile( Profile_t *Prof[], const double Center[], const double r_max_input, const double dr_min,
                         const bool LogBin, const double LogBinRatio, const bool RemoveEmpty, const long TVarBitIdx[],
                         const int NProf, const int MinLv, const int MaxLv, const PatchType_t PatchType,
                         const double PrepTimeIn );
void Aux_FindExtrema( Extrema_t *Extrema, const ExtremaMode_t Mode, const int MinLv, const int MaxLv,
                      const PatchType_t PatchType );
void Aux_FindWeightedAverageCenter( double WeightedAverageCenter[], const double Center_ref[], const double MaxR, const double MinWD,
                                    const long WeightingDensityField, const double TolErrR, const int MaxIter, double *FinaldR, int *FinalNIter );
#ifndef SERIAL
void Aux_Record_BoundaryPatch( const int lv, int *NList, int **IDList, int **PosList );
#endif
void Aux_SwapPointer( void **Ptr1, void **Ptr2 );
template <typename T> void Aux_AllocateArray2D( T** &Array, const int J, const int I );
template <typename T> void Aux_AllocateArray3D( T*** &Array, const int K, const int J, const int I );
template <typename T> void Aux_DeallocateArray2D( T** &Array );
template <typename T> void Aux_DeallocateArray3D( T*** &Array );
template <typename T> int  Aux_LoadTable( T *&Data, const char *FileName, const int NCol_Target, const int TCol[],
                                          const bool RowMajor, const bool AllocMem );
int Aux_IsFinite( const float x );
int Aux_IsFinite( const double x );
void Aux_PauseManually();


// Buffer
#ifndef SERIAL
void Buf_AllocateBufferPatch_Base( AMR_t *Tamr );
void Buf_AllocateBufferPatch( AMR_t *Tamr, const int lv );
void Buf_GetBufferData( const int lv, const int FluSg, const int MagSg, const int PotSg, const GetBufMode_t GetBufMode,
                        const long TVarCC, const long TVarFC, const int ParaBuf, const UseLBFunc_t UseLBFunc );
void Buf_RecordBoundaryFlag( const int lv );
void Buf_RecordBoundaryPatch_Base();
void Buf_RecordBoundaryPatch( const int lv );
void Buf_RecordExchangeDataPatchID( const int lv );
void Buf_RecordExchangeFluxPatchID( const int lv );
void Buf_ResetBufferFlux( const int lv );
void Buf_SortBoundaryPatch( const int NPatch, int *IDList, int *PosList );
#endif // #ifndef SERIAL


// Forward declare classes defined in GatherTree.h
class LB_PatchCount;
class LB_LocalPatchExchangeList;
class LB_GlobalPatchExchangeList;
class LB_GlobalPatch;
class LB_GlobalTree;

// Hydrodynamics
void CPU_FluidSolver( real h_Flu_Array_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                      real h_Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                      real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                      real h_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                      char h_DE_Array_Out[][ CUBE(PS2) ],
                      real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                      real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                      const double h_Corner_Array[][3],
                      const real h_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
                      const bool h_IsCompletelyRefined[],
                      const bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                      gramfe_matmul_float h_GramFE_TimeEvo[][ 2*FLU_NXT ],
                      const int NPatchGroup, const real dt, const real dh,
                      const bool StoreFlux, const bool StoreElectric,
                      const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter,
                      const real ELBDM_Eta, real ELBDM_Taylor3_Coeff, const bool ELBDM_Taylor3_Auto,
                      const double Time, const bool UsePot, const OptExtAcc_t ExtAcc, const MicroPhy_t MicroPhy,
                      const real MinDens, const real MinPres, const real MinEint,
                      const real DualEnergySwitch,
                      const bool NormPassive, const int NNorm, const int NormIdx[],
                      const bool FracPassive, const int NFrac, const int FracIdx[],
                      const bool JeansMinPres, const real JeansMinPres_Coeff,
                      const bool UseWaveFlag );
void Hydro_NormalizePassive( const real GasDens, real Passive[], const int NNorm, const int NormIdx[] );
#if ( MODEL == HYDRO )
real Hydro_Con2Pres( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinPres, const real MinPres, const real Emag,
                     const EoS_DE2P_t EoS_DensEint2Pres, const EoS_GUESS_t EoS_GuessHTilde,
                     const EoS_H2TEM_t EoS_HTilde2Temp, const double EoS_AuxArray_Flt[],
                     const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX], real *EintOut );
real Hydro_Con2Eint( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const bool CheckMinEint, const real MinEint, const real Emag,
                     const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] );
real Hydro_ConEint2Etot( const real Dens, const real MomX, const real MomY, const real MomZ, const real Eint, const real Emag );
real Hydro_Con2Temp( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinTemp, const real MinTemp, const real Emag,
                     const EoS_DE2T_t EoS_DensEint2Temp, const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                     const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] );
real Hydro_Con2Entr( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Passive[], const bool CheckMinEntr, const real MinEntr, const real Emag,
                     const EoS_DE2S_t EoS_DensEint2Entr, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX] );
real Hydro_CheckMinPres( const real InPres, const real MinPres );
real Hydro_CheckMinEint( const real InEint, const real MinEint );
real Hydro_CheckMinTemp( const real InTemp, const real MinTemp );
real Hydro_CheckMinEntr( const real InEntr, const real MinEntr );
real Hydro_CheckMinEintInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real InEngy,
                               const real MinEint, const real Emag );
bool Hydro_IsUnphysical( const IsUnphyMode_t Mode, const real Fields[], const char SingleFieldName[],
                         const real Min, const real Max, const real Emag,
                         const EoS_DE2P_t EoS_DensEint2Pres,
                         const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                         const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                         const real *const EoS_Table[EOS_NTABLE_MAX],
                         const char File[], const int Line, const char Function[], const IsUnphVerb_t Verbose );
#ifdef DUAL_ENERGY
void Hydro_DualEnergyFix( const real Dens, const real MomX, const real MomY, const real MomZ,
                          real &Etot, real &Dual, char &DE_Status, const real Gamma_m1, const real _Gamma_m1,
                          const bool CheckMinPres, const real MinPres, const real DualEnergySwitch,
                          const real Emag );
real Hydro_Con2Dual( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                     const real Emag, const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[],
                     const int EoS_AuxArray_Int[], const real *const EoS_Table[EOS_NTABLE_MAX] );
real Hydro_DensPres2Dual( const real Dens, const real Pres, const real Gamma_m1 );
real Hydro_DensDual2Pres( const real Dens, const real Dual, const real Gamma_m1,
                          const bool CheckMinPres, const real MinPres );
#endif // #ifdef DUAL_ENERGY
#endif // #if ( MODEL == HYDRO )
#ifdef SRHD
real Hydro_Con2HTilde( const real Con[], const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                       const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                       const real *const EoS_Table[EOS_NTABLE_MAX] );
#endif
int Flu_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag, const bool OverlapMPI, const bool Overlap_Sync );
void Flu_AllocateFluxArray( const int lv );
void Flu_Close( const int lv, const int SaveSg_Flu, const int SaveSg_Mag,
                real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                real h_Flu_Array_F_Out[][FLU_NOUT][ CUBE(PS2) ],
                real h_Mag_Array_F_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                char h_DE_Array_F_Out[][ CUBE(PS2) ],
                const int NPG, const int *PID0_List,
                const real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                const real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                const double dt );
void Flu_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                  real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                  real h_Pot_Array_USG_F[][ CUBE(USG_NXT_F) ],
                  double h_Corner_Array_F[][3],
                  bool h_IsCompletelyRefined[],
                  bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                  const int NPG, const int *PID0_List, LB_GlobalTree* GlobalTree );
void Flu_FixUp_Flux( const int lv, const long TVar );
void Flu_FixUp_Restrict( const int FaLv, const int SonFluSg, const int FaFluSg, const int SonMagSg, const int FaMagSg,
                         const int SonPotSg, const int FaPotSg, const long TVarCC, const long TVarFC );
void Flu_BoundaryCondition_User( real *Array, const int NVar_Flu, const int GhostSize, const int ArraySizeX, const int ArraySizeY,
                                 const int ArraySizeZ, const int Idx_Start[], const int Idx_End[],
                                 const int TFluVarIdxList[], const double Time, const double dh, const double *Corner,
                                 const long TVar, const int lv );
void Hydro_BoundaryCondition_Outflow( real *Array, const int BC_Face, const int NVar, const int GhostSize,
                                      const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                      const int Idx_Start[], const int Idx_End[] );
void Flu_CorrAfterAllSync();
#ifndef SERIAL
void Flu_AllocateFluxArray_Buffer( const int lv );
#endif
#if ( MODEL == HYDRO )
void Flu_DerivedField_DivVel( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                              const int NCellInX, const int NCellInY, const int NCellInZ,
                              const int NGhost, const double dh );
void Flu_DerivedField_Mach( real Out[], const real FluIn[], const real MagIn[], const int NFieldOut,
                            const int NCellInX, const int NCellInY, const int NCellInZ,
                            const int NGhost, const double dh );
#endif


// GAMER
void EvolveLevel( const int lv, const double dTime_FaLv );
void InvokeSolver( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld, const double dt_in,
                   const double Poi_Coeff, const int SaveSg_Flu, const int SaveSg_Mag, const int SaveSg_Pot,
                   const bool OverlapMPI, const bool Overlap_Sync );
void Prepare_PatchData( const int lv, const double PrepTime, real *OutputCC, real *OutputFC,
                        const int GhostSize, const int NPG, const int *PID0_List, long TVarCC, long TVarFC,
                        const IntScheme_t IntScheme_CC, const IntScheme_t IntScheme_FC, const PrepUnit_t PrepUnit,
                        const NSide_t NSide, const bool IntPhase, const OptFluBC_t FluBC[], const OptPotBC_t PotBC,
                        const real MinDens, const real MinPres, const real MinTemp, const real MinEntr, const bool DE_Consistency );

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
void Prepare_PatchData_HasWaveCounterpart( const int lv, bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                           const int GhostSize, const int NPG, const int *PID0_List,
                                           const NSide_t NSide, LB_GlobalTree* GlobalTree );
#endif

// Init
void End_GAMER();
void End_MemFree();
void End_MemFree_Fluid();
void End_StopManually( int &Terminate_global );
void Init_BaseLevel();
void Init_GAMER( int *argc, char ***argv );
void Init_Load_DumpTable();
void Init_Load_FlagCriteria();
void Init_Load_Parameter();
void Init_MemoryPool();
void Init_ResetParameter();
void Init_MemAllocate();
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup, const int Pot_NPatchGroup, const int Src_NPatchGroup );
void Init_Parallelization();
void Init_RecordBasePatch();
void Init_Refine( const int lv );
void Init_ByRestart();
void Init_Unit();
void Init_Reload_OldFormat();
void Init_ByFunction();
void Init_TestProb();
void Init_ByFile();
void Init_UniformGrid( const int lv, const bool FindHomePatchForPar );
void Init_Field();
FieldIdx_t AddField( const char *InputLabel, const FixUpFlux_t FixUp_Flux, const FixUpRestrict_t FixUp_Restrict,
                     const NormPassive_t Norm, const IntFracPassive_t IntFrac );
FieldIdx_t GetFieldIndex( const char *InputLabel, const Check_t Check );
#ifdef OPENMP
void Init_OpenMP();
#endif
#ifdef SUPPORT_HDF5
void Init_ByRestart_HDF5( const char *FileName );
#endif
#ifdef SUPPORT_FFTW
void End_FFTW();
void Init_FFTW();
void Patch2Slab( real *VarS, real *SendBuf_Var, real *RecvBuf_Var, long *SendBuf_SIdx, long *RecvBuf_SIdx,
                 int **List_PID, int **List_k, long *List_NSend_Var, long *List_NRecv_Var,
                 const int *List_z_start, const int local_nz, const int FFT_Size[], const int NRecvSlice,
                 const double PrepTime, const long TVar, const bool InPlacePad, const bool ForPoisson, const bool AddExtraMass );
void Slab2Patch( const real *VarS, real *SendBuf, real *RecvBuf, const int SaveSg, const long *List_SIdx,
                 int **List_PID, int **List_k, long *List_NSend, long *List_NRecv, const int local_nz, const int FFT_Size[],
                 const int NSendSlice, const long TVar, const bool InPlacePad );
#endif // #ifdef SUPPORT_FFTW
void Microphysics_Init();
void Microphysics_End();


// Interpolation
void Interpolate( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                  real FData[], const int FSize[3], const int FStart[3],
                  const int NComp, const IntScheme_t IntScheme, const bool UnwrapPhase,
                  const bool Monotonic[], const bool OppSign0thOrder, const bool AllCons,
                  const IntPrim_t IntPrim, const ReduceOrFixMonoCoeff_t ReduceMonoCoeff,
                  const real CMag_IntIter[], const real FMag_IntIter[][NCOMP_MAG] );
void Int_Table( const IntScheme_t IntScheme, int &NSide, int &NGhost );


// Miscellaneous
template <typename T> void  Mis_Idx1D2Idx3D( const int Size[], const T Idx1D, int Idx3D[] );
template <typename U, typename T> U Mis_BinarySearch( const T Array[], U Min, U Max, const T Key );
template <typename U, typename T> U Mis_BinarySearch_Real( const T Array[], U Min, U Max, const T Key );
template <typename T> T     Mis_InterpolateFromTable( const int N, const T Table_x[], const T Table_y[], const T x );
template <typename T> ulong Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
template <typename U, typename T> void  Mis_Heapsort( const U N, T Array[], U IdxTable[] );
template <typename T> int   Mis_Matching_char( const int N, const T Array[], const int M, const T Key[], char Match[] );
template <typename U, typename T> U Mis_Matching_int( const U N, const T Array[], const U M, const T Key[], U Match[] );
template <typename T> bool  Mis_CompareRealValue( const T Input1, const T Input2, const char *comment, const bool Verbose );
template <typename T> void Mis_SortByRows( T const* const* Array, long *IdxTable, const long NSort, const int *SortOrder, const int NOrder );
ulong  Mis_Idx3D2Idx1D( const int Size[], const int Idx3D[] );
double Mis_GetTimeStep( const int lv, const double dTime_SyncFaLv, const double AutoReduceDtCoeff );
double Mis_dTime2dt( const double Time_In, const double dTime_In );
void   Mis_GetTotalPatchNumber( const int lv );
double Mis_Scale2PhySize( const int Scale );
double Mis_Cell2PhySize( const int NCell, const int lv );
int    Mis_Scale2Cell( const int Scale, const int lv );
int    Mis_Cell2Scale( const int NCell, const int lv );
double dt_InvokeSolver( const Solver_t TSolver, const int lv );
void   dt_Prepare_Flu( const int lv, real h_Flu_Array_T[][FLU_NIN_T][ CUBE(PS1) ],
                       real h_Mag_Array_T[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const int NPG, const int *PID0_List );
#ifdef GRAVITY
void   dt_Prepare_Pot( const int lv, real h_Pot_Array_T[][ CUBE(GRA_NXT) ], const int NPG, const int *PID0_List,
                       const double PrepTime );
#endif
void   dt_Close( const real h_dt_Array_T[], const int NPG );
void   CPU_dtSolver( const Solver_t TSolver, real dt_Array[], const real Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                     const real Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const real Pot_Array[][ CUBE(GRA_NXT) ],
                     const double Corner_Array[][3], const int NPatchGroup, const real dh, const real Safety,
                     const MicroPhy_t MicroPhy, const real MinPres, const bool P5_Gradient,
                     const bool UsePot, const OptExtAcc_t ExtAcc, const double TargetTime );


// MPI
#ifndef SERIAL
void Init_MPI( int *argc, char ***argv );
void MPI_ExchangeBoundaryFlag( const int lv );
void MPI_ExchangeBufferPosition( int NSend[26], int NRecv[26], int *Send_PosList[26], int *Recv_PosList[26] );
void MPI_ExchangeData( const int TargetRank[2], const int SendSize[2], const int RecvSize[2],
                       real *SendBuffer[2], real *RecvBuffer[2] );
void MPI_Exit();
template <typename T> void MPI_Alltoallv_GAMER( T * SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, T *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_DataType, MPI_Comm comm );
#endif // #ifndef SERIAL


// Output
void Output_DumpData( const int Stage );
void Output_DumpData_Part( const OptOutputPart_t Part, const bool BaseOnly, const double x, const double y,
                           const double z, const char *FileName );
void Output_DumpData_Total( const char *FileName );
#ifdef SUPPORT_HDF5
void Output_DumpData_Total_HDF5( const char *FileName );
#endif
void Output_DumpManually( int &Dump_global );
void Output_FlagMap( const int lv, const int xyz, const char *comment );
void Output_Flux( const int lv, const int PID, const int Sib, const char *comment );
void Output_PatchCorner( const int lv, const char *comment );
void Output_Patch( const int lv, const int PID, const int FluSg, const int MagSg, const int PotSg, const char *comment );
void Output_PatchMap( const int lv, const int PID, const int TSg, const int Comp, const char *comment );
void Output_PreparedPatch_Fluid( const int TLv, const int TPID,
                                 const real h_Flu_Array[][FLU_NIN][ CUBE(FLU_NXT) ],
                                 const real h_Mag_Array[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                 const int NPG, const int *PID0_List, const int CLv, const char *comment );
#ifdef SUPPORT_FFTW
void Output_BasePowerSpectrum( const char *FileName, const long TVar );
#endif
void Output_L1Error( void (*AnalFunc_Flu)( real fluid[], const double x, const double y, const double z, const double Time,
                                           const int lv, double AuxArray[] ),
                     void (*AnalFunc_Mag)( real magnetic[], const double x, const double y, const double z, const double Time,
                                           const int lv, double AuxArray[] ),
                     const char *Prefix, const OptOutputPart_t Part, const double x, const double y, const double z );
#ifndef SERIAL
void Output_ExchangePatchMap( const int lv, const int xyz, const char *comment );
void Output_ExchangeFluxPatchList( const int option, const int lv, const char *comment );
void Output_ExchangeDataPatchList( const int option, const int lv, const char *comment );
void Output_BoundaryFlagList( const int option, const int lv, const char *comment );
#endif


// Refine
void FindFather( const int lv, const int Mode );
void Flag_Real( const int lv, const UseLBFunc_t UseLBFunc );
bool Flag_Check( const int lv, const int PID, const int i, const int j, const int k, const real dv,
                 const real Fluid[][PS1][PS1][PS1], const real Pot[][PS1][PS1], const real MagCC[][PS1][PS1][PS1],
                 const real Vel[][PS1][PS1][PS1], const real Pres[][PS1][PS1],  const real Lrtz[][PS1][PS1],
                 const real *Lohner_Var, const real *Lohner_Ave, const real *Lohner_Slope, const int Lohner_NVar,
                 const real ParCount[][PS1][PS1], const real ParDens[][PS1][PS1], const real JeansCoeff,
                 const real *Interf_Var, const real Spectral_Cond );
bool Flag_Lohner( const int i, const int j, const int k, const OptLohnerForm_t Form, const real *Var1D, const real *Ave1D,
                  const real *Slope1D, const int NVar, const double Threshold, const double Filter, const double Soften );
void Refine( const int lv, const UseLBFunc_t UseLBFunc );
void SiblingSearch( const int lv );
void SiblingSearch_Base();
#ifndef SERIAL
void Flag_Buffer( const int lv );
void Refine_Buffer( const int lv, const int *SonTable, const int *GrandTable );
#endif
#if ( ELBDM_SCHEME == ELBDM_HYBRID )
void Sync_UseWaveFlag( const int lv );
#endif


// SelfGravity
#ifdef GRAVITY
void CPU_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                               const real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                     real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                     real h_Flu_Array    [][GRA_NIN][PS1][PS1][PS1],
                               const double h_Corner_Array[][3],
                               const real h_Pot_Array_USG[][USG_NXT_G][USG_NXT_G][USG_NXT_G],
                               const real h_Flu_Array_USG[][GRA_NIN-1][PS1][PS1][PS1],
                                     char h_DE_Array     [][PS1][PS1][PS1],
                               const real h_Emag_Array   [][PS1][PS1][PS1],
                               const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter,
                               const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                               const int MG_NPre_Smooth, const int MG_NPost_Smooth, const real MG_Tolerated_Error,
                               const real Poi_Coeff, const IntScheme_t IntScheme, const bool P5_Gradient,
                               const real ELBDM_Eta, const real ELBDM_Lambda, const bool Poisson, const bool GraAcc,
                               const bool SelfGravity, const OptExtPot_t ExtPot, const OptExtAcc_t ExtAcc,
                               const double TimeNew, const double TimeOld, const real MinEint,
                               const bool UseWaveFlag );
void CPU_ExtPotSolver_BaseLevel( const ExtPot_t Func, const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real Table[], void **GenePtr,
                                 const double Time, const bool PotIsInit, const int SaveSg );
#ifdef SUPPORT_FFTW
void CPU_PoissonSolver_FFT( const real Poi_Coeff, const int SaveSg, const double PrepTime );
void Init_GreenFuncK();
#endif
void End_MemFree_PoissonGravity();
void Gra_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                    const int SaveSg_Flu, const int SaveSg_Pot, const bool Poisson, const bool Gravity,
                    const bool OverlapMPI, const bool Overlap_Sync, const bool Timing );
void Gra_Close( const int lv, const int SaveSg, const real h_Flu_Array_G[][GRA_NIN][PS1][PS1][PS1],
                const char h_DE_Array_G[][PS1][PS1][PS1], const real h_Emag_Array_G[][PS1][PS1][PS1],
                const int NPG, const int *PID0_List );
void Gra_Prepare_Flu( const int lv, real h_Flu_Array_G[][GRA_NIN][PS1][PS1][PS1], char h_DE_Array_G[][PS1][PS1][PS1],
                      real h_Emag_Array_G[][PS1][PS1][PS1], const int NPG, const int *PID0_List );
void Gra_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                      const int NPG, const int *PID0_List );
void Gra_Prepare_Corner( const int lv, double h_Corner_Array[][3], const int NPG, const int *PID0_List );
#ifdef UNSPLIT_GRAVITY
void Gra_Prepare_USG( const int lv, const double PrepTime,
                      real h_Pot_Array_USG_G[][USG_NXT_G][USG_NXT_G][USG_NXT_G],
                      real h_Flu_Array_USG_G[][GRA_NIN-1][PS1][PS1][PS1], const int NPG, const int *PID0_List );
#endif
void Init_ExtAccPot();
void End_ExtAccPot();
void Init_LoadExtPotTable();
void Init_MemAllocate_PoissonGravity( const int Pot_NPatchGroup );
void Init_Set_Default_MG_Parameter();
void Init_Set_Default_SOR_Parameter();
void Output_PreparedPatch_Poisson( const int TLv, const int TPID, const int TComp,
                                   const real h_Rho_Array_P   [][RHO_NXT][RHO_NXT][RHO_NXT],
                                   const real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT],
                                   const int NPG, const int *PID0_List, const int CLv, const char *comment );
void Poi_Close( const int lv, const int SaveSg, const real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                const int NPG, const int *PID0_List );
void Poi_BoundaryCondition_Extrapolation( real *Array, const int BC_Face, const int NVar, const int GhostSize,
                                          const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                          const int Idx_Start[], const int Idx_End[] );
void Poi_GetAverageDensity();
void Poi_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_In[][POT_NXT][POT_NXT][POT_NXT],
                      const int NPG, const int *PID0_List );
void Poi_Prepare_Rho( const int lv, const double PrepTime, real h_Rho_Array_P[][RHO_NXT][RHO_NXT][RHO_NXT],
                      const int NPG, const int *PID0_List );
#ifdef STORE_POT_GHOST
void Poi_StorePotWithGhostZone( const int lv, const int PotSg, const bool AllPatch );
#endif
#endif // #ifdef GRAVITY


// Tables
template <typename T> T TABLE_01( const int SibIndex, const char dim, const T w0, const T w1, const T w2 );
template <typename T> T TABLE_02( const int LocalID, const char dim, const T w0, const T w1 );
int TABLE_03( const int SibID, const int Count );
int TABLE_04( const int SibID );
int TABLE_05( const int SibID );
int TABLE_06( const int SibID, const int FlagLayer );
int TABLE_07( const int SibID, const int Count );
void TABLE_SiblingSharingSameEdge( const int EdgeID, int SibID[], int SibSibID[] );
void TABLE_GetSibPID_Delta( int NSibPatch[], int *SibPID_Delta[] );
void TABLE_GetSibPID_Based( const int lv, const int PID0, int SibPID_Based[] );


// LoadBalance
long LB_Corner2Index( const int lv, const int Corner[], const Check_t Check );
void LB_GetPID( const int GID, int& level, int& PID, int* GID_Offset );
void LB_AllgatherPatchCount( LB_PatchCount& pc );
void LB_AllgatherLBIdx( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList* gel = NULL );
void LB_FillLocalPatchExchangeList( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel );
void LB_FillGlobalPatchExchangeList( LB_PatchCount& pc, LB_LocalPatchExchangeList& lel, LB_GlobalPatchExchangeList& gel, int root );
LB_GlobalPatch* LB_ConstructGlobalTree( LB_PatchCount& pc, LB_GlobalPatchExchangeList& gel, int root );
LB_GlobalPatch* LB_GatherTree( LB_PatchCount& pc, int root );

#ifdef LOAD_BALANCE
void LB_AllocateBufferPatch_Father( const int SonLv, const bool SearchAllSon, const int NInput, int* TargetSonPID0,
                                    const bool RecordFaPID, int* NNewFaBuf0, int** NewFaBufPID0 );
void LB_AllocateBufferPatch_Sibling_Base();
void LB_AllocateBufferPatch_Sibling( const int lv );
void LB_AllocateFluxArray( const int FaLv );
void LB_ExchangeFlaggedBuffer( const int lv );
void LB_FindFather( const int SonLv, const bool SearchAllSon, const int NInput, int* TargetSonPID0, const bool ResetSonID );
void LB_FindSonNotHome( const int FaLv, const bool SearchAllFa, const int NInput, int* TargetFaPID );
void LB_GetBufferData( const int lv, const int FluSg, const int MagSg, const int PotSg, const GetBufMode_t GetBufMode,
                       const long TVarCC, const long TVarFC, const int ParaBuf );
void*LB_GetBufferData_MemAllocate_Send( const long SendSize );
void*LB_GetBufferData_MemAllocate_Recv( const long RecvSize );
void LB_GrandsonCheck( const int lv );
void LB_Init_LoadBalance( const bool Redistribute, const bool SendGridData, const double ParWeight, const bool Reset,
                          const bool SortRealPatch, const int TLv );
void LB_Init_ByFunction();
void LB_Init_Refine( const int FaLv, const bool AllocData );
void LB_SetCutPoint( const int lv, const int NPG_Total, long *CutPoint, const bool InputLBIdx0AndLoad,
                     long *LBIdx0_AllRank_Input, double *Load_AllRank_Input, const double ParWeight );
void LB_EstimateWorkload_AllPatchGroup( const int lv, const double ParWeight, double *Load_PG );
double LB_EstimateLoadImbalance();
void LB_SetCutPoint( const int lv, long *CutPoint, const bool InputLBIdx0AndLoad, long *LBIdx0_AllRank_Input,
                     double *Load_AllRank_Input, const double ParWeight );
void LB_Output_LBIdx( const int lv );
void LB_RecordExchangeDataPatchID( const int Lv, const bool AfterRefine );
void LB_RecordExchangeFixUpDataPatchID( const int Lv );
void LB_RecordExchangeRestrictDataPatchID( const int FaLv );
void LB_RecordOverlapMPIPatchID( const int Lv );
void LB_Refine( const int FaLv );
void LB_SiblingSearch( const int lv, const bool SearchAllPID, const int NInput, int *TargetPID0 );
void LB_Index2Corner( const int lv, const long Index, int Corner[], const Check_t Check );
int  LB_Index2Rank( const int lv, const long LB_Idx, const Check_t Check );
#endif // #ifdef LOAD_BALANCE


// Hydro model
#if    ( MODEL == HYDRO )
void Hydro_Aux_Check_Negative( const int lv, const int Mode, const char *comment );
void Hydro_GetTimeStep_Gravity( double &dt, double &dTime, int &MinDtLv, real &MinDtVar, const double dt_dTime );
void Hydro_GetMaxAcc( real MaxAcc[] );
void Hydro_Init_ByFunction_AssignData( const int lv );
void Hydro_BoundaryCondition_Reflecting( real *Array, const int BC_Face, const int NVar_Flu, const int GhostSize,
                                         const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                         const int Idx_Start[], const int Idx_End[], const int TFluVarIdxList[],
                                         const int NVar_Der, const long TDerVarList[] );
void Hydro_BoundaryCondition_Diode( real *Array, const int BC_Face, const int NVar_Flu, const int GhostSize,
                                    const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                    const int Idx_Start[], const int Idx_End[], const int TFluVarIdxList[],
                                    const int NVar_Der, const long TDerVarList[] );
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactorPtr );
void Hydro_Pri2Con( const real In[], real Out[], const bool FracPassive, const int NFrac, const int FracIdx[],
                    const EoS_DP2E_t EoS_DensPres2Eint, const EoS_TEM2H_t EoS_Temp2HTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], const real* const EintIn );
#ifdef MHD
void MHD_GetCellCenteredBField( real B_CC[], const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                const int Nx, const int Ny, const int Nz, const int i, const int j, const int k );
real MHD_GetCellCenteredBEnergy( const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                                 const int Nx, const int Ny, const int Nz, const int i, const int j, const int k );
void MHD_GetCellCenteredBFieldInPatch( real B[], const int lv, const int PID, const int i, const int j, const int k,
                                       const int MagSg );
real MHD_GetCellCenteredBEnergyInPatch( const int lv, const int PID, const int i, const int j, const int k,
                                        const int MagSg );
real MHD_GetCellCenteredDivBInPatch( const int lv, const int PID, const int i, const int j, const int k,
                                     const int MagSg );
void MHD_InterpolateBField( const real **CData, const int CSize[3][3], const int CStart[3][3], const int CRange[3],
                                  real **FData, const int FSize[3][3], const int FStart[3][3],
                            const real *FInterface[6], const IntScheme_t IntScheme, const bool Monotonic );
void MHD_AllocateElectricArray( const int lv );
void MHD_Aux_Check_InterfaceB( const char *comment );
void MHD_Aux_Check_DivergenceB( const bool Verbose, const char *comment );
void MHD_FixUp_Electric( const int lv );
void MHD_SameInterfaceB( const int lv );
void MHD_CopyPatchInterfaceBField( const int lv, const int PID, const int SibID, const int MagSg );
void MHD_BoundaryCondition_Outflow( real **Array, const int BC_Face, const int NVar, const int GhostSize,
                                    const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                    const int Idx_Start[], const int Idx_End[], const int TVarIdxList[] );
void MHD_BoundaryCondition_Reflecting( real **Array, const int BC_Face, const int NVar, const int GhostSize,
                                       const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                       const int Idx_Start[], const int Idx_End[], const int TVarIdxList[] );
void MHD_BoundaryCondition_Diode( real **Array, const int BC_Face, const int NVar, const int GhostSize,
                                  const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                  const int Idx_Start[], const int Idx_End[], const int TVarIdxList[] );
void MHD_BoundaryCondition_User( real **Array, const int BC_Face, const int NVar,
                                 const int ArraySizeX, const int ArraySizeY, const int ArraySizeZ,
                                 const int Idx_Start[], const int Idx_End[], const int TVarIdxList[],
                                 const double Time, const double dh, const double *Corner, const int lv );
void MHD_Init_BField_ByVecPot_File( const int B_lv );
void MHD_Init_BField_ByVecPot_Function( const int B_lv );
real MHD_ResetByUser_A2Bx( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh );
real MHD_ResetByUser_A2By( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh );
real MHD_ResetByUser_A2Bz( const real *Ax, const real *Ay, const real *Az,
                           const int i, const int j, const int k, const double dh );
#ifdef LOAD_BALANCE
void MHD_LB_EnsureBFieldConsistencyAfterRestrict( const int lv );
void MHD_LB_AllocateElectricArray( const int FaLv );
void MHD_LB_ResetBufferElectric( const int lv );
#endif
#endif // #ifdef MHD


// ELBDM model
#elif ( MODEL == ELBDM )
void   ELBDM_Init_ByFunction_AssignData( const int lv );
double ELBDM_GetTimeStep_Fluid( const int lv );
double ELBDM_GetTimeStep_Gravity( const int lv );
double ELBDM_GetTimeStep_Phase( const int lv );
#if   ( ELBDM_SCHEME == ELBDM_HYBRID )
double ELBDM_GetTimeStep_Hybrid_CFL( const int lv );
double ELBDM_GetTimeStep_Hybrid_Velocity( const int lv );
bool   ELBDM_HasWaveCounterpart( const int I, const int J, const int K, const long GID0, const long GID, const LB_GlobalTree& GlobalTree);
void   ELBDM_Aux_Record_Hybrid();
#endif

bool   ELBDM_Flag_EngyDensity( const int i, const int j, const int k, const real Real_Array[],
                               const real Imag_Array[], const double Angle_2pi, const double Eps );
bool   ELBDM_Flag_Interference( const int i, const int j, const int k, const real Var1D[], const double QPThreshold,
                                const double DensThreshold, const double LapPhaseThreshold, const bool OnlyAtExtrema );
real   ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped );
real   ELBDM_SetTaylor3Coeff( const real dt, const real dh, const real Eta );
void   ELBDM_RemoveMotionCM();
#ifdef SUPPORT_FFTW
void   CPU_ELBDMSolver_FFT( const real dt, const double PrepTime, const int SaveSg );
#endif
#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
void   ELBDM_GramFE_ComputeTimeEvolutionMatrix( gramfe_matmul_float (*output)[ 2*FLU_NXT ], const real dt, const real dh, const real Eta );
#endif


#else
#error : ERROR : unsupported MODEL !!
#endif


// GPU API
#ifdef GPU
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In[][FLU_NIN ][ CUBE(FLU_NXT) ],
                             real h_Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                             real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                             real h_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                             char h_DE_Array_Out[][ CUBE(PS2) ],
                             real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                             real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                             const double h_Corner_Array[][3],
                             real h_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
                             const bool h_IsCompletelyRefined[],
                             const bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                             const int NPatchGroup, const real dt, const real dh,
                             const bool StoreFlux, const bool StoreElectric,
                             const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter,
                             const real ELBDM_Eta, real ELBDM_Taylor3_Coeff, const bool ELBDM_Taylor3_Auto,
                             const double Time, const bool UsePot, const OptExtAcc_t ExtAcc, const MicroPhy_t MicroPhy,
                             const real MinDens, const real MinPres, const real MinEint,
                             const real DualEnergySwitch,
                             const bool NormPassive, const int NNorm,
                             const bool FracPassive, const int NFrac,
                             const bool JeansMinPres, const real JeansMinPres_Coeff,
                             const int GPU_NStream, const bool UseWaveFlag );
void CUAPI_Asyn_dtSolver( const Solver_t TSolver, real h_dt_Array[], const real h_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                          const real h_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const real h_Pot_Array[][ CUBE(GRA_NXT) ],
                          const double h_Corner_Array[][3], const int NPatchGroup, const real dh, const real Safety,
                          const MicroPhy_t MicroPhy, const real MinPres, const bool P5_Gradient, const bool UsePot,
                          const OptExtAcc_t ExtAcc, const double TargetTime, const int GPU_NStream );
void CUAPI_Asyn_SrcSolver( const real h_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
                                 real h_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
                           const real h_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
                           const double h_Corner_Array[][3],
                           const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
                           const double TimeNew, const double TimeOld,
                           const real MinDens, const real MinPres, const real MinEint,
                           const int GPU_NStream );
void CUAPI_DiagnoseDevice();
void CUAPI_MemAllocate();
void CUAPI_MemFree_Fluid( const int GPU_NStream );
void CUAPI_SetCache();
void CUAPI_SetDevice( const int Mode );
void CUAPI_SetConstMemory();
void CUAPI_SetConstMemory_EoS();
void CUAPI_Synchronize();
#ifdef GRAVITY
void CUAPI_SetConstMemory_ExtAccPot();
void CUAPI_Asyn_PoissonGravitySolver( const real h_Rho_Array    [][RHO_NXT][RHO_NXT][RHO_NXT],
                                      const real h_Pot_Array_In [][POT_NXT][POT_NXT][POT_NXT],
                                            real h_Pot_Array_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                                            real h_Flu_Array    [][GRA_NIN][PS1][PS1][PS1],
                                      const double h_Corner_Array[][3],
                                      const real h_Pot_Array_USG[][USG_NXT_G][USG_NXT_G][USG_NXT_G],
                                      const real h_Flu_Array_USG[][GRA_NIN-1][PS1][PS1][PS1],
                                            char h_DE_Array     [][PS1][PS1][PS1],
                                      const real h_Emag_Array   [][PS1][PS1][PS1],
                                      const int NPatchGroup, const real dt, const real dh, const int SOR_Min_Iter,
                                      const int SOR_Max_Iter, const real SOR_Omega, const int MG_Max_Iter,
                                      const int MG_NPre_Smooth, const int MG_NPost_Smooth,
                                      const real MG_Tolerated_Error, const real Poi_Coeff,
                                      const IntScheme_t IntScheme, const bool P5_Gradient, const real ELBDM_Eta,
                                      const real ELBDM_Lambda, const bool Poisson, const bool GraAcc,
                                      const bool SelfGravity, const OptExtPot_t ExtPot, const OptExtAcc_t ExtAcc,
                                      const double TimeNew, const double TimeOld, const real MinEint,
                                      const int GPU_NStream, const bool UseWaveFlag );
void CUAPI_SendExtPotTable2GPU( const real *h_Table );
void CUAPI_MemFree_PoissonGravity();
#endif // #ifdef GRAVITY
#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
void CUAPI_SendGramFEMatrix2GPU( gramfe_matmul_float (*h_GramFE_TimeEvo)[ 2*FLU_NXT ] );
#endif
#endif // #ifdef GPU


// Particle
#ifdef PARTICLE
void Par_Init_ByFile();
void Par_Output_TextFile( const char *comment );
void Par_Output_BinaryFile( const char *comment );
void Par_FindHomePatch_UniformGrid( const int lv, const bool OldParOnly, const long NNewPar,
                                    real_par *NewParAttFlt[PAR_NATT_FLT_TOTAL],
                                    long_par *NewParAttInt[PAR_NATT_INT_TOTAL] );
void Par_PassParticle2Son_SinglePatch( const int FaLv, const int FaPID );
void Par_PassParticle2Son_MultiPatch( const int FaLv, const ParPass2Son_t Mode, const bool TimingSendPar,
                                      const int NFaPatch, const int *FaPIDList );
void Par_PassParticle2Father( const int FaLv, const int FaPID );
void Par_Aux_Check_Particle( const char *comment );
void Par_MassAssignment( const long *ParList, const long NPar, const ParInterp_t IntScheme, real *Rho,
                         const int RhoSize, const double *EdgeL, const double dh, const bool PredictPos,
                         const double TargetTime, const bool InitZero, const bool Periodic[], const int PeriodicSize[3],
                         const bool UnitDens, const bool CheckFarAway, const bool UseInputMassPos, real_par **InputMassPos,
                         long_par **InputType );
void Par_UpdateParticle( const int lv, const double TimeNew, const double TimeOld, const ParUpStep_t UpdateStep,
                         const bool StoreAcc, const bool UseStoredAcc );
void Par_UpdateTracerParticle( const int lv, const double TimeNew, const double TimeOld,
                               const bool MapOnly );
void Par_GetTimeStep_VelAcc( double &dt_vel, double &dt_acc, const int lv );
void Par_PassParticle2Sibling( const int lv, const bool TimingSendPar );
bool Par_WithinActiveRegion( const real_par x, const real_par y, const real_par z );
int  Par_CountParticleInDescendant( const int FaLv, const int FaPID );
void Par_Aux_GetConservedQuantity( double &Mass, double &CoMX, double &CoMY, double &CoMZ,
                                   double &MomX, double &MomY, double &MomZ,
                                   double &AngMomX, double &AngMomY, double &AngMomZ, double &Ek, double &Ep );
void Par_Aux_InitCheck();
void Par_Aux_Record_ParticleCount();
void Par_CollectParticle2OneLevel( const int FaLv, const long FltAttBitIdx, const long IntAttBitIdx,
                                   const bool PredictPos, const double TargetTime, const bool SibBufPatch,
                                   const bool FaSibBufPatch, const bool JustCountNPar, const bool TimingSendPar );
void Par_CollectParticle2OneLevel_FreeMemory( const int FaLv, const bool SibBufPatch, const bool FaSibBufPatch );
int  Par_Synchronize( const double SyncTime, const ParSync_t SyncOption );
void Par_Synchronize_Restore( const double SyncTime );
void Prepare_PatchData_InitParticleDensityArray( const int lv, const double PrepTime );
void Prepare_PatchData_FreeParticleDensityArray( const int lv );
void Par_PredictPos( const long NPar, const long *ParList, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                     const double TargetTime );
void Par_Init_Attribute();
void Par_AddParticleAfterInit( const long NNewPar, real_par *NewParAttFlt[PAR_NATT_FLT_TOTAL], long_par *NewParAttInt[PAR_NATT_INT_TOTAL] );
void Par_ScatterParticleData( const long NPar_ThisRank, const long NPar_AllRank, const long FltAttBitIdx, const long IntAttBitIdx,
                              real_par *Data_Send_Flt[PAR_NATT_FLT_TOTAL], long_par *Data_Send_Int[PAR_NATT_INT_TOTAL],
                              real_par *Data_Recv_Flt[PAR_NATT_FLT_TOTAL], long_par *Data_Recv_Int[PAR_NATT_INT_TOTAL] );
void Par_MapMesh2Particles( const double EdgeL[3], const double EdgeR[3],
                            const double _dh, const int AttrSize3D, const real *Attr,
                            const int NPar, real_par *InterpParPos[3],
                            const long_par ParType[], const long ParList[],
                            const bool UseTracers, real_par ParAttr[], const bool CorrectVelocity );
void Par_Init_Attribute_Mesh();
void Par_Output_TracerParticle_Mesh();
FieldIdx_t AddParticleAttributeFlt( const char *InputLabel );
FieldIdx_t AddParticleAttributeInt( const char *InputLabel );
FieldIdx_t GetParticleAttributeFltIndex( const char *InputLabel, const Check_t Check );
FieldIdx_t GetParticleAttributeIntIndex( const char *InputLabel, const Check_t Check );
#ifdef LOAD_BALANCE
void Par_LB_CollectParticle2OneLevel( const int FaLv, const long FltAttBitIdx, const long IntAttBitIdx, const bool PredictPos,
                                      const double TargetTime, const bool SibBufPatch, const bool FaSibBufPatch,
                                      const bool JustCountNPar, const bool TimingSendPar );
void Par_LB_CollectParticle2OneLevel_FreeMemory( const int lv, const bool SibBufPatch, const bool FaSibBufPatch );
void Par_LB_CollectParticleFromRealPatch( const int lv, const long FltAttBitIdx, const long IntAttBitIdx,
                                          const int Buff_NPatchTotal, const int *Buff_PIDList, int *Buff_NPatchEachRank,
                                          const int Real_NPatchTotal, const int *Real_PIDList, int *Real_NPatchEachRank,
                                          const bool PredictPos, const double TargetTime,
                                          Timer_t *Timer, const char *Timer_Comment );
void Par_LB_ExchangeParticleBetweenPatch( const int lv,
                                          const int Send_NPatchTotal, const int *Send_PIDList, int *Send_NPatchEachRank,
                                          const int Recv_NPatchTotal, const int *Recv_PIDList, int *Recv_NPatchEachRank,
                                          Timer_t *Timer, const char *Timer_Comment );
void Par_LB_SendParticleData( const int NParAttFlt, const int NParAttInt, int *SendBuf_NPatchEachRank, int *SendBuf_NParEachPatch,
                              long *SendBuf_LBIdxEachPatch, real_par *SendBuf_ParFltDataEachPatch, long_par *SendBuf_ParIntDataEachPatch,
                              const long NSendParTotal, int *&RecvBuf_NPatchEachRank, int *&RecvBuf_NParEachPatch, long *&RecvBuf_LBIdxEachPatch,
                              real_par *&RecvBuf_ParFltDataEachPatch, long_par *&RecvBuf_ParIntDataEachPatch, int &NRecvPatchTotal, long &NRecvParTotal,
                              const bool Exchange_NPatchEachRank, const bool Exchange_LBIdxEachRank,
                              const bool Exchange_ParDataEachRank, Timer_t *Timer, const char *Timer_Comment );
void Par_LB_RecordExchangeParticlePatchID( const int MainLv );
void Par_LB_MapBuffer2RealPatch( const int lv, const int  Buff_NPatchTotal, int *&Buff_PIDList, int *Buff_NPatchEachRank,
                                                     int &Real_NPatchTotal, int *&Real_PIDList, int *Real_NPatchEachRank,
                                 const bool UseInputLBIdx, long *Buff_LBIdxList_Input );
#endif // #ifdef LOAD_BALANCE
#endif // #ifdef PARTICLE


// yt inline analysis
#ifdef SUPPORT_LIBYT
void YT_Init( int argc, char *argv[] );
void YT_End();
void YT_Inline();
#endif // #ifdef SUPPORT_LIBYT


// source terms
void Src_Init();
void Src_End();
void Src_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                    const int SaveSg_Flu, const int SaveSg_Mag, const bool OverlapMPI, const bool Overlap_Sync );
void Src_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_S_In[][FLU_NIN_S][ CUBE(SRC_NXT)           ],
                  real h_Mag_Array_S_In[][NCOMP_MAG][ SRC_NXT_P1*SQR(SRC_NXT) ],
                  double h_Corner_Array_S[][3],
                  const int NPG, const int *PID0_List );
void Src_Close( const int lv, const int SaveSg_Flu, const real h_Flu_Array_S_Out[][FLU_NOUT_S][ CUBE(PS1) ],
                const int NPG, const int *PID0_List );
void Src_WorkBeforeMajorFunc( const int lv, const double TimeNew, const double TimeOld, const double dt );
void CPU_SrcSolver( const real h_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
                          real h_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
                    const real h_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
                    const double h_Corner_Array[][3],
                    const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
                    const double TimeNew, const double TimeOld,
                    const real MinDens, const real MinPres, const real MinEint );


// Grackle
#ifdef SUPPORT_GRACKLE
void Grackle_Init();
void Grackle_Init_FieldData();
void Grackle_End();
void Init_MemAllocate_Grackle( const int Che_NPG );
void End_MemFree_Grackle();
void Grackle_Prepare( const int lv, real_che h_Che_Array[], const int NPG, const int *PID0_List );
void Grackle_Close( const int lv, const int SaveSg, const real_che h_Che_Array[], const int NPG, const int *PID0_List );
void Grackle_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt, const int SaveSg,
                        const bool OverlapMPI, const bool Overlap_Sync );
void CPU_GrackleSolver( grackle_field_data *Che_FieldData, code_units Che_Units, const int NPatchGroup, const real dt );
#endif // #ifdef SUPPORT_GRACKLE


// star formation
#ifdef STAR_FORMATION
void SF_CreateStar( const int lv, const real TimeNew, const real dt );
void SF_FreeRNG();
void SF_CreateStar_AGORA( const int lv, const real TimeNew, const real dt, RandomNumber_t *RNG,
                          const real GasDensThres, const real Efficiency, const real MinStarMass, const real MaxStarMFrac,
                          const bool DetRandom, const bool UseMetal );
#endif


// feedback
#ifdef FEEDBACK
void FB_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag );
void FB_Init();
void FB_End();
int FB_Aux_CellPatchRelPos( const int ijk[] );
#endif


// EoS in hydrodynamics
#if ( MODEL == HYDRO )
void EoS_Init();
void EoS_End();
#endif



#endif // __PROTOTYPE_H__
