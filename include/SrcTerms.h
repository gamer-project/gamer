#ifndef __SRC_TERMS_H__
#define __SRC_TERMS_H__


// forward declaration of SrcTerms_t since it is required by SrcFunc_t
// --> its content will be specified later
struct SrcTerms_t;

// typedef of the source-term function
typedef void (*SrcFunc_t)( real fluid[], const real B[],
                           const SrcTerms_t SrcTerms, const real dt, const real dh,
                           const double x, const double y, const double z,
                           const double TimeNew, const double TimeOld,
                           const real MinDens, const real MinPres, const real MinEint,
                           const EoS_DE2P_t EoS_DensEint2Pres,
                           const EoS_DP2E_t EoS_DensPres2Eint,
                           const EoS_DP2C_t EoS_DensPres2CSqr,
                           const double EoS_AuxArray_Flt[],
                           const int    EoS_AuxArray_Int[],
                           const real *const EoS_Table[EOS_NTABLE_MAX],
                           const double AuxArray_Flt[], const int AuxArray_Int[] );




//-------------------------------------------------------------------------------------------------------
// Structure   :  SrcTerms_t
// Description :  Data structure storing the source-term variables to be passed to the CPU/GPU solvers
//
// Data Member :  Any             : True if at least one of the source terms is activated
//                Deleptonization : SRC_DELEPTONIZATION
//                User            : SRC_USER
//                BoxCenter       : Simulation box center
//                Unit_*          : Code units
//                *_FuncPtr       : Major source-term functions
//                *_AuxArray_*    : Auxiliary arrays
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct SrcTerms_t
{

   bool   Any;
   bool   Deleptonization;
   bool   User;

   double BoxCenter[3];

   real   Unit_L;
   real   Unit_M;
   real   Unit_T;
   real   Unit_V;
   real   Unit_D;
   real   Unit_E;
   real   Unit_P;
#  ifdef MHD
   real   Unit_B;
#  endif

// deleptonization
   SrcFunc_t Dlep_FuncPtr;
   double    Dlep_AuxArray_Flt[SRC_NAUX_DLEP];
   int       Dlep_AuxArray_Int[SRC_NAUX_DLEP];

// user-specified source term
   SrcFunc_t User_FuncPtr;
   double    User_AuxArray_Flt[SRC_NAUX_USER];
   int       User_AuxArray_Int[SRC_NAUX_USER];

}; // struct SrcTerms_t



#endif // #ifndef __SRC_TERMS_H__
