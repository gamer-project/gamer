#ifndef __SRC_TERMS_H__
#define __SRC_TERMS_H__



#include "EoS.h"

// forward declaration of SrcTerms_t since it is required by SrcFunc_t
// --> its content will be specified later
struct SrcTerms_t;

// typedef of the source-term function
typedef void (*SrcFunc_t)( real fluid[], const real B[],
                           const SrcTerms_t *SrcTerms, const real dt, const real dh,
                           const double x, const double y, const double z,
                           const double TimeNew, const double TimeOld,
                           const real MinDens, const real MinPres, const real MinEint, const long PassiveFloor,
                           const EoS_t *EoS, const double AuxArray_Flt[], const int AuxArray_Int[] );




//-------------------------------------------------------------------------------------------------------
// Structure   :  SrcTerms_t
// Description :  Data structure storing the source-term variables to be passed to the CPU/GPU solvers
//
// Data Member :  Any                       : True if at least one of the source terms is activated
//                Deleptonization           : SRC_DELEPTONIZATION
//                User                      : SRC_USER
//                BoxCenter                 : Simulation box center
//                Unit_*                    : Code units
//                *_FuncPtr                 : Major source-term functions
//                *_CPUPtr/*_GPUPtr         : CPU/GPU function pointers to the major source-term function
//                *_AuxArrayDevPtr_*        : Auxiliary array pointers
//                                            --> For GPU, these pointers store the addresses of constant memory arrays,
//                                                which should NOT be used by host
//                Dlep_Profile_DataDevPtr   : Profile array pointer used by deleptonization
//                Dlep_Profile_RadiusDevPtr : Radius array pointer associated with Dlep_Profile_DataDevPtr[]
//                                            --> For GPU, Dlep_Profile_DataDevPtr[]/RadiusDevPtr[] store the
//                                                addresses of global memory arrays, which should NOT be used by host
//                Dlep_Profile_NBin         : Number of radial bins in Dlep_Profile_*
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
#  if ( MODEL == HYDRO )
   SrcFunc_t Dlep_FuncPtr;
   SrcFunc_t Dlep_CPUPtr;
#  ifdef GPU
   SrcFunc_t Dlep_GPUPtr;
#  endif
   double   *Dlep_AuxArrayDevPtr_Flt;
   int      *Dlep_AuxArrayDevPtr_Int;
   real    (*Dlep_Profile_DataDevPtr)[SRC_DLEP_PROF_NBINMAX];
   real     *Dlep_Profile_RadiusDevPtr;
   int       Dlep_Profile_NBin;
#  endif

// user-specified source term
   SrcFunc_t User_FuncPtr;
   SrcFunc_t User_CPUPtr;
#  ifdef GPU
   SrcFunc_t User_GPUPtr;
#  endif
   double   *User_AuxArrayDevPtr_Flt;
   int      *User_AuxArrayDevPtr_Int;

}; // struct SrcTerms_t



#endif // #ifndef __SRC_TERMS_H__
