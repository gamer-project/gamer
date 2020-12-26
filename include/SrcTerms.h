#ifndef __SRC_TERMS_H__
#define __SRC_TERMS_H__




//-------------------------------------------------------------------------------------------------------
// Structure   :  SrcTerms_t
// Description :  Data structure storing the source-term variables to be passed to the CPU/GPU solvers
//
// Data Member :  Any             : True if at least one of the source terms is activated
//                Deleptonization : SRC_DELEPTONIZATION
//                User            : SRC_USER
//                BoxCenter       : Simulation box center
//                Unit_*          : Code units
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

}; // struct SrcTerms_t



#endif // #ifndef __SRC_TERMS_H__
