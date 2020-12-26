#ifndef __SRC_TERMS_H__
#define __SRC_TERMS_H__




//-------------------------------------------------------------------------------------------------------
// Structure   :  SrcTerms_t
// Description :  Data structure storing the source-term variables to be passed to the CPU/GPU solvers
//
// Data Member :  Any             : True if at least one of the source terms is activated
//                Deleptonization : SRC_DELEPTONIZATION
//                User            : SRC_USER
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct SrcTerms_t
{

   bool Any;
   bool Deleptonization;
   bool User;

}; // struct SrcTerms_t



#endif // #ifndef __SRC_TERMS_H__
