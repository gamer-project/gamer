#ifndef __SRC_TERMS_H__
#define __SRC_TERMS_H__



void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );




//-------------------------------------------------------------------------------------------------------
// Structure   :  SrcTerms_t
// Description :  Data structure storing the source-term options to be passed to the CPU/GPU solvers
//
// Data Member :  Any             : True if at least one of the source terms is activated
//                Deleptonization : SRC_DELEPTONIZATION
//                User            : SRC_USER
//
// Method      :  SrcTerms_t : Constructor
//               ~SrcTerms_t : Destructor
//-------------------------------------------------------------------------------------------------------
struct SrcTerms_t
{

// data members
// ===================================================================================
   bool Any;
   bool Deleptonization;
   bool User;


   //===================================================================================
   // Constructor :  SrcTerms_t
   // Description :  Constructor of the structure "SrcTerms_t"
   //
   // Note        :  Initialize all data members
   //
   // Parameter   :  None
   //===================================================================================
   SrcTerms_t()
   {

      Any             = false;
      Deleptonization = false;
      User            = false;

   } // METHOD : SrcTerms_t



   //===================================================================================
   // Destructor  :  ~SrcTerms_t
   // Description :  Destructor of the structure "SrcTerms_t"
   //
   // Note        :
   //===================================================================================
   ~SrcTerms_t()
   {

   } // METHOD : ~SrcTerms_t


}; // struct SrcTerms_t



#endif // #ifndef __SRC_TERMS_H__
