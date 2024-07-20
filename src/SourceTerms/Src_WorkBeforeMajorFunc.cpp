#include "GAMER.h"



// prototypes of built-in source terms
#if ( MODEL == HYDRO )
void Src_WorkBeforeMajorFunc_Deleptonization( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                              double AuxArray_Flt[], int AuxArray_Int[] );
void Src_WorkBeforeMajorFunc_ExactCooling( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                           double AuxArray_Flt[], int AuxArray_Int[] );
#endif

// this function pointer can be set by a test problem initializer for a user-specified source term
void (*Src_WorkBeforeMajorFunc_User_Ptr)    ( const int lv, const double TimeNew, const double TimeOld, const double dt,
                                              double AuxArray_Flt[], int AuxArray_Int[] ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_WorkBeforeMajorFunc
// Description :  Work prior to the major source-term function
//
// Note        :  1. Invoked by Src_AdvanceDt()
//                2. Each source term can define its own Src_WorkBeforeMajorFunc_* function
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Target physical time to reach
//                TimeOld : Physical time before update
//                          --> The source-term function will update the system from TimeOld to TimeNew
//                dt      : Time interval to advance solution
//                          --> Physical coordinates : TimeNew - TimeOld == dt
//                              Comoving coordinates : TimeNew - TimeOld == delta(scale factor) != dt
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Src_WorkBeforeMajorFunc( const int lv, const double TimeNew, const double TimeOld, const double dt )
{

#  if ( MODEL == HYDRO )
// (1) deleptonization
   if ( SrcTerms.Deleptonization )
      Src_WorkBeforeMajorFunc_Deleptonization( lv, TimeNew, TimeOld, dt,
                                               Src_Dlep_AuxArray_Flt, Src_Dlep_AuxArray_Int );
// (2) exact cooling
   if ( SrcTerms.ExactCooling)
      Src_WorkBeforeMajorFunc_ExactCooling( lv, TimeNew, TimeOld, dt, Src_EC_AuxArray_Flt, Src_EC_AuxArray_Int );
#  endif

// (3) user-specified source term
// --> users may not define Src_WorkBeforeMajorFunc_User_Ptr
   if ( SrcTerms.User  &&  Src_WorkBeforeMajorFunc_User_Ptr != NULL )
      Src_WorkBeforeMajorFunc_User_Ptr       ( lv, TimeNew, TimeOld, dt,
                                               Src_User_AuxArray_Flt, Src_User_AuxArray_Int );

} // FUNCTION : Src_WorkBeforeMajorFunc
