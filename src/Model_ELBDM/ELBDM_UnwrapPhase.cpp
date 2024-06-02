#include "GAMER.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_UnwrapWindingNumber
// Description :  Return the integer multiple of 2*PI that needs to be added to Phase_Wrapped to ensure that the phase difference is <= PI
//
// Note        :  Phase_Ref is fixed, and Phase_Wrapped will be unwrapped
//
// Parameter   :  Phase_Ref     : Reference phase
//                Phase_Wrapped : Phase to be unwrapped
//
// Return      :  phase winding number from Phase_Ref to Phase_Wrapped as integer
//-------------------------------------------------------------------------------------------------------
static int ELBDM_UnwrapWindingNumber( const real Phase_Ref, const real Phase_Wrapped )
{

   return ROUND( (Phase_Ref - Phase_Wrapped) / real(2.0*M_PI) );

} // FUNCTION : ELBDM_UnwrapWindingNumber



//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_UnwrapPhase
// Description :  Unwrap the input two phases to ensure that the phase difference is <= PI
//
// Note        :  Phase_Ref is fixed, and Phase_Wrapped will be unwrapped
//
// Parameter   :  Phase_Ref     : Reference phase
//                Phase_Wrapped : Phase to be unwrapped
//
// Return      :  Phase_Unwrapped
//-------------------------------------------------------------------------------------------------------
real ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped )
{
   return Phase_Wrapped + ELBDM_UnwrapWindingNumber(Phase_Ref, Phase_Wrapped) * real(2.0*M_PI);

} // FUNCTION : ELBDM_UnwrapPhase



#endif // #if ( MODEL == ELBDM )
