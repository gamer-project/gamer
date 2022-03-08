#include "ExtractProfile.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_UnwrapPhase
// Description :  Unwrap the input two phases to ensure that the phase difference is <= PI
//
// Note        :  Phase_Ref is fixed, and Phase_Wrapped will be unwrapped
//
// Parameter   :  Phase_Ref      : Reference phase
//                Phase_Wrapped  : Phase to be unwrapped
//
// Return      :  Phase_Unwrapped
//-------------------------------------------------------------------------------------------------------
real ELBDM_UnwrapPhase( const real Phase_Ref, const real Phase_Wrapped )
{

   const real TwoPi = 2.0*M_PI;

   real Phase_Unwrapped = Phase_Wrapped;

   while ( Phase_Unwrapped - Phase_Ref > +M_PI )   {  Phase_Unwrapped -= TwoPi;  };
   while ( Phase_Unwrapped - Phase_Ref < -M_PI )   {  Phase_Unwrapped += TwoPi;  };

   return Phase_Unwrapped;

} // FUNCTION : ELBDM_UnwrapPhase



#endif // #if ( MODEL == ELBDM )
