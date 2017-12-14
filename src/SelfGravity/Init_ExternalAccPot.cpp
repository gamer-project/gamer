#include "GAMER.h"

#ifdef GRAVITY


// these function pointers point to "Init_ExternalAcc()" and "Init_ExternalPot()" by default
// but may be overwritten by various test problem initializers
extern void (*Init_ExternalAcc_Ptr)();
extern void (*Init_ExternalPot_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalAccPot
// Description :  Set the auxiliary CPU/GPU arrays for the external acceleration and potential
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the runtime options "OPT__GRAVITY_TYPE == 2/3" and "OPT__EXTERNAL_POT"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalAccPot()
{

// initialize the auxiliary CPU arrays
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
      if ( Init_ExternalAcc_Ptr == NULL )    Aux_Error( ERROR_INFO, "Init_ExternalAcc_Ptr == NULL !!\n" );

      Init_ExternalAcc_Ptr();
   }

   if ( OPT__EXTERNAL_POT )
   {
      if ( Init_ExternalPot_Ptr == NULL )    Aux_Error( ERROR_INFO, "Init_ExternalPot_Ptr == NULL !!\n" );

      Init_ExternalPot_Ptr();
   }


// initialize the auxiliary GPU arrays
#  ifdef GPU
   CUAPI_Init_ExternalAccPot();
#  endif

} // FUNCTION : Init_ExternalAccPot



#endif // #ifdef GRAVITY
