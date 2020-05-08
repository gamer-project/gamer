#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPotAuxArray_PointMass
// Description :  Set the auxiliary array ExtPot_AuxArray[] used by the external potential routine ExtPot_PointMass()
//
// Note        :  1. External potential can be enabled by the runtime option "OPT__EXTERNAL_POT"
//                2. To adopt this routine, link to the function pointer "Init_ExtPotAuxArray_Ptr"
//                   in a test problem initializer as follows:
//
//                      void Init_ExtPotAuxArray_PointMass( double AuxArray[] );
//
//                      ...
//
//                      Init_ExtPotAuxArray_Ptr = Init_ExtPotAuxArray_PointMass;
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                3. AuxArray[] has the size of EXT_POT_NAUX_MAX defined in CUPOT.h (default = 10)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtPotAuxArray_PointMass( double AuxArray[] )
{

// example parameters
   const double M  = 1.0;
   const double GM = NEWTON_G*M;

   AuxArray[0] = 0.5*amr->BoxSize[0];  // x coordinate of the external potential center
   AuxArray[1] = 0.5*amr->BoxSize[1];  // y ...
   AuxArray[2] = 0.5*amr->BoxSize[2];  // z ...
   AuxArray[3] = GM;                   // gravitational_constant*point_source_mass

} // FUNCTION : Init_ExtPotAuxArray_PointMass



#endif // #ifdef GRAVITY
