#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtAccAuxArray_PointMass
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by the external acceleration routine ExtAcc_PointMass()
//
// Note        :  1. External acceleration can be enabled by the runtime option "OPT__GRAVITY_TYPE = 2/3"
//                2. To adopt this routine, link to the function pointer "Init_ExtAccAuxArray_Ptr"
//                   in a test problem initializer as follows:
//
//                      void Init_ExtAccAuxArray_PointMass( double AuxArray[] );
//
//                      ...
//
//                      Init_ExtAccAuxArray_Ptr = Init_ExtAccAuxArray_PointMass;
//
//                   --> Then it will be invoked by Init_ExtAccPot()
//                3. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in CUPOT.h (default = 10)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtAccAuxArray_PointMass( double AuxArray[] )
{

// example parameters
   const double M   = 1.0;
   const double GM  = NEWTON_G*M;
   const double Eps = 0.0;

   AuxArray[0] = 0.5*amr->BoxSize[0];  // x coordinate of the external acceleration center
   AuxArray[1] = 0.5*amr->BoxSize[1];  // y ...
   AuxArray[2] = 0.5*amr->BoxSize[2];  // z ...
   AuxArray[3] = GM;                   // gravitational_constant*point_source_mass
   AuxArray[4] = Eps;                  // soften_length (<=0.0 --> disable)

} // FUNCTION : Init_ExtAccAuxArray_PointMass



#endif // #ifdef GRAVITY
