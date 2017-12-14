#include "GAMER.h"

#ifdef GRAVITY


#include "CUPOT.h"
double ExtPot_AuxArray[EXT_POT_NAUX_MAX];

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Init_ExternalPot();

// this function pointer may be overwritten by various test problem initializers
void (*Init_ExternalPot_Ptr)() = Init_ExternalPot;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExternalPot
// Description :  Set the array "ExtPot_AuxArray" used by the external potential routines
//                "CUPOT_ExternalPot.cu / CPU_ExternalPot.cpp"
//
// Note        :  1. Invoked by "Init_GAMER" using the function pointer "Init_ExternalPot_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Enabled by the runtime option "OPT__EXTERNAL_POT"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExternalPot()
{

// ExtPot_AuxArray has the size of EXT_POT_NAUX_MAX defined in CUPOT.h (default = 10)
// --> by default we set
//        ExtPot_AuxArray[0] = x coordinate of the external potential center
//        ExtPot_AuxArray[1] = y ...
//        ExtPot_AuxArray[2] = z ...
//        ExtPot_AuxArray[3] = gravitational_constant*point_source_mass
// --> to change the this default behavior, please edit "GPU_Poisson/CUPOT_ExternalPot.cu"

   /*
   const double M  = 1.0;
   const double GM = NEWTON_G*M;

   ExtPot_AuxArray[0] = 0.5*amr->BoxSize[0];
   ExtPot_AuxArray[1] = 0.5*amr->BoxSize[1];
   ExtPot_AuxArray[2] = 0.5*amr->BoxSize[2];
   ExtPot_AuxArray[3] = GM;
   */

} // FUNCTION : Init_ExternalPot



#endif // #ifdef GRAVITY
