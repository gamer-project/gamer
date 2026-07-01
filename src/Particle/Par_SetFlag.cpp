#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SetFlag
// Description :  Set the refinement flag of all particles to the specified value
//
// Note        :  1. Invoked by Init_GAMER()
//
// Parameter   :  Flag : Particle refinement flag
//
// Return      :  amr->Par->Flag
//-------------------------------------------------------------------------------------------------------
void Par_SetFlag( const int Flag )
{

// do not distinguish active and inactive particles for simplicity
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)   amr->Par->Flag[p] = (long_par)Flag;

} // FUNCTION : Par_SetFlag



#endif // #ifdef PARTICLE
