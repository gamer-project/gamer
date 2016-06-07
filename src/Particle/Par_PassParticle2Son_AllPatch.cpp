#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Son_AllPatch
// Description :  Pass particles from father to sons for all patches at the target level
//
// Note        :  1. It simply invokes Par_PassParticle2Son for all patches at the target level
//                2. It is invoked in EvolveLevel after the velocity correction in KDK
//
// Parameter   :  FaLv  : Father's refinement level
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Son_AllPatch( const int FaLv )
{

// nothing to do if there is no patch at FaLv+1
   if ( FaLv == TOP_LEVEL  ||  NPatchTotal[FaLv+1] == 0 )   return;

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->son != -1 )    Par_PassParticle2Son( FaLv, FaPID );
   }

} // FUNCTION : Par_PassParticle2Son_AllPatch 



#endif // #ifdef PARTICLE
