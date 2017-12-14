#include "GAMER.h"

#ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CountParticleInDescendant
// Description :  Count the number of particles in all descendants (sons, grandsons, ...) of the target patch
//
// Note        :  This function will search over all descendants recursively
//
// Parameter   :  FaLv  : Father patch level
//                FaPID : Father patch ID
//
// Return      :  NPar_Sum
//-------------------------------------------------------------------------------------------------------
int Par_CountParticleInDescendant( const int FaLv, const int FaPID )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "FaLv", FaLv );
#  endif


   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;

   int NPar_Sum = 0;

   if ( SonPID0 == -1 )    return 0;   // nothing to do if the target patch has no son
   else
   {
      for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      {
//       search over all descendants recursively
         if ( amr->patch[0][SonLv][SonPID]->son != -1 )  NPar_Sum += Par_CountParticleInDescendant( SonLv, SonPID );
         else                                            NPar_Sum += amr->patch[0][SonLv][SonPID]->NPar;
      }
   }

   return NPar_Sum;

} // FUNCTION : Par_CountParticleInDescendant



#endif // #ifdef PARTICLE
