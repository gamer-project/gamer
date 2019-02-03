#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_AddParticleAfterInit
// Description :  Add particles AFTER initialization
//
// Note        :  1. Unlike Par_Init_ByFile() and Par_Init_ByFunction() that are used for initialization,
//                   this function is mainly used AFTER initialization when the AMR structure and particle
//                   repository have been constructed
//                   --> One typical usage is to add new particles after restarting from a simulation
//                       snapshot
//
// Parameter   :  NNewPar   : Number of new particles to be added
//                NewParAtt : Pointer array storing the data of new particle attributes
//                            --> Format: real *NewParAtt[PAR_NATT_TOTAL]
//                            --> Must be deallocated manually after invoking this function
//
// Return      :  1. amr->Par
//                2. NPar, ParListSize, and ParList[] of all real patches on lv
//-------------------------------------------------------------------------------------------------------
void Par_AddParticleAfterInit( const long NNewPar, real *NewParAtt[PAR_NATT_TOTAL] )
{

   const bool OldParOnly_No    = false;
   const bool TimingSendPar_No = false;

// add new particles to the base level first
   Par_FindHomePatch_UniformGrid( 0, OldParOnly_No, NNewPar, NewParAtt );


// send particles to their home leaf patches
   for (int FaLv=0; FaLv<TOP_LEVEL; FaLv++)
   {
      const int NFaPatch = amr->NPatchComma[FaLv][1];

      int *FaPIDList = new int [NFaPatch];
      for (int t=0; t<NFaPatch; t++)   FaPIDList[t] = t;

      Par_PassParticle2Son_MultiPatch( FaLv, PAR_PASS2SON_GENERAL, TimingSendPar_No, NFaPatch, FaPIDList );

      delete [] FaPIDList;
   }

} // FUNCTION : Par_AddParticleAfterInit



#endif // #ifdef PARTICLE
