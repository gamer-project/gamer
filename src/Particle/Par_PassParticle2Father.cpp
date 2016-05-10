#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Father
// Description :  Pass particles from sons to father
//
// Note        :  1. After calling this function, son patches will have no particles (NPar == 0)
//                2. This function should always be called before deleting son patches within the same patch group
//                   --> invoked by "Refine & LB_Refine_AllocateNewPatch"
//
// Parameter   :  FaLv  : Father's refinement level
//                FaPID : Father's patch ID
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Father( const int FaLv, const int FaPID )
{

   const int SonPID0 = amr->patch[0][FaLv][FaPID]->son;
   const int SonLv   = FaLv + 1;


// nothing to do if father has no son
   if ( SonPID0 == -1  ||  FaLv == TOP_LEVEL )     return;


// 1. get the total number of particles in all sons
   int NParSon = 0;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)        NParSon += amr->patch[0][SonLv][SonPID]->NPar;


// nothing to do if sons have no particles
   if ( NParSon == 0 )  return;


// 2. gather the particle lists from all sons
   long *ParListSon = new long [NParSon];
   int  t = 0;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)  
   for (int p=0; p<amr->patch[0][SonLv][SonPID]->NPar; p++)    ParListSon[ t++ ] = amr->patch[0][SonLv][SonPID]->ParList[p];


// 3. add particles to father 
//###NOTE : No OpenMP since AddParticle will modify amr->Par->NPar_Lv[]
#  ifdef DEBUG_PARTICLE
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   amr->patch[0][FaLv][FaPID]->AddParticle( NParSon, ParListSon, &amr->Par->NPar_Lv[FaLv], ParPos, amr->Par->NPar, __FUNCTION__ );
#  else
   amr->patch[0][FaLv][FaPID]->AddParticle( NParSon, ParListSon, &amr->Par->NPar_Lv[FaLv] );
#  endif


// 4. remove particles in all sons
//###NOTE : No OpenMP since RemoveParticle will modify amr->Par->NPar_Lv[]
   const bool RemoveAllParticle = true;
   for (int SonPID=SonPID0; SonPID<SonPID0+8; SonPID++)
      amr->patch[0][SonLv][SonPID]->RemoveParticle( 0, NULL, &amr->Par->NPar_Lv[SonLv], RemoveAllParticle );


// free memory
   delete [] ParListSon;

} // FUNCTION : Par_PassParticle2Father



#endif // #ifdef PARTICLE
