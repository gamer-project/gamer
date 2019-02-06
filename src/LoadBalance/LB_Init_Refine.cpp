#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_Refine
// Description :  Refine patches on FaLv to create patches on FaLv+1
//
// Note        :  1. Invoked by LB_Init_ByFunction()
//                2. Only apply to real patches
//                   --> Buffer patches will be allocated later when invoking LB_Init_LoadBalance()
//                3. Must inovke Flag_Real( FaLv ) in advance
//                4. LB_Init_LoadBalance() will be invoked later to optimize load-balancing
//
// Parameter   :  FaLv : Target level to be refined (FaLv --> FaLv+1)
//-------------------------------------------------------------------------------------------------------
void LB_Init_Refine( const int FaLv )
{

   const int  SonLv         = FaLv + 1;
   const int  PScale        = PS1*amr->scale[SonLv];
   const bool AllocData_Yes = true;

   const int *FaCr;


// check
   if ( FaLv == TOP_LEVEL )
      Aux_Error( ERROR_INFO, "refine the maximum level !!\n" );

   if ( amr->num[SonLv] != 0 )
      Aux_Error( ERROR_INFO, "number of son patches on level %d = %d != 0 !!\n", SonLv, amr->num[SonLv] );


// loop over all **real** patches on FaLv
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->flag )
      {
//       construct the father->son relation
         amr->patch[0][FaLv][FaPID]->son = amr->num[SonLv];


//       allocate child patches and construct the child->father relation
         FaCr = amr->patch[0][FaLv][FaPID]->corner;

         amr->pnew( SonLv, FaCr[0],        FaCr[1],        FaCr[2],        FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1],        FaCr[2],        FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0],        FaCr[1]+PScale, FaCr[2],        FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0],        FaCr[1],        FaCr[2]+PScale, FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1]+PScale, FaCr[2],        FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0],        FaCr[1]+PScale, FaCr[2]+PScale, FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1],        FaCr[2]+PScale, FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );
         amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1]+PScale, FaCr[2]+PScale, FaPID, AllocData_Yes, AllocData_Yes, AllocData_Yes );

//       pass particles from father to sons
#        ifdef PARTICLE
         Par_PassParticle2Son_SinglePatch( FaLv, FaPID );
#        endif

      } // if ( amr->patch[0][FaLv][FaPID]->flag )
   } // for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)

   for (int m=1; m<28; m++)   amr->NPatchComma[SonLv][m] = amr->num[SonLv];

} // FUNCTION : LB_Init_Refine



#endif // #ifdef LOAD_BALANCE
