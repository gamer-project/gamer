#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Flag_TargetParticle
// Description :  Flag patches containing particles marked for refinement
//
// Note        :  1. Invoked by Flag_Real()
//                2. Support two modes controlled by FlagMode:
//                     PAR_FLAG_MUST: patches containing particles with PAR_FLAG > 0 *must* be refined to level +PAR_FLAG
//                                    --> similar to other refinement criteria (e.g., OPT__FLAG_RHO)
//                                    --> called by Flag_Real()
//                     PAR_FLAG_CAN : patches containing particles with PAR_FLAG < 0 *can*  be refined to level -PAR_FLAG
//                                    --> similar to other refinement pre-checks (e.g., OPT__FLAG_REGION)
//                                    --> called by Flag_Precheck()
//                3. Also see the notes in Flag_Precheck()
//                4. When OPT__FLAG_PAR_TARGET == PAR_FLAG_BOTH, a patch cannot be refined if it contains no particles with PAR_FLAG < 0,
//                   even if it contains particles with PAR_FLAG > 0.
//                   --> This is consistent with the current implementation of refinement pre-checks. For example, a patch cannot be refined
//                       if it fails the OPT__FLAG_REGION pre-check, even if it satisfies the OPT__FLAG_RHO criterion.
//                5. Par_CollectParticle2OneLevel() must be called before this routine
//
// Parameter   :  lv       : Target refinement level for flagging
//                PID      : Target patch ID
//                FlagMode : See the "Note" section above
//
// Return      :  true  : The target patch contains particles flagged for refinement
//                        (i.e., lv < PAR_FLAG for PAR_FLAG_MUST or lv < -PAR_FLAG for PAR_FLAG_CAN)
//                false : Otherwise
//-------------------------------------------------------------------------------------------------------
bool Par_Flag_TargetParticle( const int lv, const int PID, const ParFlag_t FlagMode )
{

// checks
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d (TOP_LEVEL = %d) !!\n", lv, TOP_LEVEL );

   if ( PID < 0  ||  PID >= amr->num[lv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (amr->num[%d] = %d) !!\n", PID, lv, amr->num[lv] );

   if ( FlagMode != PAR_FLAG_MUST  &&  FlagMode != PAR_FLAG_CAN )
      Aux_Error( ERROR_INFO, "incorrect FlagMode = %d !!\n", FlagMode );
#  endif


   int   NPar, Flag;
   long *ParList = NULL;
   bool  UseCopy;

// leaf patch
   if ( amr->patch[0][lv][PID]->son == -1 )
   {
      NPar    = amr->patch[0][lv][PID]->NPar;
      ParList = amr->patch[0][lv][PID]->ParList;
      UseCopy = false;
   }

// non-leaf patch
   else
   {
      NPar    = amr->patch[0][lv][PID]->NPar_Copy;
#     ifdef LOAD_BALANCE
      ParList = NULL;
      UseCopy = true;
#     else
      ParList = amr->patch[0][lv][PID]->ParList_Copy;
      UseCopy = false;
#     endif
   } // if ( amr->patch[0][lv][PID]->son == -1 ) ... else ...


   for (int p=0; p<NPar; p++)
   {
      if ( UseCopy )    Flag = amr->patch[0][lv][PID]->ParAttInt_Copy[PAR_FLAG][p];
      else              Flag = amr->Par->Flag[ ParList[p] ];

//    return immediately if the patch satisfies the refinement checks
      if      ( FlagMode == PAR_FLAG_MUST ) {
         if ( Flag > 0  &&  lv < +Flag )   return true;
      }

      else if ( FlagMode == PAR_FLAG_CAN ) {
         if ( Flag < 0  &&  lv < -Flag )   return true;
      }

      else {
         Aux_Error( ERROR_INFO, "unsupport FlagMode (%d) !!\n", FlagMode );
      }
   } // for (int p=0; p<NPar; p++)


// this patch does not satisfy the refinement checks
   return false;

} // FUNCTION : Par_Flag_TargetParticle



#endif // #ifdef PARTICLE
