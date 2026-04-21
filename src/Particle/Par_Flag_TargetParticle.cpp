#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Flag_TargetParticle
// Description :  Flag patches containing particles marked for refinement
//
// Note        :  1. Invoked by Flag_Real()
//                2. Support two modes controlled by FlagMode:
//                     PAR_FLAG_MUST: patches containing particles with PAR_FLAG > 0 *must* be refined to level +PAR_FLAG
//                     PAR_FLAG_CAN : patches containing particles with PAR_FLAG < 0 *can*  be refined to level -PAR_FLAG
//                3. Additional notes on PAR_FLAG_CAN:
//                   --> This behaves similarly to other refinement pre-checks (e.g., proper-nesting condition, OPT__FLAG_REGION, OPT__FLAG_ANGULAR).
//                       Even if a patch passes this pre-check, it must still satisfy at least one real refinement criterion
//                       (e.g., OPT__FLAG_RHO or PAR_FLAG_MUST) to trigger refinement.
//                   --> A patch cannot be refined if it contains no particles with PAR_FLAG < 0 (except due to flag buffers),
//                       regardless of other refinement pre-checks. Therefore, a patch is allowed to be refined only if it passes
//                       all pre-checks *simultanesouly*.
//                       --> This is consistent with the current implementation of OPT__FLAG_REGION and OPT__FLAG_ANGULAR.
//                   --> Unlike OPT__FLAG_REGION and OPT__FLAG_ANGULAR, this pre-check does not consider the spatial distribution of
//                       particles within a patch. In other words, even if the cell containing particles with PAR_FLAG < 0
//                       and the cell satisfying a real refinement criterion (e.g., OPT__FLAG_RHO) are different, this patch will
//                       still be refined. This simplifies the implementation and promote refinement.
//                4. When OPT__FLAG_PAR_TARGET == PAR_FLAG_BOTH, a patch cannot be refined if it contains no particles with PAR_FLAG < 0,
//                   even if it contains particles with PAR_FLAG > 0.
//                   --> This is consistent with the current implementation of refinement pre-checks. For example, a patch cannot be refined
//                       if it fails the OPT__FLAG_REGION pre-check, even if it satisfies the OPT__FLAG_RHO criterion.
//                5. Par_CollectParticle2OneLevel() must be called before this routine
//
// Parameter   :  FaLv     : Target refinement level for flagging
//                PID      : Target patch ID
//                FlagMode : See the "Note" section above
//
// Return      :  true  : The target patch contains particles flagged for refinement
//                        (i.e., FaLv < PAR_FLAG for PAR_FLAG_MUST or FaLv < -PAR_FLAG for PAR_FLAG_CAN)
//                false : Otherwise
//-------------------------------------------------------------------------------------------------------
bool Par_Flag_TargetParticle( const int FaLv, const int PID, const ParFlag_t FlagMode )
{

// checks
#  ifdef GAMER_DEBUG
   if ( FaLv < 0  ||  FaLv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect FaLv = %d (TOP_LEVEL = %d) !!\n", FaLv, TOP_LEVEL );

   if ( PID < 0  ||  PID >= amr->num[FaLv] )
      Aux_Error( ERROR_INFO, "incorrect PID = %d (amr->num[%d] = %d) !!\n", PID, FaLv, amr->num[FaLv] );

   if ( FlagMode != PAR_FLAG_MUST  &&  FlagMode != PAR_FLAG_CAN )
      Aux_Error( ERROR_INFO, "incorrect FlagMode = %d !!\n", FlagMode );
#  endif


   int   NPar, Flag;
   long *ParList = NULL;
   bool  UseCopy;

// leaf patch
   if ( amr->patch[0][FaLv][PID]->son == -1 )
   {
      NPar    = amr->patch[0][FaLv][PID]->NPar;
      ParList = amr->patch[0][FaLv][PID]->ParList;
      UseCopy = false;
   }

// non-leaf patch
   else
   {
      NPar    = amr->patch[0][FaLv][PID]->NPar_Copy;
#     ifdef LOAD_BALANCE
      ParList = NULL;
      UseCopy = true;
#     else
      ParList = amr->patch[0][FaLv][PID]->ParList_Copy;
      UseCopy = false;
#     endif
   } // if ( amr->patch[0][FaLv][PID]->son == -1 ) ... else ...


   for (int p=0; p<NPar; p++)
   {
      if ( UseCopy )    Flag = amr->patch[0][FaLv][PID]->ParAttInt_Copy[PAR_FLAG][p];
      else              Flag = amr->Par->Flag[ ParList[p] ];

//    return immediately if the patch satisfies the refinement checks
      if      ( FlagMode == PAR_FLAG_MUST ) {
         if ( Flag > 0  &&  FaLv < +Flag )   return true;
      }

      else if ( FlagMode == PAR_FLAG_CAN ) {
         if ( Flag < 0  &&  FaLv < -Flag )   return true;
      }

      else {
         Aux_Error( ERROR_INFO, "unsupport FlagMode (%d) !!\n", FlagMode );
      }
   } // for (int p=0; p<NPar; p++)


// this patch does not satisfy the refinement checks
   return false;

} // FUNCTION : Par_Flag_TargetParticle



#endif // #ifdef PARTICLE
