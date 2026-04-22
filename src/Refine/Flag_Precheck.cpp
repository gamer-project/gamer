#include "GAMER.h"

// defined in Flag_Check.cpp
bool Check_Angular_Max( const int i, const int j, const int k, const int lv, const int PID,
                        const double CenX, const double CenY, const double CenZ,
                        const double AngRes_Max, const double AngRes_Max_R );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Precheck
// Description :  Check whether a target patch is *allowed* to be refined
//
// Note        :  1. Called by Flag_Real()
//                2. Currently implemented pre-checks: proper-nesting condition, target particles, OPT__FLAG_REGION, OPT__FLAG_ANGULAR
//                3. Even if a patch passes all pre-checks, it must still satisfy at least one real refinement criterion
//                   (e.g., OPT__FLAG_RHO or PAR_FLAG_MUST) to trigger refinement
//                4. A patch is allowed to be refined only if it passes all pre-checks *simultanesouly*
//
// Parameter   :  lv          : Target refinement level for flagging
//                PID         : Target patch ID
//                NoRefineBnd : Boundary regions excluded from refinement (used by OPT__NO_FLAG_NEAR_BOUNDARY)
//
// Return      :  true  : All refinement pre-checks pass
//                false : Any refinement pre-check fails
//-------------------------------------------------------------------------------------------------------
bool Flag_Precheck( const int lv, const int PID, const int NoRefineBnd )
{

// pre-check 1. proper-nesting condition (all 26 siblings must exist)
// ===========================================================================================
   bool ProperNesting = true;

   for (int sib=0; sib<26; sib++)
   {
//    do not check if sibling[]<-1 to allow for refinement around boundaries
//    --> not considering OPT__NO_FLAG_NEAR_BOUNDARY yet
      if ( amr->patch[0][lv][PID]->sibling[sib] == -1 )
      {
         ProperNesting = false;
         break;
      }
   }

// check further if refinement around boundaries is forbidden
   if ( OPT__NO_FLAG_NEAR_BOUNDARY && ProperNesting )
   {
      for (int d=0; d<3; d++)
      {
         const int CornerL = amr->patch[0][lv][PID]->corner[d];
         const int CornerR = CornerL + Mis_Cell2Scale( PS1, lv );

         if ( CornerL <= 0                + NoRefineBnd  ||
              CornerR >= amr->BoxScale[d] - NoRefineBnd    )
         {
            ProperNesting = false;
            break;
         }
      }
   }

   if ( ! ProperNesting )  return false;


// pre-check 2. particle flag
// ===========================================================================================
#  ifdef PARTICLE
   if ( OPT__FLAG_PAR_TARGET == PAR_FLAG_CAN  ||  OPT__FLAG_PAR_TARGET == PAR_FLAG_BOTH )
   {
      if (  ! Par_Flag_TargetParticle( lv, PID, PAR_FLAG_CAN )  )
         return false;
   }
#  endif


// pre-check 3. whether this patch is within the regions allowed to be refined
// ===========================================================================================
   if ( OPT__FLAG_REGION )
   {
      if ( Flag_Region_Ptr == NULL )   Aux_Error( ERROR_INFO, "Flag_Region_Ptr == NULL for OPT__FLAG_REGION !!\n" );

      bool InRefineRegion = false;

      for (int k=0; k<PS1 && !InRefineRegion; k++)
      for (int j=0; j<PS1 && !InRefineRegion; j++)
      for (int i=0; i<PS1 && !InRefineRegion; i++)
         InRefineRegion = Flag_Region_Ptr( i, j, k, lv, PID ); // use = instead of |= as the loop exits immediately once InRefineRegion is true

      if ( ! InRefineRegion )    return false;
   } // if ( OPT__FLAG_REGION )


// pre-check 4. maximum angular resolution
// ===========================================================================================
   if ( OPT__FLAG_ANGULAR )
   {
      bool PassAngCheck = false;

      for (int k=0; k<PS1 && !PassAngCheck; k++)
      for (int j=0; j<PS1 && !PassAngCheck; j++)
      for (int i=0; i<PS1 && !PassAngCheck; i++)
         PassAngCheck = Check_Angular_Max( i, j, k, lv, PID, FLAG_ANGULAR_CEN_X, FLAG_ANGULAR_CEN_Y, FLAG_ANGULAR_CEN_Z,
                                           FlagTable_Angular[lv][0], FlagTable_Angular[lv][2] );   // also use = instead of |= here

      if ( ! PassAngCheck )   return false;
   } // if ( OPT__FLAG_ANGULAR )


// pass all pre-checks
   return true;

} // FUNCTION : Flag_Precheck
