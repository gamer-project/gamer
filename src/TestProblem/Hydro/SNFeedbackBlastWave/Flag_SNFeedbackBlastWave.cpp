#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined FEEDBACK )



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_SNFeedbackBlastWave
// Description :  Flag cells for refinement for the SNFeedbackBlastWave test problem
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by Init_TestProb_Hydro_SNFeedbackBlastWave()
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the targeted element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the targeted patch
//                PID       : ID of the targeted patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_SNFeedbackBlastWave( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh = amr->dh[lv];
   const double dx = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh - amr->BoxCenter[0];
   const double dy = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh - amr->BoxCenter[1];
   const double dz = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh - amr->BoxCenter[2];

   bool Flag = false;


// flag cells within the central box (excluding the upper-right corner)
   const double RefineHalfWidth = Threshold[0]*dh;
   const bool   CentralCells    = ( fabs(dx) < RefineHalfWidth  &&  fabs(dy) < RefineHalfWidth  &&  fabs(dz) < RefineHalfWidth );
   const bool   UpperRight      = (      dx  >               0  &&       dy  >               0  &&       dz  >               0 );

   Flag = ( CentralCells  &&  !UpperRight );

   return Flag;

} // FUNCTION : Flag_SNFeedbackBlastWave



#endif // #if ( MODEL == HYDRO  &&  defined FEEDBACK )
