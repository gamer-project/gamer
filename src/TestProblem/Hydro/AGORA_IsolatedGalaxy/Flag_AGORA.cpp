#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_AGORA
// Description :  Flag cells for refinement for the AGORA isolated galaxy simulation
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by "Init_TestProb_Hydro_AGORA_IsolatedGalaxy()"
//                   to replace "Flag_User()"
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the targeted element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the targeted patch
//                PID         : ID of the targeted patch
//                Threshold   : Useless here
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_AGORA( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

// flag cells within the target region [Threshold ... BoxSize-Threshold]
   const double EdgeL = Threshold;
   const double EdgeR = amr->BoxSize[0]-Threshold;    // here we have assumed a cubic box

   bool Flag;

   if (  Pos[0] >= EdgeL  &&  Pos[0] < EdgeR  &&
         Pos[1] >= EdgeL  &&  Pos[1] < EdgeR  &&
         Pos[2] >= EdgeL  &&  Pos[2] < EdgeR     )
      Flag = true;

   else
      Flag = false;


   return Flag;

} // FUNCTION : Flag_AGORA

