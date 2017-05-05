#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_UserCriteria
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  Users can put their favorite flag criteria in this function
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the
//                              file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_UserCriteria( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
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

} // FUNCTION : Flag_UserCriteria

