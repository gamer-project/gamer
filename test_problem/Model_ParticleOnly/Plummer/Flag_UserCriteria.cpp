#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_UserCriteria
// Description :  Check if the element (i,j,k) of the input data satisfies the user-defined flag criteria
//
// Note        :  Users can put their favorite flag criteria in this function
//
// Parameter   :  i,j,k       : Indices of the targeted element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the targeted patch
//                PID         : ID of the targeted patch
//                Threshold   : User-provided threshold for the flag operation, which is loaded from the 
//                              file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_UserCriteria( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{

   bool Flag = false;

   const int NParMin = 200;

   if ( amr->patch[0][lv][PID]->son == -1 )     Flag = amr->patch[0][lv][PID]->NPar > NParMin;
   else                                         Flag = Par_CountParticleInDescendant( lv, PID ) > NParMin;

   return Flag;

} // FUNCTION : Flag_UserCriteria

