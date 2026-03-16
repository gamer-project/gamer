#include "GAMER.h"


extern FieldIdx_t RefineFieldIdx;

//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_ClusterMerger
// Description :  Flag cells for refinement for the galaxy cluster merger simulation
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by "Init_TestProb_ClusterMerger()"
//                   to replace "Flag_User()"
//                2. Please turn on the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//                Threshold   : Useless here
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_ClusterMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];           // density
   const real (*Scal)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[RefineFieldIdx]; // passive scalar

   // flag if the mass fraction of the scalar exceeds the given threshold
   const real Frac = Scal[k][j][i] / Rho[k][j][i];
   bool Flag = Frac > Threshold[0];

   return Flag;

} // FUNCTION : Flag_ClusterMerger

