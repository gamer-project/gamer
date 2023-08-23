#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_BlastWave
// Description :  User-defined flag criteria
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_BlastWave( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dr[3]  = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
   const double Radius = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );

   bool Flag = Radius < Threshold[0];


   return Flag;

} // FUNCTION : Flag_BlastWave
