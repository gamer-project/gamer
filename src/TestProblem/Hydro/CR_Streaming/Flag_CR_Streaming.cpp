#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_CR_Streaming
// Description :  Flag cells within a spherical shell around the box center
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER" with "OPT__FLAG_USER_NUM 2"
//                3. Depends only on position --> static refinement hierarchy
//                4. Set Threshold[1] <= Threshold[0] on a given level to disable refinement there
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold loaded from "Input__Flag_User"
//                            --> Threshold[0] = inner shell radius
//                                Threshold[1] = outer shell radius
//
// Return      :  "true"  if Threshold[0] <= r < Threshold[1]
//                "false" otherwise
//-------------------------------------------------------------------------------------------------------
bool Flag_CR_Streaming( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double dr[3]  = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
   const double Radius = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );

   const double R_in   = Threshold[0];
   const double R_out  = Threshold[1];

   if ( R_out <= R_in )    return false;

   bool Flag = ( Radius >= R_in  &&  Radius < R_out );


   return Flag;

} // FUNCTION : Flag_CR_Streaming
