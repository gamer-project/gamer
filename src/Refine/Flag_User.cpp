#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static bool Flag_User_Template( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );

// this function pointer must be set by a test problem initializer
bool (*Flag_User_Ptr)( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_User_Template
// Description :  Template of user-defined flag criteria
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
bool Flag_User_Template( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   /*
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const real (*Rho )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];  // density
   const real (*MomX)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX];  // momentum x
   const real (*MomY)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY];  // momentum y
   const real (*MomZ)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ];  // momentum z
   const real (*Egy )[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY];  // total energy
#  ifdef GRAVITY
   const real (*Pot )[PS1][PS1] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;          // potential
#  endif
   */

   bool Flag = false;


// put your favorite flag criteria here
// ##########################################################################################################

   /*
// Example 1 : flag if the velocity exceeds the given threshold
   const real Vel = sqrt(  ( MomX[k][j][i]*MomX[k][j][i] + MomY[k][j][i]*MomY[k][j][i] +
                             MomZ[k][j][i]*MomZ[k][j][i] )  ) / Rho[k][j][i];
   Flag = Vel > Threshold;
   */


   /*
// Example 2 : flag if the grid is within a sphere with the radius eqaul to the input "Threshold" and the origin
//             in the center of the simulation box
   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dr[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double Radius    = sqrt( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );

   Flag = Radius < Threshold[0];
   */
// ##########################################################################################################


   return Flag;

} // FUNCTION : Flag_User_Template
