#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
extern int      Merger_Coll_NumBHs;
extern double   R_acc;
extern double (*CM_ClusterCen)[3];
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_ClusterMerger
// Description :  Flag cells for refinement for the cluster merger test problem
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by Init_TestProb_Hydro_ClusterMerger()
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
bool Flag_ClusterMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh        = amr->dh[lv];
   const double Pos[3]    = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double FlagNCell = Threshold[0]; // flag cells within "FlagRFac" times the accretion radius R_acc, and
   const double FlagRFac  = Threshold[1]; // if R_acc is not resolved with "FlagNCell" cells

   bool Flag = false;

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      if (  DIST_SQR_3D( Pos, CM_ClusterCen[c] ) <= SQR( FlagRFac*R_acc )  &&  R_acc/dh <= FlagNCell  )
      {
         Flag = true;
         return Flag;
      } // if ( R_SQR <= SQR(FlagRFac*R_acc)  &&  R_acc/dh <= FlagNCell )
   } // for (int c=0; c<Merger_Coll_NumBHs; c++)

   return Flag;

} // FUNCTION : Flag_ClusterMerger
