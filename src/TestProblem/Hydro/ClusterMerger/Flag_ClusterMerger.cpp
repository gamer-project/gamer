#include "GAMER.h"




#if ( MODEL == HYDRO  &&  defined GRAVITY )

static bool   FirstTime = true;
extern int    Merger_Coll_NumHalos;
extern double R_acc;                // the radius to compute the accretion rate
extern double ClusterCen[3][3];


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

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Flag = false;

// flag cells within the target radius, and if the radius is not resolved with a specific number (Threshold[0]) of cells
   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      double R_SQR = SQR(Pos[0]-ClusterCen[c][0]) + SQR(Pos[1]-ClusterCen[c][1]) + SQR(Pos[2]-ClusterCen[c][2]);
      if ( R_SQR <= SQR(25*R_acc)  &&  R_acc/dh <= Threshold[0] )
      {
         Flag = true;
         return Flag;
      } // if ( R_SQR <= SQR(25*R_acc)  &&  R_acc/dh <= Threshold[0] )
   } // for (int c=0; c<Merger_Coll_NumHalos; c++)

   if ( FirstTime )
   {
      const double dh_max = amr->dh[MAX_LEVEL];
      if ( R_acc/dh_max <= Threshold[0] )
         Aux_Message( stderr, "WARNING : MAX_LEVEL is less than the desired refinement level set in Input__Flag_User!! dh_max = %13.7e\n", dh_max );

      FirstTime = false;
   } // if ( FirstTime )

   return Flag;

} // FUNCTION : Flag_ClusterMerger



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
