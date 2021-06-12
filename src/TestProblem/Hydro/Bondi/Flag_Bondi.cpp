#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined GRAVITY )



extern double Bondi_RefineRadius0;
extern bool   Bondi_HalfMaxLvRefR;




//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Bondi
// Description :  Flag cells for refinement for the Bondi accretion test problem
//
// Note        :  1. Linked to the function pointer "Flag_User_Ptr" by Init_TestProb_Hydro_Bondi()
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
bool Flag_Bondi( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh     = amr->dh[lv];
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Flag = false;


// flag cells within the target radius
   const double dr[3]   = { Pos[0]-amr->BoxCenter[0], Pos[1]-amr->BoxCenter[1], Pos[2]-amr->BoxCenter[2] };
   const double Radius  = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
   const double RefineR = (Bondi_HalfMaxLvRefR && lv==MAX_LEVEL-1) ? ( Bondi_RefineRadius0 / double(1<<(lv+1)) )
                                                                   : ( Bondi_RefineRadius0 / double(1<<(lv  )) );

   Flag = Radius < RefineR;

   return Flag;

} // FUNCTION : Flag_Bondi



#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )
