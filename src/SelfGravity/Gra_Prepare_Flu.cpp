#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Prepare_Flu
// Description :  Fill up the input array "h_Flu_Array_G" with fluid variables for the Gravity solver
//                --> When DUAL_ENERGY is on, this function also prepares the dual-energy status array "h_DE_Array_G"
//
// Note        :  1. This function always prepare the latest FluSg data
//
// Parameter   :  lv             : Target refinement level
//                h_Flu_Array_G  : Host array to store the prepared data
//                h_DE_Array_G   : Host array to store the dual-energy status
//                NPG            : Number of patch groups prepared at a time
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Prepare_Flu( const int lv, real h_Flu_Array_G[][GRA_NIN][PS1][PS1][PS1], char h_DE_Array_G[][PS1][PS1][PS1],
                      const int NPG, const int *PID0_List )
{

   int N, PID, PID0;

#  pragma omp parallel for private( N, PID, PID0 ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

#        if ( MODEL == HYDRO  ||  MODEL == MHD )
//       all active fields (including density) are sent into the hydro/MHD gravity solver
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
            h_Flu_Array_G[N][v][k][j][i] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

//       dual-energy status (which is always stored in Sg=0)
#        ifdef DUAL_ENERGY
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
            h_DE_Array_G[N][k][j][i] = amr->patch[0][lv][PID]->de_status[k][j][i];
#        endif

#        elif ( MODEL == ELBDM )
//       density field is useless in the ELBDM gravity solver
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PS1; k++)
         for (int j=0; j<PS1; j++)
         for (int i=0; i<PS1; i++)
            h_Flu_Array_G[N][v][k][j][i] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v+1][k][j][i];

#        else
#        error : unsupported MODEL !!
#        endif // MODEL
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Gra_Prepare_Flu



#endif // #ifdef GRAVITY
