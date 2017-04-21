#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Close
// Description :  Copy the momentum and energy density stored in the h_Flu_Array_G array back into the
//                patch pointers
//
// Note        :  The updated data will be stored in the same Sg as the current fluid data
//
// Parameter   :  lv             : Target refinement level
//                SaveSg         : Sandglass to store the updated data
//                h_Flu_Array_G  : Host array storing the updated fluid variables
//                h_DE_Array_G   : Host array storing the dual-energy status
//                NPG            : Number of patch groups to store the updated data
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Close( const int lv, const int SaveSg, const real h_Flu_Array_G[][GRA_NIN][PS1][PS1][PS1],
                const char h_DE_Array_G[][PS1][PS1][PS1], const int NPG, const int *PID0_List )
{

   int N, PID, PID0;

#  pragma omp parallel for private( N, PID, PID0 ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

#        if ( MODEL == HYDRO  ||  MODEL == MHD )
//       density field is sent in and out but NOT updated in the hydro/MHD gravity solver
         for (int v=1; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            amr->patch[SaveSg][lv][PID]->fluid[v][k][j][i] = h_Flu_Array_G[N][v][k][j][i];

//       dual-energy status (which is always stored in Sg=0)
#        ifdef DUAL_ENERGY
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            amr->patch[0][lv][PID]->de_status[k][j][i] = h_DE_Array_G[N][k][j][i];
#        endif

#        elif ( MODEL == ELBDM )
//       density field is NOT sent in and out in the ELBDM gravity solver
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            amr->patch[SaveSg][lv][PID]->fluid[v+1][k][j][i] = h_Flu_Array_G[N][v][k][j][i];

#        else
#        error : unsupported MODEL !!
#        endif // MODEL
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Gra_Close



#endif // #ifdef GRAVITY
