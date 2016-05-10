#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Prepare_Flu
// Description :  Fill up the input array "h_Flu_Array_G" with fluid variables for the Gravity solver  
//
// Parameter   :  lv             : Targeted refinement level
//                h_Flu_Array_G  : Host array to store the prepared data
//                NPG            : Number of patch groups prepared at a time
//                PID0_List      : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Prepare_Flu( const int lv, real h_Flu_Array_G[][GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE],
                      const int NPG, const int *PID0_List )
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

//       density field is useless in the ELBDM gravity solver
#        if ( MODEL == ELBDM )
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            h_Flu_Array_G[N][v][k][j][i] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v+1][k][j][i];

#        else
         for (int v=0; v<GRA_NIN; v++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
            h_Flu_Array_G[N][v][k][j][i] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];
#        endif 
      }
   }

} // FUNCTION : Gra_Prepare_Flu



#endif // #ifdef GRAVITY
