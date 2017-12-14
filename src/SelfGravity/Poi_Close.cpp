#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Close
// Description :  Copy the potential stored in the h_Pot_Array_P_Out array back into the patch pointers
//
// Note        :  The potential will be stored in the same Sg as the current fluid data
//
// Parameter   :  lv                : Target refinement level
//                SaveSg            : Sandglass to store the updated data
//                h_Pot_Array_P_Out : Host array storing the updated potential
//                NPG               : Number of patch groups to store the updated data
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Close( const int lv, const int SaveSg, const real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                const int NPG, const int *PID0_List )
{

   int ii, jj, kk, N, PID, PID0;

#  pragma omp parallel for private( ii, jj, kk, N, PID, PID0 ) schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

         for (int k=0; k<PATCH_SIZE; k++)    {  kk = k + GRA_GHOST_SIZE;
         for (int j=0; j<PATCH_SIZE; j++)    {  jj = j + GRA_GHOST_SIZE;
         for (int i=0; i<PATCH_SIZE; i++)    {  ii = i + GRA_GHOST_SIZE;

            amr->patch[SaveSg][lv][PID]->pot[k][j][i] = h_Pot_Array_P_Out[N][kk][jj][ii];

         }}}

#        ifdef STORE_POT_GHOST
         memcpy( amr->patch[SaveSg][lv][PID]->pot_ext, h_Pot_Array_P_Out[N], CUBE(GRA_NXT)*sizeof(real) );
#        endif
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Poi_Close



#endif // #ifdef GRAVITY
