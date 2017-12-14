#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_Prepare_Flu
// Description :  Fill up the input array "h_Flu_Array_T" with fluid variables for estimating the
//                evolution time-step
//
// Note        :  1. Always prepare the latest FluSg data
//                2. Prepare NCOMP_FLUID variables
//                3. Use patches instead of patch groups as the basic unit
//                4. No ghost zones
//
// Parameter   :  lv            : Target refinement level
//                h_Flu_Array_T : Host array to store the prepared fluid data
//                NPG           : Number of patch groups prepared at a time
//                PID0_List     : List recording the target patch indicies with LocalID==0
//-------------------------------------------------------------------------------------------------------
void dt_Prepare_Flu( const int lv, real h_Flu_Array_T[][NCOMP_FLUID][ CUBE(PS1) ], const int NPG, const int *PID0_List )
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

//       exclude passive scalars
         for (int v=0; v<NCOMP_FLUID; v++)
         for (int t=0; t<CUBE(PS1); t++)
            h_Flu_Array_T[N][v][t] = *( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][0][0] + t );
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : dt_Prepare_Flu
