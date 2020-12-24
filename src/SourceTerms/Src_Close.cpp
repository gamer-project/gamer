#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Close
// Description :  Closing step for the source terms
//
// Note        :  1. Store all fluid variables for now
//                   --> Should remove unused fields in the future
//                2. Do not store B field
//                3. Use the same array for input and output
//                4. Always update the latest FluSg data
//
// Parameter   :  lv            : Target refinement level
//                h_Flu_Array_S : Host array storing the updated fluid data
//                NPG           : Number of patch groups updated at a time
//                PID0_List     : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void Src_Close( const int lv, const real h_Flu_Array_S[][NCOMP_TOTAL][ CUBE(PS1) ],
                const int NPG, const int *PID0_List )
{

#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         const int N   = 8*TID + LocalID;

//       update all fluid variables for now
         memcpy( amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[0][0][0], h_Flu_Array_S[N][0],
                 NCOMP_TOTAL*CUBE(PS1)*sizeof(real) );
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Src_Close
