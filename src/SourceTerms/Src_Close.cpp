#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Close
// Description :  Closing step for the source terms
//
// Note        :  1. Store FLU_NOUT_S fluid variables
//                   --> FLU_NOUT_S is fixed to NCOMP_TOTAL for now
//                   --> Should remove unused fields in the future
//                2. Do not store B field
//                3. Use the same array for input and output
//
// Parameter   :  lv                : Target refinement level
//                SaveSg_Flu        : Sandglass to store the updated fluid data
//                h_Flu_Array_S_Out : Host array storing the updated fluid data
//                NPG               : Number of patch groups updated at a time
//                PID0_List         : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void Src_Close( const int lv, const int SaveSg_Flu, const real h_Flu_Array_S_Out[][FLU_NOUT_S][ CUBE(PS1) ],
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
         memcpy( amr->patch[SaveSg_Flu][lv][PID]->fluid[0][0][0], h_Flu_Array_S_Out[N][0],
                 FLU_NOUT_S*CUBE(PS1)*sizeof(real) );
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Src_Close
