#include "Copyright.h"
#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_Close
// Description :  Copy the specific internal energy updated by the CPU/GPU Grackle solver back to the
//                patch pointers
//
// Note        :  1. Use SaveSg to determine where to store the data
//                   --> Currently it's set to the same Sg as the fluid data when calling
//                       Grackle_AdvanceDt() in EvolveLevel()
//
// Parameter   :  lv          : Target refinement level
//                SaveSg      : Sandglass to store the updated data
//                h_Che_Array : Host array storing the updated data
//                NPG         : Number of patch groups to store the updated data
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Grackle_Close( const int lv, const int SaveSg, const real h_Che_Array[][CHE_NPREP][ CUBE(PS1) ],
                    const int NPG, const int *PID0_List )
{

   const int Idx_Dens  = 0;
   const int Idx_sEint = 1;
   const int Idx_Ek    = 2;

   int N, PID, PID0;


#  pragma omp parallel for private( N, PID, PID0 ) schedule( runtime )
   for (int TID=0; TID<NPG; TID++)
   {
      PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;
         N   = 8*TID + LocalID;

         for (int t=0; t<CUBE(PS1); t++)
            *( amr->patch[SaveSg][lv][PID]->fluid[ENGY][0][0] + t )
               = h_Che_Array[N][Idx_sEint][t]*h_Che_Array[N][Idx_Dens][t] + h_Che_Array[N][Idx_Ek][t];
      }
   }

} // FUNCTION : Grackle_Close



#endif // #ifdef SUPPORT_GRACKLE
