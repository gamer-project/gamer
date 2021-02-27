#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Prepare_Corner
// Description :  Prepare h_Corner_Array_PGT[] for the Poisson, Gravity, and dt solvers
//                --> NOT just gravity solver (so the function name is a bit misleading)
//
// Note        :  1. For calculating the external acceleration and/or potential
//                2. Corner coordinates are defined as the central coordinates of the first cell located
//                   at the bottom left corner
//                   --> Excluding ghost zones
//
// Parameter   :  lv                 : Target refinement level
//                h_Corner_Array_PGT : Host array to store the prepared data
//                NPG                : Number of patch groups prepared at a time
//                PID0_List          : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Prepare_Corner( const int lv, double h_Corner_Array_PGT[][3], const int NPG, const int *PID0_List )
{

   const double dh_half = 0.5*amr->dh[lv];

#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         const int N   = 8*TID + LocalID;

         for (int d=0; d<3; d++)    h_Corner_Array_PGT[N][d] = amr->patch[0][lv][PID]->EdgeL[d] + dh_half;
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Gra_Prepare_Corner



#endif // #ifdef GRAVITY
