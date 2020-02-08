#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_Prepare_Flu
// Description :  Fill up h_Flu_Array_T[] with fluid variables and h_Mag_Array_T[] with B field for
//                estimating the evolution time-step
//
// Note        :  1. Always prepare the latest FluSg and MagSg data
//                2. Prepare NCOMP_FLUID fluid variables and NCOMP_MAG B field components
//                3. Use patches instead of patch groups as the basic unit
//                4. No ghost zones
//
// Parameter   :  lv            : Target refinement level
//                h_Flu_Array_T : Host array to store the prepared fluid data
//                h_Mag_Array_T : Host array to store the prepared B field data
//                NPG           : Number of patch groups prepared at a time
//                PID0_List     : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void dt_Prepare_Flu( const int lv, real h_Flu_Array_T[][NCOMP_FLUID][ CUBE(PS1) ],
                     real h_Mag_Array_T[][NCOMP_MAG][ PS1P1*SQR(PS1) ], const int NPG, const int *PID0_List )
{

#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         const int N   = 8*TID + LocalID;

//       fluid variables (excluding passive scalars)
         memcpy( h_Flu_Array_T[N][0], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[0][0][0],
                 NCOMP_FLUID*CUBE(PS1)*sizeof(real) );

//       B field
#        ifdef MHD
         memcpy( h_Mag_Array_T[N][0], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[0],
                 NCOMP_MAG*PS1P1*SQR(PS1)*sizeof(real) );
#        endif
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : dt_Prepare_Flu
