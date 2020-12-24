#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Prepare
// Description :  Fill up h_Flu_Array_S[] with fluid variables and h_Mag_Array_S[] with B field for
//                source terms
//
// Note        :  1. Always prepare the latest FluSg and MagSg data
//                2. Prepare all (i.e., NCOMP_TOTAL) fluid variables and NCOMP_MAG B field components
//                   --> Should remove unused fields in the future
//                3. Use patches instead of patch groups as the basic unit
//                4. No ghost zones
//                   --> Should support ghost zones in the future
//
// Parameter   :  lv            : Target refinement level
//                h_Flu_Array_S : Host array to store the prepared fluid data
//                h_Mag_Array_S : Host array to store the prepared B field data
//                NPG           : Number of patch groups prepared at a time
//                PID0_List     : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void Src_Prepare( const int lv,
                  real h_Flu_Array_S[][NCOMP_TOTAL][ CUBE(PS1)      ],
                  real h_Mag_Array_S[][NCOMP_MAG  ][ PS1P1*SQR(PS1) ],
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

//       fluid variables (include all fields for now)
         memcpy( h_Flu_Array_S[N][0], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[0][0][0],
                 NCOMP_TOTAL*CUBE(PS1)*sizeof(real) );

//       B field
#        ifdef MHD
         memcpy( h_Mag_Array_S[N][0], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[0],
                 NCOMP_MAG*PS1P1*SQR(PS1)*sizeof(real) );
#        endif
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Src_Prepare
