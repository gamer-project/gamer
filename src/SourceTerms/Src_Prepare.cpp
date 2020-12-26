#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Src_Prepare
// Description :  Prepare the input arrays h_Flu_Array_S_In[], h_Mag_Array_S_In[], and h_Corner_Array_S[]
//                for source terms
//
// Note        :  1. Always prepare the latest FluSg and MagSg data
//                2. Prepare FLU_NIN_S fluid variables and NCOMP_MAG B field components
//                   --> FLU_NIN_S is fixed to NCOMP_TOTAL for now
//                   --> Should remove unused fields in the future
//                3. Use patches instead of patch groups as the basic unit
//                4. No ghost zones
//                   --> Should support ghost zones in the future
//                5. Corner coordinates are defined as the central coordinates of the first cell located
//                   at the bottom left corner
//                   --> Excluding ghost zones
//                   --> Implementation is the same as Gra_Prepare_Corner()
//
// Parameter   :  lv               : Target refinement level
//                h_Flu_Array_S_In : Host array to store the prepared fluid   data
//                h_Mag_Array_S_In : Host array to store the prepared B field data
//                h_Corner_Array_S : Host array to store the prepared corner  data
//                NPG              : Number of patch groups prepared at a time
//                PID0_List        : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void Src_Prepare( const int lv,
                  real h_Flu_Array_S_In[][FLU_NIN_S][ CUBE(PS1)      ],
                  real h_Mag_Array_S_In[][NCOMP_MAG][ PS1P1*SQR(PS1) ],
                  double h_Corner_Array_S[][3],
                  const int NPG, const int *PID0_List )
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

//       fluid variables (include all fields for now)
         memcpy( h_Flu_Array_S_In[N][0], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[0][0][0],
                 FLU_NIN_S*CUBE(PS1)*sizeof(real) );

//       B field
#        ifdef MHD
         memcpy( h_Mag_Array_S_In[N][0], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[0],
                 NCOMP_MAG*PS1P1*SQR(PS1)*sizeof(real) );
#        endif

//       corner coordinates
         for (int d=0; d<3; d++)    h_Corner_Array_S[N][d] = amr->patch[0][lv][PID]->EdgeL[d] + dh_half;
      }
   } // for (int TID=0; TID<NPG; TID++)

} // FUNCTION : Src_Prepare
