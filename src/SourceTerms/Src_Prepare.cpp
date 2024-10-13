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
//                PrepTime         : Target physical time to prepare the data
//                h_Flu_Array_S_In : Host array to store the prepared fluid   data
//                h_Mag_Array_S_In : Host array to store the prepared B field data
//                h_Corner_Array_S : Host array to store the prepared corner  data
//                NPG              : Number of patch groups prepared at a time
//                PID0_List        : List recording the target patch indices with LocalID==0
//-------------------------------------------------------------------------------------------------------
void Src_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_S_In[][FLU_NIN_S][ CUBE(SRC_NXT)           ],
                  real h_Mag_Array_S_In[][NCOMP_MAG][ SRC_NXT_P1*SQR(SRC_NXT) ],
                  double h_Corner_Array_S[][3],
                  const int NPG, const int *PID0_List )
{

// check
// assuming we are always preparing the lastest data
#  ifdef GAMER_DEBUG
   if ( PrepTime != amr->FluSgTime[lv][ amr->FluSg[lv] ] )
      Aux_Error( ERROR_INFO, "PrepTime (%13.7e) != FluSgTime (%13.7e) !!\n",
                 PrepTime, amr->FluSgTime[lv][ amr->FluSg[lv] ] );

#  ifdef MHD
   if ( PrepTime != amr->MagSgTime[lv][ amr->MagSg[lv] ] )
      Aux_Error( ERROR_INFO, "PrepTime (%13.7e) != MagSgTime (%13.7e) !!\n",
                 PrepTime, amr->MagSgTime[lv][ amr->MagSg[lv] ] );
#  endif
#  endif // #ifdef GAMER_DEBUG


   const double dh_half = 0.5*amr->dh[lv];

#  pragma omp parallel for schedule( static )
   for (int TID=0; TID<NPG; TID++)
   {
      const int PID0 = PID0_List[TID];

      for (int LocalID=0; LocalID<8; LocalID++)
      {
         const int PID = PID0 + LocalID;
         const int N   = 8*TID + LocalID;

//       1. fast version for zero ghost zone
#        if ( SRC_GHOST_SIZE == 0 )
//       fluid variables (include all fields for now)
         memcpy( h_Flu_Array_S_In[N][0], amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[0][0][0],
                 FLU_NIN_S*CUBE(SRC_NXT)*sizeof(real) );

//       B field
#        ifdef MHD
         memcpy( h_Mag_Array_S_In[N][0], amr->patch[ amr->MagSg[lv] ][lv][PID]->magnetic[0],
                 NCOMP_MAG*SRC_NXT_P1*SQR(SRC_NXT)*sizeof(real) );
#        endif
#        endif // #if ( SRC_GHOST_SIZE == 0 )

//       corner coordinates
         for (int d=0; d<3; d++)    h_Corner_Array_S[N][d] = amr->patch[0][lv][PID]->EdgeL[d] + dh_half;
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int TID=0; TID<NPG; TID++)


// 2. slower version for non-zero ghost zones
#  if ( SRC_GHOST_SIZE > 0 )
#  ifndef MHD
   const int    OPT__MAG_INT_SCHEME = INT_NONE;
#  endif
   const bool   IntPhase_No         = false;
   const real   MinDens_No          = -1.0;
   const real   MinPres_No          = -1.0;
   const real   MinTemp_No          = -1.0;
   const real   MinEntr_No          = -1.0;
   const bool   DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? false : true;
   const real   MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;

// always prepare all fields for now (so FLU_NIN_S = NCOMP_TOTAL)
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_S_In[0][0], h_Mag_Array_S_In[0][0],
                      SRC_GHOST_SIZE, NPG, PID0_List, _TOTAL, _MAG,
                      OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCH, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency );
#  endif // #if ( SRC_GHOST_SIZE > 0 )

} // FUNCTION : Src_Prepare
