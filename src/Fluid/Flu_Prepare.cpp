#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Prepare
// Description :  Prepare input arrays for the fluid solver
//
// Note        :  Invoke Prepare_PatchData()
//
// Parameter   :  lv                   : Target refinement level
//                PrepTime             : Target physical time to prepare the coarse-grid data
//                h_Flu_Array_F_In     : Host array to store the prepared fluid data
//                h_Mag_Array_F_In     : Host array to store the prepared B field (for MHD onlhy)
//                h_Pot_Array_USG_F    : Host array to store the prepared potential data (for UNSPLIT_GRAVITY only)
//                h_Corner_Array_USG_F : Host array to store the prepared corner data (for UNSPLIT_GRAVITY only)
//                NPG                  : Number of patch groups to be prepared at a time
//                PID0_List            : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Flu_Prepare( const int lv, const double PrepTime,
                  real h_Flu_Array_F_In[][FLU_NIN][ CUBE(FLU_NXT) ],
                  real h_Mag_Array_F_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                  real h_Pot_Array_USG_F[][ CUBE(USG_NXT_F) ],
                  double h_Corner_Array_F[][3], const int NPG, const int *PID0_List )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if (  ( OPT__GRAVITY_TYPE == GRAVITY_SELF || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&
         ( h_Pot_Array_USG_F == NULL )  )
      Aux_Error( ERROR_INFO, "h_Pot_Array_USG_F == NULL !!\n" );

   if (  ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL || OPT__GRAVITY_TYPE == GRAVITY_BOTH || OPT__EXTERNAL_POT )  &&
         ( h_Corner_Array_F == NULL )  )
      Aux_Error( ERROR_INFO, "h_Corner_Array_F == NULL !!\n" );
#  endif
#  endif


#  if ( MODEL != HYDRO )
   const double MIN_DENS            = -1.0;  // set to an arbitrarily negative value to disable it
   const double MIN_PRES            = -1.0;  // ...
#  endif
#  ifndef MHD
   const int    OPT__MAG_INT_SCHEME = INT_NONE;
#  endif
   const bool   IntPhase_No         = false;
   const real   MinDens_No          = -1.0;
   const real   MinPres_No          = -1.0;
   const bool   DE_Consistency_Yes  = true;
   const bool   DE_Consistency_No   = false;
   const bool   DE_Consistency      = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
   const real   MinDens             = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;


// prepare the fluid array
#  if ( MODEL == ELBDM )
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], NULL,
                      FLU_GHOST_SIZE, NPG, PID0_List, _REAL|_IMAG|_PASSIVE, _NONE,
                      OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE,
                      OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, DE_Consistency_No );
#  else
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In[0][0], h_Mag_Array_F_In[0][0],
                      FLU_GHOST_SIZE, NPG, PID0_List, _TOTAL, _MAG,
                      OPT__FLU_INT_SCHEME, OPT__MAG_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens,    MinPres_No, DE_Consistency );
#  endif

#  ifdef UNSPLIT_GRAVITY
// prepare the potential array
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   Prepare_PatchData( lv, PrepTime, h_Pot_Array_USG_F[0], NULL,
                      USG_GHOST_SIZE_F, NPG, PID0_List, _POTE, _NONE,
                      OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

// prepare the corner array
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
   {
      const double dh_half = 0.5*amr->dh[lv];

      int PID0;

//#     pragma omp parallel for private( PID0 ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

         for (int d=0; d<3; d++)    h_Corner_Array_F[TID][d] = amr->patch[0][lv][PID0]->EdgeL[d] + dh_half;
      } // for (int TID=0; TID<NPG; TID++)
   }
#  endif // #ifdef UNSPLIT_GRAVITY

} // FUNCTION : Flu_Prepare
