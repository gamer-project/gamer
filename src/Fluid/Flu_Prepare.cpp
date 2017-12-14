#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_Prepare
// Description :  Prepare the input array "Flu_Array_F_In" for the fluid solver
//
// Note        :  Invoke the function "Prepare_PatchData"
//
// Parameter   :  lv                   : Target refinement level
//                PrepTime             : Target physical time to prepare the coarse-grid data
//                h_Flu_Array_F_In     : Host array to store the prepared fluid data
//                h_Pot_Array_USG_F    : Host array to store the prepared potential data (for UNSPLIT_GRAVITY only)
//                h_Corner_Array_USG_F : Host array to store the prepared corner data (for UNSPLIT_GRAVITY only)
//                NPG                  : Number of patch groups to be prepared at a time
//                PID0_List            : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Flu_Prepare( const int lv, const double PrepTime, real h_Flu_Array_F_In[], real h_Pot_Array_USG_F[],
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


#  if ( MODEL != HYDRO  &&  MODEL != MHD )
   const double MIN_DENS           = -1.0;    // set to an arbitrarily negative value to disable it
   const double MIN_PRES           = -1.0;    // ...
#  endif

   const bool   IntPhase_No        = false;
   const real   MinDens_No         = -1.0;
   const real   MinPres_No         = -1.0;
   const bool   DE_Consistency_Yes = true;
   const bool   DE_Consistency_No  = false;
   const bool   DE_Consistency     = ( OPT__OPTIMIZE_AGGRESSIVE ) ? DE_Consistency_No : DE_Consistency_Yes;
   const real   MinDens            = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinDens_No : MIN_DENS;
   const real   MinPres            = ( OPT__OPTIMIZE_AGGRESSIVE ) ? MinPres_No : MIN_PRES;


// prepare the fluid array
#  if ( MODEL == ELBDM )
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In,  FLU_GHOST_SIZE, NPG, PID0_List, _REAL|_IMAG|_PASSIVE,
                      OPT__FLU_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, OPT__INT_PHASE,
                      OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, DE_Consistency_No );
#  else
   Prepare_PatchData( lv, PrepTime, h_Flu_Array_F_In,  FLU_GHOST_SIZE, NPG, PID0_List, _TOTAL,
                      OPT__FLU_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, BC_POT_NONE, MinDens, MinPres, DE_Consistency );
#  endif

#  ifdef UNSPLIT_GRAVITY
// prepare the potential array
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   Prepare_PatchData( lv, PrepTime, h_Pot_Array_USG_F, USG_GHOST_SIZE, NPG, PID0_List,
                      _POTE,         OPT__GRA_INT_SCHEME, UNIT_PATCHGROUP, NSIDE_26, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

// prepare the corner array
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
   {
      const double dh_half = 0.5*amr->dh[lv];

      int PID0;

//#     pragma omp parallel for private(  PID0 ) schedule( runtime )
      for (int TID=0; TID<NPG; TID++)
      {
         PID0 = PID0_List[TID];

         for (int d=0; d<3; d++)    h_Corner_Array_F[TID][d] = amr->patch[0][lv][PID0]->EdgeL[d] + dh_half;
      } // for (int TID=0; TID<NPG; TID++)
   }
#  endif // #ifdef UNSPLIT_GRAVITY

} // FUNCTION : Flu_Prepare
