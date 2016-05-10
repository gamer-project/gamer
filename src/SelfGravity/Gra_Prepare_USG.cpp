#include "Copyright.h"
#include "GAMER.h"

#if (defined GRAVITY  &&  defined UNSPLIT_GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Prepare_USG
// Description :  Prepare the input arrays "h_Pot_Array_USG_G" and "h_Flu_Array_USG_G" for the Gravity solver 
//                when UNSPLIT_GRAVITY is adopted
//
// Note        :  1. Invoke the function "Prepare_PatchData"
//                2. Prepare density and potential at the **previous** time-step at Lv=lv
//                   --> Use the PREP_OLD option in Prepare_PatchData
//                   --> Data at the **current** time-step should already be prepared by the original Gravity solver
//                3. Still need "PrepTime" to determine whether temporal interpolation is required for the
//                   **Lv=lv-1** data
//
// Parameter   :  lv                : Targeted refinement level
//                PrepTime          : Targeted physical time to prepare the coarse-grid data
//                h_Pot_Array_USG_G : Host array to store the prepared potential (size = USG_NXT_G^3)
//                h_Flu_Array_USG_G : Host array to store the prepared density   (size =       PS1^3)
//                NPG               : Number of patch groups prepared at a time
//                PID0_List         : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Prepare_USG( const int lv, const double PrepTime,
                      real h_Pot_Array_USG_G[][USG_NXT_G][USG_NXT_G][USG_NXT_G], 
                      real h_Flu_Array_USG_G[][GRA_NIN-1][PS1][PS1][PS1], const int NPG, const int *PID0_List )
{

   const OptFluBC_t *FluBC_None = NULL;
   const bool IntPhase_No       = false;
   const bool GetTotDens_Yes    = true;
   const bool GetTotDens_No     = false;

// prepare potential
   if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
      if ( lv == 0 )
      Prepare_PatchData( lv, PrepTime, &h_Pot_Array_USG_G[0][0][0][0], USG_GHOST_SIZE, NPG, PID0_List, _POTE,
                         OPT__GRA_INT_SCHEME, UNIT_PATCH, NSIDE_06, IntPhase_No, FluBC_None,  OPT__BC_POT, GetTotDens_No );
      else
      Prepare_PatchData( lv, PrepTime, &h_Pot_Array_USG_G[0][0][0][0], USG_GHOST_SIZE, NPG, PID0_List, _POTE,
                         OPT__GRA_INT_SCHEME, UNIT_PATCH, NSIDE_06, IntPhase_No, FluBC_None,  OPT__BC_POT, GetTotDens_No );
   }

// prepare density + momentum
      Prepare_PatchData( lv, PrepTime,  h_Flu_Array_USG_G[0][0][0][0], 0,              NPG, PID0_List, _DENS|_MOMX|_MOMY|_MOMZ,
                         INT_NONE,            UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE, GetTotDens_No );

} // FUNCTION : Gra_Prepare_USG



#endif // #if (defined GRAVITY  &&  defined UNSPLIT_GRAVITY )
