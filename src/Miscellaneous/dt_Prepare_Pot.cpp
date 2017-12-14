#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  dt_Prepare_Pot
// Description :  Fill up the input array "h_Pot_Array_T" with gravitational potential for estimating the
//                evolution time-step
//
// Note        :  1. Always prepare the latest PotSg data with GRA_GHOST_SIZE ghost zones
//                2. Use patches instead of patch groups as the basic unit
//                3. It does NOT include the external potential, which will be added in CPU/CUPOT_dtSolver_XXX()
//
// Parameter   :  lv            : Target refinement level
//                h_Pot_Array_T : Host array to store the prepared gravitational potential
//                NPG           : Number of patch groups prepared at a time
//                PID0_List     : List recording the target patch indicies with LocalID==0
//                PrepTime      : Target physical time
//-------------------------------------------------------------------------------------------------------
void dt_Prepare_Pot( const int lv, real h_Pot_Array_T[][ CUBE(GRA_NXT) ], const int NPG, const int *PID0_List,
                     const double PrepTime )
{

// nothing to do if self-gravity is disabled
   if ( OPT__GRAVITY_TYPE != GRAVITY_SELF  &&  OPT__GRAVITY_TYPE != GRAVITY_BOTH )  return;


   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;

   Prepare_PatchData( lv, PrepTime, &h_Pot_Array_T[0][0], GRA_GHOST_SIZE, NPG, PID0_List, _POTE,
                      OPT__GRA_INT_SCHEME, UNIT_PATCH, (GRA_GHOST_SIZE==0)?NSIDE_00:NSIDE_06, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

} // FUNCTION : dt_Prepare_Pot



#endif // #ifdef GRAVITY
