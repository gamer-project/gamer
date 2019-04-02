#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_Prepare_Pot
// Description :  Prepare the input array "h_Pot_Array_P_Out" for the Gravity solver
//
// Note        :  Invoke Prepare_PatchData()
//
// Parameter   :  lv                : Target refinement level
//                PrepTime          : Target physical time to prepare the coarse-grid data
//                h_Pot_Array_P_Out : Host array to store the prepared data
//                NPG               : Number of patch groups prepared at a time
//                PID0_List         : List recording the patch indices with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Gra_Prepare_Pot( const int lv, const double PrepTime, real h_Pot_Array_P_Out[][GRA_NXT][GRA_NXT][GRA_NXT],
                      const int NPG, const int *PID0_List )
{

   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinDens_No        = -1.0;
   const real MinPres_No        = -1.0;

   Prepare_PatchData( lv, PrepTime, &h_Pot_Array_P_Out[0][0][0][0], NULL, GRA_GHOST_SIZE, NPG, PID0_List, _POTE, _NONE,
                      OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, (GRA_GHOST_SIZE==0)?NSIDE_00:NSIDE_06, IntPhase_No,
                      OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

} // FUNCTION : Gra_Prepare_Pot



#endif // #ifdef GRAVITY
