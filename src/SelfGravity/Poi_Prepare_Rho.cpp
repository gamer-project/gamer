#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_Prepare_Rho
// Description :  Prepare h_Rho_Array_P[] for the Poisson solver
//
// Note        :  1. Invoke the function "Prepare_PatchData"
//                2. Minimum density threshold (MIN_DENS) is applied
//
// Parameter   :  lv            : Target refinement level
//                PrepTime      : Target physical time to prepare the coarse-grid data
//                h_Rho_Array_P : Host array to store the prepared data
//                NPG           : Number of patch groups to be prepared at a time
//                PID0_List     : List recording the patch indicies with LocalID==0 to be udpated
//-------------------------------------------------------------------------------------------------------
void Poi_Prepare_Rho( const int lv, const double PrepTime, real h_Rho_Array_P[][RHO_NXT][RHO_NXT][RHO_NXT],
                      const int NPG, const int *PID0_List )
{

   const bool IntPhase_No       = false;
   const bool DE_Consistency_No = false;
   const real MinPres_No        = -1.0;

   Prepare_PatchData( lv, PrepTime, &h_Rho_Array_P[0][0][0][0], RHO_GHOST_SIZE, NPG, PID0_List, _TOTAL_DENS,
                      OPT__RHO_INT_SCHEME, UNIT_PATCH, NSIDE_26, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                      MIN_DENS, MinPres_No, DE_Consistency_No );

// subtract the background density to be consistent with the periodic (and/or comoving) base-level FFT solver
#  ifdef GAMER_DEBUG
   if ( OPT__BC_POT == BC_POT_PERIODIC  &&  AveDensity_Init <= 0.0 )
      Aux_Error( ERROR_INFO, "AveDensity_Init (%14.7e) hasn't been set properly !!\n", AveDensity_Init );
#  endif

#  ifdef COMOVING
   const bool Comoving = true;
#  else
   const bool Comoving = false;
#  endif

   real RhoSubtract = NULL_REAL;
   if      ( OPT__BC_POT == BC_POT_PERIODIC )   RhoSubtract = AveDensity_Init;
   else if ( Comoving )                         RhoSubtract = (real)1.0;

   if ( OPT__BC_POT == BC_POT_PERIODIC  ||  Comoving )
   {
#     pragma omp parallel for schedule( static )
      for (int TID=0; TID<8*NPG; TID++)
      for (int k=0; k<RHO_NXT; k++)
      for (int j=0; j<RHO_NXT; j++)
      for (int i=0; i<RHO_NXT; i++)
         h_Rho_Array_P[TID][k][j][i] -= RhoSubtract;
   }

} // FUNCTION : Poi_Prepare_Rho



#endif // #ifdef GRAVITY
