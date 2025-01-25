#include "GAMER.h"
#include "CUFLU.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Hybrid_CFL
// Description :  Estimate the evolution time-step from the CFL condition of the Hamilton-Jacobi solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. CFL constant purely empirical, a value of C_CFL = 0.4 seems to work for the 3rd-order RK phase scheme
//                   a fourth-order RK method allows a value of up to C_CFL = 0.8, the second-order scheme requires a value of C_CFL = 0.05
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Hybrid_CFL( const int lv )
{

   const double dh = amr->dh[lv];
   double dt;

   dt = ELBDM_ETA*SQR(dh);

   dt *= (Step==0) ? DT__HYBRID_CFL_INIT : DT__HYBRID_CFL;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Hybrid_CFL



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
