#include "GAMER.h"
#include "CUFLU.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Hybrid
// Description :  Estimate the evolution time-step from the CFL condition of the Hamilton-Jacobi solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. CFL constant purely empirical, a value of C_CFL = 0.4 seems to work for the 3rd-order RK phase schem
//                   a fourth-order RK method allows a value of up to C_CFL = 0.8, the second-order scheme requires a value of C_CFL = 0.05
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Hybrid( const int lv )
{

   const double dh = amr->dh[lv];
   double dt;

   dt = ELBDM_ETA*SQR(dh);

   dt *= DT__HYBRID;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Hybrid



#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )
