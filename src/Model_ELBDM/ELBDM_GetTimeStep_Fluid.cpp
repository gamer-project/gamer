#include "GAMER.h"
#include "CUFLU.h"

#if ( MODEL == ELBDM )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_GetTimeStep_Fluid
// Description :  Estimate the evolution time-step from the ELBDM kinematic energy solver
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Time-step is set to restrict the 1D k-space rotation angle to be "DT__FLUID*2*pi"
//
// Parameter   :  lv : Target refinement level
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double ELBDM_GetTimeStep_Fluid( const int lv )
{

   const double dh = amr->dh[lv];
   double dt;

   dt = 4.0/M_PI*ELBDM_ETA*SQR(dh);
/*
#  ifdef GRAVITY
   dt = 4.0/3.0/M_PI*ELBDM_ETA*SQR(dh);      // 3D k-space rotation angle
#  else
   dt = 0.5*sqrt(3.0)*ELBDM_ETA*SQR(dh);     // CFL condition
#  endif
*/

   dt *= (Step==0) ? DT__FLUID_INIT : DT__FLUID;

   return dt;

} // FUNCTION : ELBDM_GetTimeStep_Fluid



#endif // #if ( MODEL == ELBDM )
