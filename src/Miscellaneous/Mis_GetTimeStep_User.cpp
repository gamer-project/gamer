#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_User
// Description :  User-defined criteria to estimate the evolution time-step and physical time interval
//
// Note        :  1. Physical coordinates : dTime == dt
//                   Comoving coordinates : dTime == dt*(Hubble parameter)*(scale factor)^3 == delta(scale factor)
//                2. Invoked by "Mis_GetTimeStep_Check" using the function pointer "Mis_GetTimeStep_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  dt       : Time interval to advance solution
//                dTime    : Time interval to update physical time
//                dt_dTime : dt/dTime (== 1.0 if COMOVING is off)
//
// Return      :  dt, dTime
//-------------------------------------------------------------------------------------------------------
void Mis_GetTimeStep_User( double &dt, double &dTime, const double dt_dTime )
{

// put your favorite time-step criteria here
// ##########################################################################################################

// Example 1 : set upper limit for the time interval to advance solution (per sub-step)
   /*
   double dt_user = 1.e-2;

// return 2*dt for the individual time-step since at the base level each step actually includes two sub-steps
#  ifdef INDIVIDUAL_TIMESTEP
   dt_user *= 2.0;
#  endif

   dt    = dt_user;
   dTime = dt / dt_dTime;
   */


// Example 2 : set upper limit for the time interval to update the scale factor in COMOVING
   /*
   double dTime_user = 2.e-5;

   dTime = dTime_user;
   dt    = dTime * dt_dTime;
   */

// ##########################################################################################################

} // FUNCTION : Mis_GetTimeStep_User

