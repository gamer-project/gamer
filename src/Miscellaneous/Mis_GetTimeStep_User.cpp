#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
double Mis_GetTimeStep_User_Template( const int lv, const double dTime_dt );

// this function pointer must be set by a test problem initializer
double (*Mis_GetTimeStep_User_Ptr)( const int lv, const double dTime_dt ) = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_User_Template
// Description :  Template of user-defined criteria to estimate the evolution time-step
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep() using the function pointer "Mis_GetTimeStep_User_Ptr",
//                   which must be set by a test problem initializer
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_User_Template( const int lv, const double dTime_dt )
{

// put your favorite time-step criteria here
// ##########################################################################################################
   double dt_user = HUGE_NUMBER;

// Example 1 : set upper limit for the time interval to advance solution (per sub-step)
   /*
   dt_user = 1.0e-2;
   */


// Example 2 : set upper limit for the time interval to update the scale factor in COMOVING
   /*
   double dTime_user = 1.0e-5;
   dt_user = dTime_user / dTime_dt;
   */
// ##########################################################################################################

   return dt_user;

} // FUNCTION : Mis_GetTimeStep_User_Template

