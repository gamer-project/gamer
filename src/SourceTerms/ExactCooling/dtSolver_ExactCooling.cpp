#include "GAMER.h"




#ifdef EXACT_COOLING
//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_ExactCooling
// Description :  Estimate the evolution time-step constrained by the ExactCooling source term
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
double Mis_GetTimeStep_ExactCooling( const int lv, const double dTime_dt )
{

   if ( !SrcTerms.ExactCooling )   return HUGE_NUMBER;


// allocate memory for per-thread arrays
#  ifdef OPENMP
   const int NT = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int NT = 1;
#  endif

   double  dt_EC     = HUGE_NUMBER;
   double *OMP_dt_EC;


   return dt_EC;

} // FUNCTION : Mis_GetTimeStep_ExactCooling
#endif // #ifdef EXACT_COOLING
