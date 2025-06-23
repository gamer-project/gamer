#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_AdvanceDt
// Description :  Update the internal energy by the various cooling and heating mechanisms in Grackle
//
// Note        :  1. Invoke InvokeSolver()
//                2. Invoked by EvolveLevel()
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                TimeOld      : Physical time before update
//                               --> This function updates physical time from TimeOld to TimeNew
//                dt           : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//                SaveSg       : Sandglass to store the updated data
//                OverlapMPI   : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync : true  --> Advance the patches which cannot be overlapped with MPI communication
//                               false --> Advance the patches which can    be overlapped with MPI communication
//                               (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void Grackle_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt, const int SaveSg,
                        const bool OverlapMPI, const bool Overlap_Sync )
{
#  ifdef COMOVING
// convert dt from the comoving time interval to the physical time interval
// --> see Equation (15) of Schive, Tsai, & Chiueh (2010)
   const double dt_grackle = dt * SQR(TimeOld);
#  else
   const double dt_grackle = dt;
#  endif

   InvokeSolver( GRACKLE_SOLVER, lv, TimeNew, TimeOld, dt_grackle, NULL_REAL, SaveSg, NULL_INT, NULL_INT, OverlapMPI, Overlap_Sync );

} // FUNCTION : Grackle_AdvanceDt



#endif // #ifdef SUPPORT_GRACKLE
