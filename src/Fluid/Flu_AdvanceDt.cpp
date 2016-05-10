#include "Copyright.h"
#include "GAMER.h"
#include <climits>




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AdvanceDt
// Description :  Advance the fluid attributes by the flux gradients
//
// Note        :  a. Invoke the function "InvokeSolver"
//                b. Currently the updated data can only be stored in the different sandglass from the 
//                   input fluid data
//
// Parameter   :  lv             : Targeted refinement level 
//                TimeNew        : Targeted physical time to reach
//                TimeOld        : Physical time before update
//                                 --> This function updates physical time from TimeOld to TimeNew
//                dt             : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//                SaveSg         : Sandglass to store the updated data 
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void Flu_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt, const int SaveSg,
                    const bool OverlapMPI, const bool Overlap_Sync )
{

   InvokeSolver( FLUID_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg, NULL_INT, OverlapMPI, Overlap_Sync );

   if ( OPT__FIXUP_FLUX )  Buf_ResetBufferFlux( lv );

#  ifndef GRAVITY
   if ( OPT__RESET_FLUID ) Flu_ResetByUser( lv, SaveSg, TimeNew );
#  endif

} // FUNCTION : Flu_AdvanceDt
