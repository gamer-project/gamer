#include "Copyright.h"
#include "GAMER.h"

extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AdvanceDt
// Description :  Advance the fluid attributes by the flux gradients
//
// Note        :  a. Invoke the function "InvokeSolver"
//                b. Currently the updated data can only be stored in the different sandglass from the
//                   input fluid data
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
void Flu_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt, const int SaveSg,
                    const bool OverlapMPI, const bool Overlap_Sync )
{

   InvokeSolver( FLUID_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg, NULL_INT, OverlapMPI, Overlap_Sync );

   if ( OPT__FIXUP_FLUX )  Buf_ResetBufferFlux( lv );

// call Flu_ResetByUser_API_Ptr() here only if both GRAVITY and GRACKLE are disabled
#  ifdef GRAVITY
   if ( false )
#  endif
#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_MODE == GRACKLE_MODE_NONE )
#  endif
   if ( OPT__RESET_FLUID  &&  Flu_ResetByUser_API_Ptr != NULL )   Flu_ResetByUser_API_Ptr( lv, SaveSg, TimeNew );

} // FUNCTION : Flu_AdvanceDt
