#include "GAMER.h"

#ifdef SUPPORT_GRACKLE

extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const double TTime );




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

   InvokeSolver( GRACKLE_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg, NULL_INT, OverlapMPI, Overlap_Sync );

// always call Flu_ResetByUser_API_Ptr() here
// --> when Grackle is enabled, we do not invoke Flu_ResetByUser_API_Ptr() in either Flu_AdvanceDt or Gra_AdvanceDt
// --> we want to invoke Flu_ResetByUser_API_Ptr() before calling Buf_GetBufferData() to reduce the MPI communication
   if ( OPT__RESET_FLUID  &&  Flu_ResetByUser_API_Ptr != NULL )   Flu_ResetByUser_API_Ptr( lv, SaveSg, TimeNew );

} // FUNCTION : Grackle_AdvanceDt



#endif // #ifdef SUPPORT_GRACKLE
