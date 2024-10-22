#include "GAMER.h"



// flag for checking whether the tcool field is initialized
bool IsInit_tcool[NLEVEL] = { false };


//-------------------------------------------------------------------------------------------------------
// Function    :  Src_AdvanceDt
// Description :  Add various local source terms
//
// Note        :  1. Invoke InvokeSolver()
//                2. Invoked by EvolveLevel()
//                3. Grackle library is treated separately
//                4. Invoke Src_WorkBeforeMajorFunc()
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                TimeOld      : Physical time before update
//                               --> This function updates physical time from TimeOld to TimeNew
//                dt           : Time interval to advance solution
//                SaveSg_Flu   : Sandglass to store the updated fluid data
//                SaveSg_Mag   : Sandglass to store the updated B field (NOT SUPPORTED YET)
//                OverlapMPI   : true --> Overlap MPI time with CPU/GPU computation (NOT SUPPORTED YET)
//                Overlap_Sync : true  --> Advance the patches which cannot be overlapped with MPI communication
//                               false --> Advance the patches which can    be overlapped with MPI communication
//                               (useful only if "OverlapMPI == true")
//
// Return      : fluid[] in all patches
//-------------------------------------------------------------------------------------------------------
void Src_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                    const int SaveSg_Flu, const int SaveSg_Mag, const bool OverlapMPI, const bool Overlap_Sync )
{

// nothing to do if no source term is activated
   if ( ! SrcTerms.Any )  return;

// work before calling the major source-term function
   Src_WorkBeforeMajorFunc( lv, TimeNew, TimeOld, dt );

// major source-term function
   InvokeSolver( SRC_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, SaveSg_Mag, NULL_INT,
                 OverlapMPI, Overlap_Sync );

   if ( SrcTerms.ExactCooling )   IsInit_tcool[lv] = true;

} // FUNCTION : Src_AdvanceDt
