#include "GAMER.h"


// status of the fluid solver used by AUTO_REDUCE_DT
int FluStatus_ThisRank;


// defined in Flu_ManageFixUpTempArray.cpp
void Flu_SwapFixUpTempArray( const int lv );
void Flu_InitFixUpTempArray( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AdvanceDt
// Description :  Advance the fluid attributes by the flux gradients
//
// Note        :  1. Invoke InvokeSolver()
//                2. Currently the updated data can only be stored in the different sandglass from the
//                   input data
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                TimeOld      : Physical time before update
//                               --> This function updates physical time from TimeOld to TimeNew
//                dt           : Time interval to advance solution (can be different from TimeNew-TimeOld in COMOVING)
//                SaveSg_Flu   : Sandglass to store the updated fluid data
//                SaveSg_Mag   : Sandglass to store the updated B field
//                OverlapMPI   : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync : true  --> Advance the patches which cannot be overlapped with MPI communication
//                               false --> Advance the patches which can    be overlapped with MPI communication
//                               (useful only if "OverlapMPI == true")
//
// Return      : GAMER_SUCCESS / GAMER_FAILED
//               --> Mainly used for the option "AUTO_REDUCE_DT"
//-------------------------------------------------------------------------------------------------------
int Flu_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const int SaveSg_Flu, const int SaveSg_Mag, const bool OverlapMPI, const bool Overlap_Sync )
{

// initialize flux_tmp[] (and electric_tmp[] in MHD) on the parent level for AUTO_REDUCE_DT
   if ( AUTO_REDUCE_DT  &&  lv != 0 )  Flu_InitFixUpTempArray( lv-1 );


// initialize patch->ele_corrected[] for correcting the coarse-grid electric field
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC  &&  lv != 0 )
   {
      const int FaLv = lv - 1;

      for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][19]; FaPID++)
      for (int e=0; e<12; e++)
         amr->patch[0][FaLv][FaPID]->ele_corrected[e] = false;
   }
#  endif


// invoke the fluid solver
   FluStatus_ThisRank = GAMER_SUCCESS;
#  if ( MODEL == ELBDM  &&  defined SUPPORT_FFTW )
   if ( lv == 0  &&  ELBDM_BASE_SPECTRAL )
   {
      CPU_ELBDMSolver_FFT( dt, TimeOld, SaveSg_Flu );
   }
   else
#  endif
   InvokeSolver( FLUID_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, SaveSg_Mag, NULL_INT, OverlapMPI, Overlap_Sync );


// collect the fluid solver status from all ranks (only necessary for AUTO_REDUCE_DT)
   int FluStatus_AllRank;

// the parenthesis enclosing MPI_Allreduce() is necessary for the serial mode
// note that AUTO_REDUCE_DT may deteriorate the performance of OPT__MINIMIZE_MPI_BARRIER because of the following extra MPI call
// --> use MPI_BAND instead of MPI_BOR since GAMER_SUCCESS=1
   if ( AUTO_REDUCE_DT )   { MPI_Allreduce( &FluStatus_ThisRank, &FluStatus_AllRank, 1, MPI_INT, MPI_BAND, MPI_COMM_WORLD ); }
   else                    FluStatus_AllRank = GAMER_SUCCESS;


// note that we always have FluStatus == GAMER_SUCCESS if AUTO_REDUCE_DT is disabled
   if ( FluStatus_AllRank == GAMER_SUCCESS )
   {
//    reset the fluxes and electric field in the buffer patches at lv as zeros
//    --> for accumulating the coarse-fine fluxes and electric field later when evolving lv+1
      if ( OPT__FIXUP_FLUX )  Buf_ResetBufferFlux( lv );

#     if ( defined MHD  &&  defined LOAD_BALANCE )
      if ( OPT__FIXUP_ELECTRIC )    MHD_LB_ResetBufferElectric( lv );
#     endif


//    swap the flux (and electric in MHD) pointers on the parent level if the fluid solver works successfully
      if ( AUTO_REDUCE_DT  &&  lv != 0 )  Flu_SwapFixUpTempArray( lv-1 );
   }


   return FluStatus_AllRank;

} // FUNCTION : Flu_AdvanceDt
