#include "GAMER.h"

#ifdef LOAD_BALANCE



void LB_Refine_GetNewRealPatchList( const int FaLv, int &NNew_Home, int *&NewPID_Home, int &NNew_Away,
                                    ulong *&NewCr1D_Away, real *&NewCData_Away, int &NDel_Home, int *&DelPID_Home,
                                    int &NDel_Away, ulong *&DelCr1D_Away,
                                    int &RefineF2S_Send_NPatchTotal, int *&RefineF2S_Send_PIDList );
void LB_Refine_AllocateNewPatch( const int FaLv, int NNew_Home, int *NewPID_Home, int NNew_Away,
                                 ulong *NewCr1D_Away, real *NewCData_Away, int NDel_Home, int *DelPID_Home,
                                 int NDel_Away, ulong *DelCr1D_Away,
                                 int &RefineS2F_Send_NPatchTotal, int *&RefineS2F_Send_PIDList );
#ifdef PARTICLE
void Par_LB_Refine_SendParticle2Father( const int FaLv, const int RefineS2F_Send_NPatchTotal, int *RefineS2F_Send_PIDList );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine
// Description :  Construct patches at FaLv+1 according to the flagging results at level FaLv
//
// Note        :  1. This function will also construct buffer patches at both FaLv and FaLv+1
//                2. Data of all sibling-buffer patches must be prepared in advance for creating new
//                   patches at FaLv+1 by spatial interpolation
//
// Parameter   :  FaLv : Target refinement level to be refined
//-------------------------------------------------------------------------------------------------------
void LB_Refine( const int FaLv )
{

   const int SonLv = FaLv + 1;

// check
   if ( FaLv == NLEVEL-1 )
   {
      Aux_Message( stderr, "WARNING : function <%s> should NOT be applied to the finest level !!\n",
                   __FUNCTION__ );
      return;
   }

// check the synchronization
   if ( NPatchTotal[SonLv] != 0 )   Mis_CompareRealValue( Time[FaLv], Time[SonLv], __FUNCTION__, true );


// 0. initialize variables for exchanging particles
// ==========================================================================================
   int   RefineS2F_Send_NPatchTotal = 0;
   int  *RefineS2F_Send_PIDList     = NULL;
   int   RefineF2S_Send_NPatchTotal = 0;
   int  *RefineF2S_Send_PIDList     = NULL;

#  ifdef PARTICLE
   RefineS2F_Send_PIDList   = new int  [ amr->NPatchComma[FaLv][3] - amr->NPatchComma[FaLv][1] ];
   RefineF2S_Send_PIDList   = new int  [ amr->NPatchComma[FaLv][1] ];
#  endif


// 1. construct LB_CutPoint for the newly-created level
// ==========================================================================================
//###NOTE : here we have assumed that the newly-created son patches will have LB_Idx ~ 8*(father LB_Idx)
//          --> hold for Hilbert curve method, but may not hold for other methods
   if ( NPatchTotal[SonLv] == 0 )
      for (int r=0; r<MPI_NRank+1; r++)   amr->LB->CutPoint[SonLv][r] = amr->LB->CutPoint[FaLv][r]*8;


// 2. get the unsorted lists of father patches at FaLv to allocate/deallocate son patches at FaLv+1
// ==========================================================================================
   int    NNew_Home, NDel_Home, NNew_Away, NDel_Away;  // Home/Away : for target patches at home/not at home
   int   *NewPID_Home=NULL, *DelPID_Home=NULL;
   ulong *NewCr1D_Away=NULL, *DelCr1D_Away=NULL;
   real  *NewCData_Away=NULL;

   LB_Refine_GetNewRealPatchList( FaLv, NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCData_Away,
                                  NDel_Home, DelPID_Home, NDel_Away, DelCr1D_Away,
                                  RefineF2S_Send_NPatchTotal, RefineF2S_Send_PIDList );


// 3. allocate/deallocate son patches at FaLv+1
// ==========================================================================================
   LB_Refine_AllocateNewPatch( FaLv, NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCData_Away,
                               NDel_Home, DelPID_Home, NDel_Away, DelCr1D_Away,
                               RefineS2F_Send_NPatchTotal, RefineS2F_Send_PIDList );


// 4. construct the MPI send and recv data list
// ==========================================================================================
// 4.1 list for exchanging hydro and potential data in buffer patches (also construct the "SibDiff" lists)
   LB_RecordExchangeDataPatchID(  FaLv, true );
   LB_RecordExchangeDataPatchID( SonLv, true );

// 4.2 list for exchanging restricted hydro data
//     --> note that even when OPT__FIXUP_RESTRICT is off we still need to do data restriction in several places
//         (e.g., restart, and OPT__CORR_AFTER_ALL_SYNC)
//     --> for simplicity and sustainability, we always invoke LB_RecordExchangeRestrictDataPatchID()
   LB_RecordExchangeRestrictDataPatchID(  FaLv );
   LB_RecordExchangeRestrictDataPatchID( SonLv );

// 4.3 list for exchanging hydro fluxes (also allocate flux arrays)
   if ( amr->WithFlux )
   {
      LB_AllocateFluxArray(  FaLv );
      LB_AllocateFluxArray( SonLv );
   }

// 4.4 list for exchanging hydro data after the fix-up operation
//     --> for simplicity and sustainability, we always invoke LB_RecordExchangeFixUpDataPatchID()
//     --> see the comments 4.2 above
   LB_RecordExchangeFixUpDataPatchID(  FaLv );
   LB_RecordExchangeFixUpDataPatchID( SonLv );

// 4.5 list for overlapping MPI time with CPU/GPU computation
   if ( OPT__OVERLAP_MPI )
   {
      LB_RecordOverlapMPIPatchID(  FaLv );
      LB_RecordOverlapMPIPatchID( SonLv );
   }

// 4.6 list for exchanging particles
#  ifdef PARTICLE
   Par_LB_RecordExchangeParticlePatchID( SonLv );

   if ( SonLv < MAX_LEVEL )
   Par_LB_RecordExchangeParticlePatchID( SonLv+1 );

// only for reconstructing the amr->Par->B2R_Real/Buffer_NPatchTotal[FaLv][0] lists
   if ( FaLv >= 0 )
   Par_LB_RecordExchangeParticlePatchID( FaLv );
#  endif


// 5. send particles to leaf real patches
// ==========================================================================================
#  ifdef PARTICLE
// 5.1 send particles from buffer patches at FaLv to their corresponding real patches
   Par_LB_Refine_SendParticle2Father( FaLv, RefineS2F_Send_NPatchTotal, RefineS2F_Send_PIDList );

// 5.2 send particles from real patches at FaLv to their real son patches living abroad
   const bool TimingSendPar_No = false;
   Par_PassParticle2Son_MultiPatch( FaLv, PAR_PASS2SON_GENERAL, TimingSendPar_No,
                                    RefineF2S_Send_NPatchTotal, RefineF2S_Send_PIDList );
#  endif


// 6. miscellaneous
// ==========================================================================================
// 6.1 reset LB_CutPoint to the default values (-1) if SonLv has been totally removed
   if ( NPatchTotal[SonLv] == 0 )
      for (int r=0; r<MPI_NRank+1; r++)   amr->LB->CutPoint[SonLv][r] = -1;

// 6.2 free memory
   if ( NewPID_Home   == NULL  &&  NNew_Home != 0 )
      Aux_Error( ERROR_INFO, "%s has not been allocated !!\n", "NewPID_Home"   );

   if ( DelPID_Home   == NULL  &&  NDel_Home != 0 )
      Aux_Error( ERROR_INFO, "%s has not been allocated !!\n", "DelPID_Home"   );

   if ( NewCr1D_Away == NULL )
      Aux_Error( ERROR_INFO, "%s has not been allocated !!\n", "NewCr1D_Away" );

   if ( NewCData_Away == NULL )
      Aux_Error( ERROR_INFO, "%s has not been allocated !!\n", "NewCData_Away" );

   if ( DelCr1D_Away == NULL )
      Aux_Error( ERROR_INFO, "%s has not been allocated !!\n", "DelCr1D_Away" );

   free( NewPID_Home );
   free( DelPID_Home );
   delete [] NewCr1D_Away;
   delete [] NewCData_Away;
   delete [] DelCr1D_Away;

#  ifdef PARTICLE
   delete [] RefineS2F_Send_PIDList;
   delete [] RefineF2S_Send_PIDList;
#  endif

} // FUNCTION : LB_Refine



#endif // #ifdef LOAD_BALANCE
