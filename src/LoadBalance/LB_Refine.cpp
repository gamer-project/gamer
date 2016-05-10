#include "Copyright.h"
#include "GAMER.h"

#ifdef LOAD_BALANCE



void LB_Refine_GetNewRealPatchList( const int FaLv, int &NNew_Home, int *&NewPID_Home, int &NNew_Away, 
                                    ulong *&NewCr1D_Away, real *&NewCData_Away, int &NDel_Home, int *&DelPID_Home, 
                                    int &NDel_Away, ulong *&DelCr1D_Away );
void LB_Refine_AllocateNewPatch( const int FaLv, int NNew_Home, int *NewPID_Home, int NNew_Away, 
                                 ulong *NewCr1D_Away, real *NewCData_Away, int NDel_Home, int *DelPID_Home, 
                                 int NDel_Away, ulong *DelCr1D_Away );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine
// Description :  Construct patches at FaLv+1 according to the flagging results at level FaLv 
//
// Note        :  1. This function will also construct buffer patches at both FaLv and FaLv+1
//                2. Data of all sibling-buffer patches must be prepared in advance for creating new 
//                   patches at FaLv+1 by spatial interpolation
//
// Parameter   :  FaLv  : Targeted refinement level to be refined
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


// 1. construct LB_CutPoint for the newly-created level
// ==========================================================================================
//###NOTE : here we have assumed that the newly-created son patches will have LB_Idx ~ 8*(father LB_Idx)
//          --> hold for Hilbert curve method, but may not hold for other methods
   if ( NPatchTotal[SonLv] == 0 )
      for (int r=0; r<MPI_NRank+1; r++)   amr->LB->CutPoint[SonLv][r] = amr->LB->CutPoint[FaLv][r]*8;


// 2. get the unsorted lists of father patches at FaLv to allocate/deallocate son patches at FaLv+1
// ==========================================================================================
   int    NNew_Home, NDel_Home, NNew_Away, NDel_Away;  // Home/Away : for targeted patches at home/not at home
   int   *NewPID_Home=NULL, *DelPID_Home=NULL;
   ulong *NewCr1D_Away=NULL, *DelCr1D_Away=NULL;
   real  *NewCData_Away=NULL;

   LB_Refine_GetNewRealPatchList( FaLv, NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCData_Away, 
                                  NDel_Home, DelPID_Home, NDel_Away, DelCr1D_Away );


// 3. allocate/deallocate son patches at FaLv+1
// ==========================================================================================
   LB_Refine_AllocateNewPatch( FaLv, NNew_Home, NewPID_Home, NNew_Away, NewCr1D_Away, NewCData_Away, 
                               NDel_Home, DelPID_Home, NDel_Away, DelCr1D_Away );


// 4. construct the MPI send and recv data list 
// ==========================================================================================
// 4.1 list for exchanging hydro and potential data in buffer patches (also construct the "SibDiff" lists)
   LB_RecordExchangeDataPatchID(  FaLv, true );
   LB_RecordExchangeDataPatchID( SonLv, true );

// 4.2 list for exchanging restricted hydro data
#  ifndef GAMER_DEBUG
   if ( OPT__FIXUP_RESTRICT )    
#  endif
   {
      LB_RecordExchangeRestrictDataPatchID(  FaLv );
      LB_RecordExchangeRestrictDataPatchID( SonLv );
   }

// 4.3 list for exchanging hydro fluxes (also allocate flux arrays)
   if ( amr->WithFlux )        
   {
      LB_AllocateFluxArray(  FaLv ); 
      LB_AllocateFluxArray( SonLv ); 
   }

// 4.4 list for exchanging hydro data after the fix-up operation
   if ( OPT__FIXUP_RESTRICT  ||  OPT__FIXUP_FLUX )
   {
      LB_RecordExchangeFixUpDataPatchID(  FaLv );
      LB_RecordExchangeFixUpDataPatchID( SonLv );
   }

// 4.5 list for overlapping MPI time with CPU/GPU computation
   if ( OPT__OVERLAP_MPI )
   {
      LB_RecordOverlapMPIPatchID(  FaLv );
      LB_RecordOverlapMPIPatchID( SonLv );
   }


// 5. miscellaneous 
// ==========================================================================================
// 5.1 reset LB_CutPoint to the default values (-1) if SonLv has been totally removed 
   if ( NPatchTotal[SonLv] == 0 )
      for (int r=0; r<MPI_NRank+1; r++)   amr->LB->CutPoint[SonLv][r] = -1;

// 5.2 free memory
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

} // FUNCTION : LB_Refine



#endif // #ifdef LOAD_BALANCE
