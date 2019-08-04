#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_LB_EnsureBFieldConsistencyAfterRestrict
// Description :  Ensure consistency of B field between leaf buffer patches and their sibling real/buffer
//                non-leaf patches after applying the restrict operation
//
// Note        :  1. Only necessary on the coarse-fine boundaries
//                   --> Use non-leaf sibling real/buffer patches to correct leaf buffer patches
//                       --> Non-leaf sibling **buffer** patches should have received the restricted data
//                           in advance when invoking LB_GetBufferData(DATA_AFTER_FIXUP) in EvolveLevel()
//                       --> Non-leaf sibling **real** patches should have received the restricted data
//                           in advance when invoking either Flu_FixUp_Restrict() (if father patches are
//                           also real patches) or LB_GetBufferData(DATA_RESTRICT) (if father patches are
//                           buffer patches) in EvolveLevel()
//                   --> Note that we must still correct leaf buffer patches even if the fathers of their
//                       sibling non-leaf patches are **real** patches
//                       --> Although the consistency in this case had been ensured in Flu_FixUp_Restrict() when
//                           copying data from father patches to father-sibling patches, it is broken again when
//                           we invoke another LB_GetBufferData(DATA_GENERAL) within LB_GetBufferData(DATA_RESTRICT)
//                           with zero ghost zone
//                2. Invoked by EvolveLevel()
//
// Parameter   :  lv : AMR level
//
// Return      :  Longitudinal B field on the coarse-fine interfaces of leaf buffer patches
//-------------------------------------------------------------------------------------------------------
void MHD_LB_EnsureBFieldConsistencyAfterRestrict( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );
#  endif


// nothing to do on the top level
   if ( lv == TOP_LEVEL )     return;


   const int MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
   const int MagSg        = amr->MagSg[lv];

// iterate over all buffer patches to be corrected
   for (int DesBufPID=amr->NPatchComma[lv][1]; DesBufPID<amr->NPatchComma[lv][3]; DesBufPID++)
   {
//    for the **destination** patches (i.e., patches to be corrected), skip the **non-leaf** buffer patches
      if ( amr->patch[MagSg][lv][DesBufPID]->magnetic != NULL  &&  amr->patch[0][lv][DesBufPID]->son == -1 )
      {
         for (int s=0; s<6; s++)
         {
            const int SrcBufPID = amr->patch[0][lv][DesBufPID]->sibling[s];

//          for the **source** patches (i.e., patches to correct others), skip the **leaf** patches
//          --> note that SrcBufPID can be either real or buffer patches
            if ( SrcBufPID >= 0  &&  amr->patch[MagSg][lv][SrcBufPID]->magnetic != NULL  &&
                 amr->patch[0][lv][SrcBufPID]->son != -1 )
               MHD_CopyPatchInterfaceBField( lv, SrcBufPID, MirrorSib[s], MagSg );
         }
      }
   } // for (int DesBufPID=NReal; DesBufPID<NTotal; DesBufPID++)

} // FUNCTION : MHD_LB_EnsureBFieldConsistencyAfterRestrict



#endif // #if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )
