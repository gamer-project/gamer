#include "GAMER.h"
#include "CUFLU.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp
// Description :  Correct the coarse-grid data with the fine-grid data
//
// Note        :  1. There are three types of correction
//                   (1) Spatial averaging --> Flu_FixUp_Restrict()
//                   (2) Electric field    --> MHD_FixUp_Electric()
//                   (3) Flux              --> Flu_FixUp_Flux()
//                2. Invoked by EvolveLevel()
//
// Parameter   :  lv : Target coarse level
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp( const int lv )
{

// 1. use the average data at level "lv+1" to correct the data at level "lv"
//    --> we do not correct the potential data
   if ( OPT__FIXUP_RESTRICT )    Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], amr->MagSg[lv+1], amr->MagSg[lv],
                                                     NULL_INT, NULL_INT, _TOTAL, _MAG );


// 2. use the fine-grid electric field on the coarse-fine boundaries to correct the coarse-grid magnetic field
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )    MHD_FixUp_Electric( lv );
#  endif


// 3. use the fine-grid fluxes across the coarse-fine boundaries to correct the coarse-grid data
//    --> apply this correction AFTER all other checks since it will check negative pressure as well
//        (which requires the coarse-grid B field updated by Flu_FixUp_Restrict() and MHD_FixUp_Electric())
   if ( OPT__FIXUP_FLUX )  Flu_FixUp_Flux( lv );

} // FUNCTION : Flu_FixUp
