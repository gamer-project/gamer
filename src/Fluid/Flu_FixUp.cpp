#include "GAMER.h"
#include "CUFLU.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_FixUp
// Description :  Correct the coarse-grid data with the fine-grid data
//
// Note        :  1. There are three types of correction
//                   (1) Flux              --> Flu_FixUp_Flux()
//                   (2) Spatial averaging --> Flu_FixUp_Restrict()
//                   (3) Electric field    --> MHD_FixUp_Electric()
//                2. Invoked by EvolveLevel()
//
// Parameter   :  lv : Target coarse level
//-------------------------------------------------------------------------------------------------------
void Flu_FixUp( const int lv )
{

// 1. use the fine-grid fluxes across the coarse-fine boundaries to correct the coarse-grid data
   if ( OPT__FIXUP_FLUX )  Flu_FixUp_Flux( lv );


// 2. use the average data at level "lv+1" to correct the data at level "lv"
//    --> we do not correct the potential data
   if ( OPT__FIXUP_RESTRICT )
   {
#     ifdef MHD
      const int SonMagSg = amr->MagSg[lv+1];
      const int  FaMagSg = amr->MagSg[lv  ];
#     else
      const int SonMagSg = NULL_INT;
      const int  FaMagSg = NULL_INT;
#     endif

      Flu_FixUp_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], NULL_INT, NULL_INT, SonMagSg, FaMagSg, _TOTAL, _MAG );
   }

// 3. use the fine-grid electric field on the coarse-fine boundaries to correct the coarse-grid magnetic field
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )    MHD_FixUp_Electric( lv );
#  endif

} // FUNCTION : Flu_FixUp
