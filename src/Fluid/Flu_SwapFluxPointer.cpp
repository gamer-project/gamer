#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_SwapFluxPointer
// Description :  Swap the flux and flux_tmp pointers for all patches in the target level
//
// Note        :  1. Work on both real and buffer patches
//                2. Work for the option "AUTO_REDUCE_DT", for which Flu_Close()->CorrectFlux()
//                   will store the updated (corrected) coarse-grid fluxes in the temporary array flux_tmp[]
//                   instead of flux[]
//                3. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flu_SwapFluxPointer( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect lv (%d) !!\n", lv );
#  endif


// loop over the flux arrays in all patches
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      for (int s=0; s<6; s++)
      {
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL )
         {
            Aux_SwapPointer( (void**)&amr->patch[0][lv][PID]->flux    [s],
                             (void**)&amr->patch[0][lv][PID]->flux_tmp[s] );
         }
      }
   }

} // FUNCTION : Flu_SwapFluxPointer



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_InitTempFlux
// Description :  Initialize the flux_tmp array
//
// Note        :  1. Copy flux to flux_tmp
//                2. Work for the option "AUTO_REDUCE_DT", for which Flu_SwapFluxPointer() will swap the
//                   flux and flux_tmp pointers
//                3. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Flu_InitTempFlux( const int lv )
{

// check
#  ifdef GAMER_DEBUG
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect lv (%d) !!\n", lv );
#  endif


// loop over the flux_tmp arrays in all patches
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      for (int s=0; s<6; s++)
      {
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL )
            memcpy( amr->patch[0][lv][PID]->flux_tmp[s], amr->patch[0][lv][PID]->flux[s], SQR(PS1)*NFLUX_TOTAL*sizeof(real) );
      }
   }

} // FUNCTION : Flu_InitTempFlux
