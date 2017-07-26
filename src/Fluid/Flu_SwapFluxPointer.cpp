#include "Copyright.h"
#include "GAMER.h"

#define FLUX_UNUSED     NULL_REAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_SwapFluxPointer
// Description :  Swap the flux and flux_tmp pointers for all patches in the target level
//
// Note        :  1. Work on both real and buffer patches
//                2. Work for the option "OPT__AUTO_REDUCE_DT", for which Flu_Close()->CorrectFlux()
//                   will store the updated (corrected) coarse-grid fluxes in the temporary array flux_tmp[]
//                   instead of flux[]
//                3. Invoked by Flu_AdvanceDt()
//                4. Only swap the pointers if flux_tmp has been updated
//                   --> Check if the [0][0][0] element is still FLUX_UNUSED or not
//                   --> We must not swap the pointers if flux_tmp has not been updated by Flu_Close()->CorrectFlux()
//                       since it's values are incorrect and ill-defined
//                       --> This may happen, for example, if the coarse and fine patches adjacent to a
//                           coarse-fine boundary are at different MPI ranks
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
//       swap pointers only if flux_tmp[s][0][0][0] has been updated
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL  &&  amr->patch[0][lv][PID]->flux_tmp[s][0][0][0] != FLUX_UNUSED )
         {
            Aux_SwapPointer( (void**)&amr->patch[0][lv][PID]->flux    [s],
                             (void**)&amr->patch[0][lv][PID]->flux_tmp[s] );
         }
      }
   }

} // FUNCTION : Flu_SwapFluxPointer



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_InitTempFlux
// Description :  Initialize the flux_tmp array so that one can check later whether this array has been updated
//                or not
//
// Note        :  1. Set the [0][0][0] element as FLUX_UNUSED
//                2. Work for the option "OPT__AUTO_REDUCE_DT", for which Flu_SwapFluxPointer() will swap the
//                   flux and flux_tmp pointers if flux_tmp has been updated
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
//#  pragma omp parallel for private( tmp_ptr ) schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][27]; PID++)
   {
      for (int s=0; s<6; s++)
      {
         if ( amr->patch[0][lv][PID]->flux_tmp[s] != NULL )
            amr->patch[0][lv][PID]->flux_tmp[s][0][0][0] = FLUX_UNUSED;
      }
   }

} // FUNCTION : Flu_InitTempFlux
