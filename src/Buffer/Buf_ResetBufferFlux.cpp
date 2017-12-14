#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_ResetBufferFlux
// Description :  Reset all fluxes in the buffer patches as zero
//
// Note        :  1. Invoked by Flu_AdvanceDt()
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Buf_ResetBufferFlux( const int lv )
{

// check
   if ( !amr->WithFlux )
   {
      Aux_Message( stderr, "WARNING : invoking %s is useless since no flux is required !!\n", __FUNCTION__ );
      return;
   }


   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;

#  pragma omp parallel for private( FluxPtr ) schedule( runtime )
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][27]; PID++)
   for (int s=0; s<6; s++)
   {
      FluxPtr = amr->patch[0][lv][PID]->flux[s];

      if ( FluxPtr != NULL )
      {
         for(int v=0; v<NFLUX_TOTAL; v++)
         for(int m=0; m<PATCH_SIZE; m++)
         for(int n=0; n<PATCH_SIZE; n++)
            FluxPtr[v][m][n] = 0.0;
      }
   }

} // FUNCTION : Buf_ResetBufferFlux



#endif // #ifndef SERIAL
