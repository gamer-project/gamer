#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AllocateFluxArray
// Description :  Allocate flux arrays for the coarse-grid patches (at level lv ) adjacent to the 
//                coarse-fine boundaries (including the buffer patches)
//
// Parameter   :  lv : Coarse-grid level
//-------------------------------------------------------------------------------------------------------
void Flu_AllocateFluxArray( const int lv )
{

// check
   if ( !amr->WithFlux )    
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


// deallocate the flux arrays allocated previously
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][7]; PID++)   amr->patch[0][lv][PID]->fdelete();


// allocate flux arrays for the real patches
   int SibPID;

   if ( amr->NPatchComma[lv+1][7] != 0 )
   {
#     pragma omp parallel for private( SibPID ) schedule( runtime )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         if ( amr->patch[0][lv][PID]->son == -1 )
         {
            for (int s=0; s<6; s++)
            {
               if (  ( SibPID = amr->patch[0][lv][PID]->sibling[s] ) >= 0  )
                  if ( amr->patch[0][lv][SibPID]->son != -1 )  amr->patch[0][lv][PID]->fnew( s, AUTO_REDUCE_DT );
            }
         }
      }
   }


// allocate flux arrays for the buffer patches
   if ( amr->NPatchComma[lv+1][7] != 0 )  Flu_AllocateFluxArray_Buffer( lv );

   
// get the PIDs for sending/receiving fluxes to/from neighboring ranks
   Buf_RecordExchangeFluxPatchID( lv );

} // Flu_AllocateFluxArray
