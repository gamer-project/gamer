#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_AllocateFluxArray_Buffer
// Description :  Allocate flux arrays for the coarse-grid buffer patches (at level lv) adjacent to the
//                coarse-fine boundaries
//
// Parameter   :  lv : Target coarse-grid level
//-------------------------------------------------------------------------------------------------------
void Flu_AllocateFluxArray_Buffer( const int lv )
{

// check
   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


   const int MirrorSib[6] = { 1, 0, 3, 2, 5, 4 };
   int PID, SonPID, SibPID, Table[4];

   for (int s=0; s<6; s++)
   {
      for (int t=0; t<4; t++)    Table[t] = TABLE_03(s,t);

#     pragma omp parallel for private( PID, SonPID, SibPID ) schedule( runtime )
      for (int PID0=amr->NPatchComma[lv][s+1]; PID0<amr->NPatchComma[lv][s+2]; PID0+=8)
      for (int t=0; t<4; t++)
      {
         PID    = PID0 + Table[t];
         SonPID = amr->patch[0][lv][PID]->son;

         if ( SonPID == -1 )
         {
            SibPID = amr->patch[0][lv][PID]->sibling[ MirrorSib[s] ];

//          SibPID, if exists, must be a real patch
#           ifdef GAMER_DEBUG
            if ( SibPID < -1 )   Aux_Error( ERROR_INFO, "lv %d, SibPID (%d) < -1 !!\n", lv, SibPID );
#           endif

            if ( SibPID != -1 )
            if ( amr->patch[0][lv][SibPID]->son != -1 )
               amr->patch[0][lv][PID]->fnew( MirrorSib[s], AUTO_REDUCE_DT );
         }
      }
   } // for (int s=0; s<6; s++)

} // FUNCTION : Flu_AllocateFluxArray_Buffer



#endif // #ifndef SERIAL
