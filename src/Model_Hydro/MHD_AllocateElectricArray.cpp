#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_AllocateElectricArray
// Description :  Allocate electric arrays for the coarse-grid real patches (on level lv) adjacent to the
//                coarse-fine boundaries
//
// Note        :  1. Do nothing on the top level
//                2. Only used in the serial mode
//
// Parameter   :  lv : Coarse-grid level
//-------------------------------------------------------------------------------------------------------
void MHD_AllocateElectricArray( const int lv )
{

// check
   if ( !amr->WithElectric )
      Aux_Error( ERROR_INFO, "amr->WithElectric is off !!\n" );

   if ( lv < 0  ||  lv > TOP_LEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter lv = %d (either < 0 or > %d) !!\n", lv, TOP_LEVEL );


// nothing to do on the top level
   if ( lv == TOP_LEVEL )  return;


// deallocate the electric arrays allocated previously
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][19]; PID++)  amr->patch[0][lv][PID]->edelete();


// nothing to do if there are no sibling-son patches in any target direction
   if ( amr->NPatchComma[lv+1][19] == 0 )    return;


// allocate electric arrays for the real patches
#  pragma omp parallel for schedule( runtime )
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      if ( amr->patch[0][lv][PID]->son == -1 )
      {
         for (int EdgeID=0; EdgeID<6; EdgeID++)
         {
            const int SibPID = amr->patch[0][lv][PID]->sibling[EdgeID];
            if ( SibPID >= 0  &&  amr->patch[0][lv][SibPID]->son != -1 )
               amr->patch[0][lv][PID]->enew( EdgeID, AUTO_REDUCE_DT );
         }

         for (int EdgeID=6; EdgeID<18; EdgeID++)
         {
            int SibID[3];
            TABLE_SiblingSharingSameEdge( EdgeID, SibID, NULL );

            for (int s=0; s<3; s++)
            {
               const int SibPID = amr->patch[0][lv][PID]->sibling[ SibID[s] ];
               if ( SibPID >= 0  &&  amr->patch[0][lv][SibPID]->son != -1 )
               {
                  amr->patch[0][lv][PID]->enew( EdgeID, AUTO_REDUCE_DT );
                  break;
               }
            }
         }
      }
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

} // MHD_AllocateElectricArray



#endif // #if ( MODEL == HYDRO  &&  defined MHD )
