#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_RecordExchangeFluxPatchID
// Description :  Record the information of patches for sending and receiving fluxes between neighbor ranks
//                in the variable "amr->ParaVar"
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void Buf_RecordExchangeFluxPatchID( const int lv )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL-1 )   Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );

   if ( !amr->WithFlux )
   {
      Aux_Message( stderr, "WARNING : invoking %s is useless since no flux is required !!\n", __FUNCTION__ );
      return;
   }


   const int MirrorSib[6] = { 1,0,3,2,5,4 };
   int PID, SibPID;
   real (*FluxPtr)[PATCH_SIZE][PATCH_SIZE] = NULL;

   for (int s=0; s<6; s++)
   {
//    initialize counters as zero
      amr->ParaVar->SendF_NList[lv][s] = 0;
      amr->ParaVar->RecvF_NList[lv][s] = 0;


//    deallocate memory
      if ( amr->ParaVar->SendF_IDList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->SendF_IDList[lv][s];
         amr->ParaVar->SendF_IDList[lv][s] = NULL;
      }

      if ( amr->ParaVar->RecvF_IDList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->RecvF_IDList[lv][s];
         amr->ParaVar->RecvF_IDList[lv][s] = NULL;
      }


//    nothing to do if there are no boundary patches at the target direction
      if ( amr->ParaVar->BounP_NList[lv][s] == 0 )    continue;


//    nothing to do if the target direction lies outside the simulation domain for the non-periodic B.C.
      if ( MPI_SibRank[s] < 0 )  continue;


//    allocate the maximum necessary memory
      amr->ParaVar->SendF_IDList [lv][s] = new int [ amr->ParaVar->BounP_NList[lv][s] ];
      amr->ParaVar->RecvF_IDList [lv][s] = new int [ amr->ParaVar->BounP_NList[lv][s] ];


//    a. set up the SendF_IDList
//    --------------------------
      for (int ID=0; ID<amr->ParaVar->BounP_NList[lv][s]; ID++)
      {
         PID = amr->ParaVar->BounP_IDList[lv][s][ID];

         if ( PID < amr->NPatchComma[lv][1] )
         {
            SibPID = amr->patch[0][lv][PID]->sibling[s];

            if ( SibPID >= 0 )
            {
               FluxPtr = amr->patch[0][lv][SibPID]->flux[ MirrorSib[s] ];

               if ( FluxPtr != NULL )
               {
                  amr->ParaVar->SendF_IDList[lv][s][ amr->ParaVar->SendF_NList[lv][s] ] = SibPID;
                  amr->ParaVar->SendF_NList [lv][s] ++;
               }
            }
         }
      }


//    b. set up the RecvF_IDList
//    --------------------------
      for (int ID=0; ID<amr->ParaVar->BounP_NList[lv][s]; ID++)
      {
         PID = amr->ParaVar->BounP_IDList[lv][s][ID];

         if ( PID < amr->NPatchComma[lv][1] )
         {
            FluxPtr = amr->patch[0][lv][PID]->flux[s];

            if ( FluxPtr != NULL )
            {
               amr->ParaVar->RecvF_IDList[lv][s][ amr->ParaVar->RecvF_NList[lv][s] ] = PID;
               amr->ParaVar->RecvF_NList [lv][s] ++;
            }
         }
      }


//    deallocate memory if no flux data needed to be sent or received
      if ( amr->ParaVar->SendF_NList[lv][s] == 0 )
      {
         delete [] amr->ParaVar->SendF_IDList[lv][s];
         amr->ParaVar->SendF_IDList [lv][s] = NULL;
      }

      if ( amr->ParaVar->RecvF_NList[lv][s] == 0 )
      {
         delete [] amr->ParaVar->RecvF_IDList[lv][s];
         amr->ParaVar->RecvF_IDList [lv][s] = NULL;
      }
   } // for (int s=0; s<6; s++)

} // FUNCTION : Buf_RecordExchangeFluxPatchID




#endif // #ifndef SERIAL
