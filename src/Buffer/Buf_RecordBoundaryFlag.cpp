#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_RecordBoundaryFlag
// Description :  Record the flags of boundary patches in amr->ParaVar->BounFlag_PosList[] for
//                MPI_ExchangeBoundaryFlag()
//
// Note        :  1. Invoked by Flag_Real()
//                2. No OpenMP directives are applied in this function since the counter
//                   amr->ParaVar->BounFlag_NList[]
//
// Parameter   :  lv : Target refinement level to be flagged
//-------------------------------------------------------------------------------------------------------
void Buf_RecordBoundaryFlag( const int lv )
{

// check
   if ( lv < 0  ||  lv >= TOP_LEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


   int FlagPos, PID, SibPID, FlagLayer, Sib;

// begin the main loop of Buf_RecordBoundaryFlag()
   for (int s=0; s<26; s++)
   {
//    initialize the counter as zero
      amr->ParaVar->BounFlag_NList[lv][s] = 0;


//    deallocate memory
      if ( amr->ParaVar->BounFlag_PosList[lv][s] != NULL )
      {
         delete [] amr->ParaVar->BounFlag_PosList[lv][s];
         amr->ParaVar->BounFlag_PosList[lv][s] = NULL;
      }


//    nothing to do if there are no boundary patches at the target direction
      if ( amr->ParaVar->BounP_NList[lv][s] == 0 )     continue;


//    nothing to do if the target direction lies outside the simulation domain for the non-periodic B.C.
      if ( MPI_SibRank[s] < 0 )  continue;


//    set up the FlagLayer
      FlagLayer = TABLE_05( s );


//    allocate the maximum necessary memory (which will be deallocated by MPI_ExchangeBoundaryFlag())
      amr->ParaVar->BounFlag_PosList[lv][s] = new int [ FlagLayer*amr->ParaVar->BounP_NList[lv][s] ];


//    fill up BounFlag_PosList[lv][s][]
      for (int ID=0; ID<amr->ParaVar->BounP_NList[lv][s]; ID++)
      {
         FlagPos = amr->ParaVar->BounP_PosList[lv][s][ID];
         PID     = amr->ParaVar->BounP_IDList [lv][s][ID];

#        ifdef GAMER_DEBUG
         if ( PID < 0  ||  PID >= amr->num[lv] )   Aux_Error( ERROR_INFO, "incorrect PID = %d !!", PID );
#        endif

         if ( amr->patch[0][lv][PID]->flag )
         {
//          record the flag of layer 0
            amr->ParaVar->BounFlag_PosList[lv][s][ amr->ParaVar->BounFlag_NList[lv][s] ] = FlagPos;

//          record the flags of layers > 0 (only when the layer 0 is flagged)
            for (int n=1; n<FlagLayer; n++)
            {
               Sib    = TABLE_06( s, n );
               SibPID = amr->patch[0][lv][PID]->sibling[Sib];

#              ifdef GAMER_DEBUG
               if ( SibPID <= SIB_OFFSET_NONPERIODIC )   Aux_Error( ERROR_INFO, "incorrect SibPID = %d !!\n", SibPID );
#              endif

//             set BounFlag_PosList = BUFFER_IS_FLAGGED to indicate that this layer is flagged (only for layer > 0)
               if ( SibPID >= 0  &&  amr->patch[0][lv][SibPID]->flag )
                  amr->ParaVar->BounFlag_PosList[lv][s][ amr->ParaVar->BounFlag_NList[lv][s]+n ] = BUFFER_IS_FLAGGED;
               else
                  amr->ParaVar->BounFlag_PosList[lv][s][ amr->ParaVar->BounFlag_NList[lv][s]+n ] = -1;
            }

            amr->ParaVar->BounFlag_NList[lv][s] += FlagLayer;

         } // if ( amr->patch[0][lv][PID]->flag )
      } // for (int ID=0; ID<amr->ParaVar->BounP_NList[lv][s]; ID++)
   } // for (int s=0; s<26; s++)

} // FUNCTION : Buf_RecordBoundaryFlag



#endif // #ifndef SERIAL
