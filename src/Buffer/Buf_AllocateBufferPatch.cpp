#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_AllocateBufferPatch
// Description :  Allocate the buffer patches at level "lv"
//
// Note        :  a. Currently it only works for the functions "Init_Reload" and "Init_BaseLevel"
//                b. For the base level, no data transfer is required
//
// Parameter   :  Tamr : Target AMR_t pointer
//                lv   : Target refinement lv to allocate buffer patches
//-------------------------------------------------------------------------------------------------------
void Buf_AllocateBufferPatch( AMR_t *Tamr, const int lv )
{

   if ( lv < 0  ||  lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


// invoke the function "Buf_AllocateBufferPatch_Base" for the base level
   if ( lv == 0 )
   {
      Buf_AllocateBufferPatch_Base( Tamr );

      return;
   }


   int  BounPID, BuffPID, SonPID, RefinePos, TargetID, NSend[26], NRecv[26];
   int *Send_PosList[26], *Recv_PosList[26];
   bool AllocData[8];   // allocate data or not

   for (int s=0; s<26; s++)
   {
      Send_PosList[s] = NULL;
      Recv_PosList[s] = NULL;
   }

// a. set up the Send_PosList by scanning over the BounP_IDList at level "lv-1"
   for (int s=0; s<26; s++)
   {
//    initialize counter as zero
      NSend[s] = 0;


//    nothing to do if the target direction lies outside the simulation domain for the non-periodic B.C.
      if ( MPI_SibRank[s] < 0 )  continue;


//    allocate the maximum necessary memory
      Send_PosList[s] = new int [ Tamr->ParaVar->BounP_NList[lv-1][s] ];

      for (int ID=0; ID<Tamr->ParaVar->BounP_NList[lv-1][s]; ID++)
      {
         RefinePos = Tamr->ParaVar->BounP_PosList[lv-1][s][ID];
         BounPID   = Tamr->ParaVar->BounP_IDList [lv-1][s][ID];
         SonPID    = Tamr->patch[0][lv-1][BounPID]->son;

         if ( SonPID != -1 )
         {
#           ifdef GAMER_DEBUG
            if ( BounPID >= Tamr->NPatchComma[lv-1][1] )
               Aux_Error( ERROR_INFO, "Boundary patch %d is a buffer patch at level %d !!\n", BounPID, lv-1 );
#           endif

            Send_PosList[s][ NSend[s] ] = RefinePos;
            NSend[s] ++;
         }
      }
   } // for (int s=0; s<26; s++)


// b. get the position of buffer patches to be refined at level "lv"
   MPI_ExchangeBufferPosition( NSend, NRecv, Send_PosList, Recv_PosList );



// c. allocate the buffer patches at level "lv"
   const int scale = Tamr->scale[lv];
   const int Disp  = PATCH_SIZE*scale;
   int *Corner;

   for (int s=0; s<26; s++)
   {
//    nothing to do if the target direction lies outside the simulation domain for the non-periodic B.C.
      if ( MPI_SibRank[s] < 0 )  continue;


//    determine the buffer patches that must store physical data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         if (  TABLE_01( s, 'x', -1, 0, 1 ) == TABLE_02( LocalID, 'x', -1, 1 )  ||
               TABLE_01( s, 'y', -1, 0, 1 ) == TABLE_02( LocalID, 'y', -1, 1 )  ||
               TABLE_01( s, 'z', -1, 0, 1 ) == TABLE_02( LocalID, 'z', -1, 1 )      )
            AllocData[LocalID] = false;

         else
            AllocData[LocalID] = true;
      }


//    find the target buffer patches to be refined
      TargetID = 0;

      for (int ID=0; ID<NRecv[s]; ID++)
      {
         RefinePos = Recv_PosList[s][ID];

         while ( Tamr->ParaVar->BounP_PosList[lv-1][s][TargetID] != RefinePos )
         {
            TargetID++;

//          the BounPID must exist due to the proper-nesting condition
#           ifdef GAMER_DEBUG
            if ( TargetID >= Tamr->ParaVar->BounP_NList[lv-1][s] )
               Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TargetID", TargetID );
#           endif
         }

         BounPID = Tamr->ParaVar->BounP_IDList[lv-1][s][TargetID];
         BuffPID = Tamr->patch[0][lv-1][BounPID]->sibling[s];


//       the BuffPID must exist and be a buffer patch since that the flagging status should be the same over all processes
#        ifdef GAMER_DEBUG
         if ( BuffPID < Tamr->NPatchComma[lv-1][1]  ||  BuffPID >= Tamr->num[lv-1] )
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "BuffPID", BuffPID );
#        endif


//       create buffer patches
         Corner = Tamr->patch[0][lv-1][BuffPID]->corner;

         Tamr->patch[0][lv-1][BuffPID]->son = Tamr->num[lv];

         Tamr->pnew( lv, Corner[0],      Corner[1],      Corner[2],      BuffPID, AllocData[0], AllocData[0] );
         Tamr->pnew( lv, Corner[0]+Disp, Corner[1],      Corner[2],      BuffPID, AllocData[1], AllocData[1] );
         Tamr->pnew( lv, Corner[0],      Corner[1]+Disp, Corner[2],      BuffPID, AllocData[2], AllocData[2] );
         Tamr->pnew( lv, Corner[0],      Corner[1],      Corner[2]+Disp, BuffPID, AllocData[3], AllocData[3] );
         Tamr->pnew( lv, Corner[0]+Disp, Corner[1]+Disp, Corner[2],      BuffPID, AllocData[4], AllocData[4] );
         Tamr->pnew( lv, Corner[0],      Corner[1]+Disp, Corner[2]+Disp, BuffPID, AllocData[5], AllocData[5] );
         Tamr->pnew( lv, Corner[0]+Disp, Corner[1],      Corner[2]+Disp, BuffPID, AllocData[6], AllocData[6] );
         Tamr->pnew( lv, Corner[0]+Disp, Corner[1]+Disp, Corner[2]+Disp, BuffPID, AllocData[7], AllocData[7] );

//       record the number of buffer patches in each sibling direction
         Tamr->NPatchComma[lv][s+2] += 8;

      } // for (int ID=0; ID<NRecv[s]; ID++)

      for (int n=s+3; n<28; n++)    Tamr->NPatchComma[lv][n] = Tamr->num[lv];

   } // for (int s=0; s<26; s++)


   if ( Tamr->NPatchComma[lv][27] != Tamr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][27] (%d) != amr->num[%d] (%d) !!\n",
                 lv, Tamr->NPatchComma[lv][27], lv, Tamr->num[lv] );


   for (int s=0; s<26; s++)
   {
      if ( Send_PosList[s] != NULL )   delete [] Send_PosList[s];
      if ( Recv_PosList[s] != NULL )   delete [] Recv_PosList[s];
   }

} // FUNCTION : Buf_AllocateBufferPatch



#endif // #ifndef SERIAL
