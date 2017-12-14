#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Refine_AllocateBufferPatch_Sibling
// Description :  Allocate sibling-buffer patches for all real patches at SonLv
//
// Note        :  1. This function is invoked by the function "LB_Refine_AllocateNewPatch"
//                   --> Alternative function of LB_AllocateBufferPatch_Sibling
//                   --> Faster
//                2. This function should be invoked BEFORE LB_AllocateBufferPatch_Father
//                   --> There should be no father-buffer patches at SonLv
//                   --> amr->num[SonLv] should be equal to amr->NPatchComma[SonLv][2]
//                3. This function assumes that all real and buffer patches at SonLv-1 have been allocated,
//                   and the sibling relation at SonLv-1 and the father <-> son relation between SonLv
//                   and SonLv-1 have been constructed
//                4. This function should NOT be applied to the base level
//
// Parameter   :  SonLv : Target refinement level to allocate sibling-buffer patches
//-------------------------------------------------------------------------------------------------------
void LB_Refine_AllocateBufferPatch_Sibling( const int SonLv )
{

// check
   if ( SonLv == 0 )
      Aux_Error( ERROR_INFO, "this function should NOT be applied to the base level !!\n" );

   if ( amr->NPatchComma[SonLv][1] != amr->num[SonLv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d)\" !!\n",
                 SonLv, amr->NPatchComma[SonLv][1], SonLv, amr->num[SonLv] );


// 1. initialize the NotAllocateList (to ensure no duplicate patches at SonLv)
// ==========================================================================================
   const int FaLv     = SonLv - 1;
   const int NFaPatch = amr->num[FaLv];
   const int PScale   = PATCH_SIZE*amr->scale[SonLv];

   bool *NotAllocateList = new bool [NFaPatch];

   for (int FaPID=0; FaPID<NFaPatch; FaPID++)   NotAllocateList[FaPID] = true;



// 2. allocate sibling-buffer patches for all real patches at SonLv
// ==========================================================================================
   int  FaPID, FaSibPID, FaSibSonPID;
   int *FaCr;

   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)
   {
      FaPID = amr->patch[0][SonLv][SonPID0]->father;

#     ifdef GAMER_DEBUG
      if ( FaPID == -1 )   Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d, FaPID == -1 !!\n", SonLv, SonPID0 );
#     endif

      for (int s=0; s<26; s++)
      {
         FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#        ifdef GAMER_DEBUG
         if ( FaSibPID == -1 )   Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d, FaPID %d, FaSibPID == -1 !!\n",
                                            SonLv, SonPID0, FaPID, FaSibPID );
#        endif

//       ensure that it also works with the non-periodic B.C.
         if ( FaSibPID >= 0 )
         {
            FaSibSonPID = amr->patch[0][FaLv][FaSibPID]->son;

//          if son patches exist but are not home --> allocate sibling-buffer patches
            if ( FaSibSonPID <= SON_OFFSET_LB  &&  NotAllocateList[FaSibPID] )
            {
               FaCr = amr->patch[0][FaLv][FaSibPID]->corner;

//             1. data array is not allocated here
//             2. father indices of sibling-buffer patches are always set to -1
//             3. we do not reset the son indices for father patches with sons not home
               amr->pnew( SonLv, FaCr[0],        FaCr[1],        FaCr[2],        -1, false, false );
               amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1],        FaCr[2],        -1, false, false );
               amr->pnew( SonLv, FaCr[0],        FaCr[1]+PScale, FaCr[2],        -1, false, false );
               amr->pnew( SonLv, FaCr[0],        FaCr[1],        FaCr[2]+PScale, -1, false, false );
               amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1]+PScale, FaCr[2],        -1, false, false );
               amr->pnew( SonLv, FaCr[0],        FaCr[1]+PScale, FaCr[2]+PScale, -1, false, false );
               amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1],        FaCr[2]+PScale, -1, false, false );
               amr->pnew( SonLv, FaCr[0]+PScale, FaCr[1]+PScale, FaCr[2]+PScale, -1, false, false );

               amr->NPatchComma[SonLv][2] += 8;

               NotAllocateList[FaSibPID] = false;

            } // if ( FaSibSonPID <= SON_OFFSET_LB  &&  NotAllocateList[FaSibPID] )
         } // if ( FaSibPID >= 0 )
      } // for (int s=0; s<26; s++)
   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)


// set NPatchComma
   for (int m=3; m<28; m++)   amr->NPatchComma[SonLv][m] = amr->NPatchComma[SonLv][2];

   if ( amr->NPatchComma[SonLv][2] != amr->num[SonLv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][2] (%d) != amr->num[%d] (%d)\" !!\n",
                 SonLv, amr->NPatchComma[SonLv][2], SonLv, amr->num[SonLv] );



// 3. record the padded 1D corner coordinates (which can be overwritten by "LB_AllocateBufferPatch_Father")
// ==========================================================================================
   const int NPatch = amr->num[SonLv];

   amr->LB->PaddedCr1DList         [SonLv] = (ulong*)realloc( amr->LB->PaddedCr1DList         [SonLv],
                                                              NPatch*sizeof(ulong) );
   amr->LB->PaddedCr1DList_IdxTable[SonLv] = (int*  )realloc( amr->LB->PaddedCr1DList_IdxTable[SonLv],
                                                              NPatch*sizeof(int ) );

   for (int SonPID=0; SonPID<NPatch; SonPID++)
      amr->LB->PaddedCr1DList[SonLv][SonPID] = amr->patch[0][SonLv][SonPID]->PaddedCr1D;

   Mis_Heapsort( NPatch, amr->LB->PaddedCr1DList[SonLv], amr->LB->PaddedCr1DList_IdxTable[SonLv] );

// check : no duplicate patches at SonLv
#  ifdef GAMER_DEBUG
   for (int t=1; t<NPatch; t++)
   {
      if ( amr->LB->PaddedCr1DList[SonLv][t] == amr->LB->PaddedCr1DList[SonLv][t-1] )
         Aux_Error( ERROR_INFO, "duplicate patches at SonLv %d, PaddedCr1D %lu, PID = %d and %d !!\n",
                    SonLv, amr->LB->PaddedCr1DList[SonLv][t], amr->LB->PaddedCr1DList_IdxTable[SonLv][t],
                    amr->LB->PaddedCr1DList_IdxTable[SonLv][t-1] );
   }
#  endif


// free memory
   delete [] NotAllocateList;

} // FUNCTION : LB_Refine_AllocateBufferPatch_Sibling



#endif // #ifdef LOAD_BALANCE
