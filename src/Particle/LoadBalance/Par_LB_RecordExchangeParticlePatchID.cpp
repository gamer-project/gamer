#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )



// defined in LB_RecordExchangeDataPatchID.cpp
extern void SetTargetLocalID( int NTLocalID[], int *TLocalID[] );
extern void SetTargetSibPID0( const int lv, const int PID0, int SibPID0_List[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_RecordExchangeParticlePatchID
// Description :  Record the patch indices for exchanging particles between different ranks
//
// Note        :  1. R2B : send particles from real to buffer patches
//                   1-1. Include buffer patches at both Lv and Lv-1 surrounding real patches at Lv
//                   1-2. Construct R2B_Real/Buff_NPatchTotal[lv][0/1], R2B_Real/Buff_NPatchEachRank[lv][0/1],
//                        R2B_Real/Buff_PIDList[lv][0/1]
//                        --> They all have the dimension [NLEVEL][2], where [MainLv][0/1] is for receiving particles
//                            at MainLv/MainLv-1 buffer patches adjacent to real patches at MainLv
//                        --> Mainly for the Poisson Solver at MainLv (i.e., for calculating the total density field at MainLv)
//                        --> More specific,
//                            [MainLv][0] is for receiving the particles of sibling-buffer patches at MainLv
//                            adjacent to real patches at MainLv
//                            [MainLv][1] is for receiving the particles of father-sibling-buffer patches at MainLv-1
//                            adjacent to real patches at MainLv
//
//                2. B2R : send particles from buffer to real patches
//                   2-1. Include buffer patches at both Lv and Lv-1 surrounding real patches at Lv
//                   2-2. Construct B2R_Send/Buff_NPatchTotal[lv][0/1], B2R_Send/Buff_NPatchEachRank[lv][0/1],
//                        B2R_Send/Buff_PIDList[lv][0/1]
//                        --> Mainly for sending particles temporarily stored in the buffer patches after updating
//                            their position to their corresponding real patches
//                        --> Similar to the R2B case. But instead of considering **all** real patches, it only
//                            considers **leaf** real patches (since only these patches may have particles)
//                        --> B2R list is a subset of the R2B list
//
// Parameter   :  MainLv   : Main target refinement level
//
// Return      :  R2B_Real/Buff_NPatchTotal[lv][0/1], R2B_Real/Buff_NPatchEachRank[lv][0/1], R2B_Real/Buff_PIDList[lv][0/1]
//                B2R_Send/Buff_NPatchTotal[lv][0/1], B2R_Send/Buff_NPatchEachRank[lv][0/1], B2R_Send/Buff_PIDList[lv][0/1]
//-------------------------------------------------------------------------------------------------------
void Par_LB_RecordExchangeParticlePatchID( const int MainLv )
{

// nothing to do for levels above MAX_LEVEL
   if ( MainLv > MAX_LEVEL )  return;


   const int FaLv         = MainLv - 1;
   const int RelatedLv[2] = { MainLv, FaLv };
   const int NLv          = ( MainLv > 0 ) ? 2 : 1;

   int lv, NReal[2], NBuff[2], MemUnit[2], MemSize[2];
   int FaPID, FaSibPID, SibPID, SibPID0, SibPID0_List[26], NTLocalID[26], *TLocalID[26], Buff_NPatchTotal_Dup;


// 1. initialize arrays
   for (int t=0; t<NLv; t++)
   {
      lv         = RelatedLv[t];
      NReal  [t] = amr->NPatchComma[lv][1];
      NBuff  [t] = amr->NPatchComma[lv][3] - amr->NPatchComma[lv][1];
      MemUnit[t] = MAX( 1, NBuff[t]/4 );  // set arbitrarily (but must > 0)
      MemSize[t] = MemUnit[t];

      if ( amr->Par->R2B_Buff_PIDList[MainLv][t] != NULL )  free( amr->Par->R2B_Buff_PIDList[MainLv][t] );

      amr->Par->R2B_Buff_NPatchTotal[MainLv][t] = 0;
      amr->Par->R2B_Buff_PIDList    [MainLv][t] = (int*)malloc( MemSize[t]*sizeof(int) );
   }

// set up the target local indices
   SetTargetLocalID( NTLocalID, TLocalID );


// 2. get the buffer patches at levels MainLv and MainLv-1 to receive particles
// loop over all "real" patches at MainLv with LocalID == 0
   for (int PID0=0; PID0<NReal[0]; PID0+=8)
   {
      SetTargetSibPID0( MainLv, PID0, SibPID0_List );

      for (int s=0; s<26; s++)
      {
         SibPID0 = SibPID0_List[s];

//       check if SibPID0 exists and is a buffer patch
         if ( SibPID0 >= NReal[0] )    // work for both periodic and non-periodic boundary conditions
         {
            for (int Count=0; Count<NTLocalID[s]; Count++)
            {
               SibPID = SibPID0 + TLocalID[s][Count];

//             allocate enough memory for the PID array
               if ( amr->Par->R2B_Buff_NPatchTotal[MainLv][0] >= MemSize[0] )
               {
                  MemSize[0] += MemUnit[0];
                  amr->Par->R2B_Buff_PIDList[MainLv][0] = (int*)realloc( amr->Par->R2B_Buff_PIDList[MainLv][0],
                                                                         MemSize[0]*sizeof(int) );
               }

//             store the target sibling-buffer patch index (note that there may be duplicate PID at this point)
               amr->Par->R2B_Buff_PIDList[MainLv][0][ amr->Par->R2B_Buff_NPatchTotal[MainLv][0] ++ ] = SibPID;
            } // for (int Count=0; Count<NTLocalID[s]; Count++)
         } // if ( SibPID0 >= NReal[0] )

         else if ( SibPID0 == -1 )  // work for both periodic and non-periodic boundary conditions
         {
#           ifdef DEBUG_PARTICLE
            if ( MainLv == 0 )
               Aux_Error( ERROR_INFO, "Root-level PID0 %d has no sibling patch along sib %d !!\n", PID0, s );
#           endif

            FaPID = amr->patch[0][MainLv][PID0]->father;

#           ifdef DEBUG_PARTICLE
            if ( FaPID < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID0 %d has no father patch (FaPID %d) !!\n", MainLv, PID0, FaPID );
#           endif

            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#           ifdef DEBUG_PARTICLE
            if ( FaSibPID < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID0 %d, FaPID %d has no sibling [%d] (FaSibPID = %d) !!\n",
                          MainLv, PID0, FaPID, s, FaSibPID );
#           endif

//          check if FaSibPID is a buffer patch
            if ( FaSibPID >= NReal[1] )   // work for both periodic and non-periodic boundary conditions
            {
//             allocate enough memory for the PID array
               if ( amr->Par->R2B_Buff_NPatchTotal[MainLv][1] >= MemSize[1] )
               {
                  MemSize[1] += MemUnit[1];
                  amr->Par->R2B_Buff_PIDList[MainLv][1] = (int*)realloc( amr->Par->R2B_Buff_PIDList[MainLv][1],
                                                                         MemSize[1]*sizeof(int) );
               }

//             store the target father-sibling-buffer patch index (note that there may be duplicate PID at this point)
               amr->Par->R2B_Buff_PIDList[MainLv][1][ amr->Par->R2B_Buff_NPatchTotal[MainLv][1] ++ ] = FaSibPID;
            }
         } // if ( SibPID0 >= NReal[0] ) ... else if ( SibPID0 == -1 )
      } // for (int s=0; s<26; s++)
   } // for (int PID0=0; PID0<NReal[0]; PID0+=8)


// 3. sort the candidate list and remove duplicates (defined as the patches with the same PID)
//    --> note that patches with different PID may still have the same physical coordinates (i.e., PaddedCr1D)
   for (int t=0; t<NLv; t++)
   {
      Buff_NPatchTotal_Dup = amr->Par->R2B_Buff_NPatchTotal[MainLv][t];

      Mis_Heapsort( Buff_NPatchTotal_Dup, amr->Par->R2B_Buff_PIDList[MainLv][t], NULL );

      amr->Par->R2B_Buff_NPatchTotal[MainLv][t] = ( Buff_NPatchTotal_Dup > 0 ) ? 1 : 0;

      for (int p=1; p<Buff_NPatchTotal_Dup; p++)
         if ( amr->Par->R2B_Buff_PIDList[MainLv][t][p] != amr->Par->R2B_Buff_PIDList[MainLv][t][p-1] )
            amr->Par->R2B_Buff_PIDList[MainLv][t][ amr->Par->R2B_Buff_NPatchTotal[MainLv][t] ++ ]
               = amr->Par->R2B_Buff_PIDList[MainLv][t][p];

#     ifdef DEBUG_PARTICLE
      for (int p=1; p<amr->Par->R2B_Buff_NPatchTotal[MainLv][t]; p++)
         if ( amr->Par->R2B_Buff_PIDList[MainLv][t][p] <= amr->Par->R2B_Buff_PIDList[MainLv][t][p-1] )
            Aux_Error( ERROR_INFO, "Duplicate PID (Lv %d, t %d, p %d, PID %d, next PID %d) !!\n",
                       MainLv, t, p, amr->Par->R2B_Buff_PIDList[MainLv][t][p-1], amr->Par->R2B_Buff_PIDList[MainLv][t][p] );
#     endif
   } // for (int t=0; t<NLv; t++)


// 4. map buffer patches to real patches
   for (int t=0; t<NLv; t++)
   {
      lv = RelatedLv[t];

      Par_LB_MapBuffer2RealPatch( lv, amr->Par->R2B_Buff_NPatchTotal   [MainLv][t],
                                      amr->Par->R2B_Buff_PIDList       [MainLv][t],
                                      amr->Par->R2B_Buff_NPatchEachRank[MainLv][t],
                                      amr->Par->R2B_Real_NPatchTotal   [MainLv][t],
                                      amr->Par->R2B_Real_PIDList       [MainLv][t],
                                      amr->Par->R2B_Real_NPatchEachRank[MainLv][t] );
   }


// 5. free memory
   for (int s=0; s<26; s++)   delete [] TLocalID [s];

} // FUNCTION : Par_LB_RecordExchangeParticlePatchID



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
