#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )



// defined in LB_RecordExchangeDataPatchID.cpp
extern void SetTargetLocalID( int NTLocalID[], int *TLocalID[] );
extern void SetTargetSibPID0( const int lv, const int PID0, int SibPID0_List[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_RecordExchangeParticlePatchID
// Description :  Record the buffer patches required to receive particles from their corresponding real patches
//
// Note        :  1. Include buffer patches at both Lv and Lv-1 surrounding real patches at Lv
//                2. Construct RecvPar_NPatch and RecvPar_PIDList at Lv.
//                   Both RecvPar_NPatch and RecvPar_PIDList have the dimension [NLEVEL][2].
//                   [Lv][0/1] is for receiving particles at Lv/Lv-1 buffer patches adjacent to real patches at Lv
//                   --> Mainly for the Poisson Solver at Lv (i.e., for calculating the total density field at Lv)
//                   --> More specific,
//                       [Lv][0] is for receiving the particles of sibling-buffer patches at Lv
//                       adjacent to real patches at Lv
//                       [Lv][1] is for receiving the particles of father-sibling-buffer patches at Lv-1
//                       adjacent to real patches at Lv
//
// Parameter   :  Lv : Target refinement level
//
// Return      :  amr->Par->RecvPar_NPatch[Lv][0/1], amr->Par->RecvPar_PIDList[Lv][0/1]
//-------------------------------------------------------------------------------------------------------
void Par_LB_RecordExchangeParticlePatchID( const int Lv )
{

   const int FaLv     = Lv - 1;
   const int ParLv[2] = { Lv, FaLv };
   const int NParLv   = ( Lv > 0 ) ? 2 : 1;

   int lv, NReal[2], NBuff[2], MemUnit[2], MemSize[2];
   int FaPID, FaSibPID, SibPID, SibPID0, SibPID0_List[26], NTLocalID[26], *TLocalID[26], NPatch_Dup;


// 1. initialize arrays
   for (int t=0; t<NParLv; t++)
   {
      lv         = ParLv[t];
      NReal  [t] = amr->NPatchComma[lv][1];
      NBuff  [t] = amr->NPatchComma[lv][3] - amr->NPatchComma[lv][1];
      MemUnit[t] = MAX( 1, NBuff[t]/4 );  // set arbitrarily (but must > 0)
      MemSize[t] = MemUnit[t];

      if ( amr->Par->RecvPar_PIDList[Lv][t] != NULL )    free( amr->Par->RecvPar_PIDList[Lv][t] );

      amr->Par->RecvPar_NPatch [Lv][t] = 0;
      amr->Par->RecvPar_PIDList[Lv][t] = (int*)malloc( MemSize[t]*sizeof(int) );
   }

// set up the target local indices
   SetTargetLocalID( NTLocalID, TLocalID );


// 2. get the patches at levels Lv and Lv-1 to receive particles
// loop over all "real" patches at Lv with LocalID == 0
   for (int PID0=0; PID0<NReal[0]; PID0+=8)
   {
      SetTargetSibPID0( Lv, PID0, SibPID0_List );

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
               if ( amr->Par->RecvPar_NPatch[Lv][0] >= MemSize[0] )
               {
                  MemSize[0] += MemUnit[0];
                  amr->Par->RecvPar_PIDList[Lv][0] = (int*)realloc( amr->Par->RecvPar_PIDList[Lv][0],
                                                                    MemSize[0]*sizeof(int) );
               }

//             store the target sibling-buffer patch index (note that there may be duplicate PID at this point)
               amr->Par->RecvPar_PIDList[Lv][0][ amr->Par->RecvPar_NPatch[Lv][0] ++ ] = SibPID;
            } // for (int Count=0; Count<NTLocalID[s]; Count++)
         } // if ( SibPID0 >= NReal[0] )

         else if ( SibPID0 == -1  &&  Lv > 0 ) // work for both periodic and non-periodic boundary conditions
         {
            FaPID = amr->patch[0][Lv][PID0]->father;

#           ifdef DEBUG_PARTICLE
            if ( FaPID < 0 )  Aux_Error( ERROR_INFO, "Lv %d, PID0 %d has no father patch (FaPID %d) !!\n", Lv, PID0, FaPID );
#           endif

            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#           ifdef DEBUG_PARTICLE
            if ( FaSibPID < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID0 %d, FaPID %d has no sibling [%d] (FaSibPID = %d) !!\n",
                          Lv, PID0, FaPID, s, FaSibPID );
#           endif

//          check if FaSibPID is a buffer patch
            if ( FaSibPID >= NReal[1] )   // work for both periodic and non-periodic boundary conditions
            {
//             allocate enough memory for the PID array
               if ( amr->Par->RecvPar_NPatch[Lv][1] >= MemSize[1] )
               {
                  MemSize[1] += MemUnit[1];
                  amr->Par->RecvPar_PIDList[Lv][1] = (int*)realloc( amr->Par->RecvPar_PIDList[Lv][1],
                                                                    MemSize[1]*sizeof(int) );
               }

//             store the target father-sibling-buffer patch index (note that there may be duplicate PID at this point)
               amr->Par->RecvPar_PIDList[Lv][1][ amr->Par->RecvPar_NPatch[Lv][1] ++ ] = FaSibPID;
            }
         } // else if ( SibPID0 == -1  &&  Lv > 0 )
      } // for (int s=0; s<26; s++)
   } // for (int PID0=0; PID0<NReal[0]; PID0+=8)


// 3. sort the candidate list and remove duplicates (defined as the patches with the same PID)
//    --> note that patches with different PID may still have the same physical coordinates (i.e., PaddedCr1D)
   for (int t=0; t<NParLv; t++)
   {
      NPatch_Dup = amr->Par->RecvPar_NPatch[Lv][t];

      Mis_Heapsort( NPatch_Dup, amr->Par->RecvPar_PIDList[Lv][t], NULL );

      amr->Par->RecvPar_NPatch[Lv][t] = ( NPatch_Dup > 0 ) ? 1 : 0;

      for (int p=1; p<NPatch_Dup; p++)
         if ( amr->Par->RecvPar_PIDList[Lv][t][p] != amr->Par->RecvPar_PIDList[Lv][t][p-1] )
            amr->Par->RecvPar_PIDList[Lv][t][ amr->Par->RecvPar_NPatch[Lv][t] ++ ] = amr->Par->RecvPar_PIDList[Lv][t][p];

#     ifdef DEBUG_PARTICLE
      for (int p=1; p<amr->Par->RecvPar_NPatch[Lv][t]; p++)
         if ( amr->Par->RecvPar_PIDList[Lv][t][p] <= amr->Par->RecvPar_PIDList[Lv][t][p-1] )
            Aux_Error( ERROR_INFO, "Duplicate PID (Lv %d, t %d, p %d, PID %d, next PID %d) !!\n",
                       Lv, t, p, amr->Par->RecvPar_PIDList[Lv][t][p-1], amr->Par->RecvPar_PIDList[Lv][t][p] );
#     endif
   } // for (int t=0; t<NParLv; t++)


// 4. free memory
   for (int s=0; s<26; s++)   delete [] TLocalID [s];

} // FUNCTION : Par_LB_RecordExchangeParticlePatchID



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
