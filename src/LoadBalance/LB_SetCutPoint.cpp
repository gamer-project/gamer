#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SetCutPoint
// Description :  Set the range of LB_Idx for distributing patches to different ranks
//
// Note        :  1. Set the input array CutPoint[]
//                2. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be sent to rank "r"
//                3. Option "InputLBIdx0AndLoad" is useful during RESTART where we have very limited information
//                   (e.g., we don't know the number of patches in each rank, amr->NPatchComma, and any
//                   particle information yet ...)
//                   --> See the description of "InputLBIdx0AndLoad, LBIdx0_AllRank_Input, and
//                       Load_AllRank_Input" below
//
// Parameter   :  lv                   : Target refinement level
//                NPG_Total            : Total number of patch groups on level "lv"
//                CutPoint             : Cut point array to be set
//                InputLBIdx0AndLoad   : Provide both LBIdx0_AllRank_Input[] and Load_AllRank_Input[] directly
//                                       so that they don't have to be collected from all ranks again
//                                       --> Useful during RESTART
//                LBIdx0_AllRank_Input : LBIdx of all patch groups in all ranks
//                                       --> Useful only when InputLBIdx0AndLoad == true
//                                       --> Only need the **minimum** LBIdx in each patch group
//                                       --> Only rank 0 needs to provide this list
//                                       --> Can be unsorted
//                Load_AllRank_Input   : Load-balance weighting of all patch groups in all ranks
//                                       --> Useful only when InputLBIdx0AndLoad == true
//                                       --> Please provide the **sum** of all patches within each patch group
//                                       --> Only rank 0 needs to provide this list
//                                       --> Must be in the same order as LBIdx0_AllRank_Input
//                ParWeight            : Relative load-balance weighting of particles
//                                       --> Weighting of each patch is estimated as "PATCH_SIZE^3 + NParThisPatch*ParWeight"
//                                       --> <= 0.0 : do not consider particle weighting
//
// Return      :  CutPoint[]
//-------------------------------------------------------------------------------------------------------
void LB_SetCutPoint( const int lv, const int NPG_Total, long *CutPoint, const bool InputLBIdx0AndLoad,
                     long *LBIdx0_AllRank_Input, double *Load_AllRank_Input, const double ParWeight )
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "      %s at Lv %2d ...\n", __FUNCTION__, lv );


// check
   if ( MPI_Rank == 0  &&  InputLBIdx0AndLoad  &&  ( LBIdx0_AllRank_Input == NULL || Load_AllRank_Input == NULL )  )
      Aux_Error( ERROR_INFO, "LBIdx0_AllRank_Input/Load_AllRank_Input == NULL when InputLBIdx0AndLoad is on !!\n" );

   if ( NPG_Total < 0 )
      Aux_Error( ERROR_INFO, "NPG_Total (%d) < 0 !!\n", NPG_Total );


// 1. collect the load-balance weighting and LB_Idx of all patch groups from all ranks
   long   *LBIdx0_AllRank = NULL;
   double *Load_AllRank   = NULL;
   int    *IdxTable       = ( MPI_Rank == 0 ) ? new int [NPG_Total] : NULL;

// use the input tables directly
// --> useful during RESTART, where we have very limited information
//     (e.g., we don't know the number of patches in each rank, amr->NPatchComma, and any particle information yet ...)
   if ( InputLBIdx0AndLoad )
   {
      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = LBIdx0_AllRank_Input;
         Load_AllRank   = Load_AllRank_Input;
      }
   }

   else
   {
//    allocate memory
      int NPG_ThisRank = amr->NPatchComma[lv][1] / 8;

      long   *LBIdx0_ThisRank = new long   [ NPG_ThisRank ];
      double *Load_ThisRank   = new double [ NPG_ThisRank ];
      int    *NPG_EachRank    = NULL;
      int    *Recv_Disp       = NULL;

      if ( MPI_Rank == 0 )
      {
         NPG_EachRank   = new int    [ MPI_NRank ];
         Recv_Disp      = new int    [ MPI_NRank ];
         LBIdx0_AllRank = new long   [ NPG_Total ];
         Load_AllRank   = new double [ NPG_Total ];
      }


//    collect the number of patch groups in each rank
      MPI_Gather( &NPG_ThisRank, 1, MPI_INT, NPG_EachRank, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         Recv_Disp[0] = 0;
         for (int r=0; r<MPI_NRank-1; r++)   Recv_Disp[r+1] = Recv_Disp[r] + NPG_EachRank[r];
      }


//    collect the minimum LBIdx in each patch group
//    --> assuming patches within the same patch group have consecutive LBIdx
      for (int t=0; t<NPG_ThisRank; t++)
      {
         const int PID0 = t*8;

         LBIdx0_ThisRank[t]  = amr->patch[0][lv][PID0]->LB_Idx;
         LBIdx0_ThisRank[t] -= LBIdx0_ThisRank[t] % 8;         // get the **minimum** LBIdx in this patch group
      }

      MPI_Gatherv( LBIdx0_ThisRank, NPG_ThisRank, MPI_LONG, LBIdx0_AllRank, NPG_EachRank, Recv_Disp,
                   MPI_LONG, 0, MPI_COMM_WORLD );


//    collect the load-balance weighting in each patch group
      LB_EstimateWorkload_AllPatchGroup( lv, ParWeight, Load_ThisRank );

      MPI_Gatherv( Load_ThisRank, NPG_ThisRank, MPI_DOUBLE, Load_AllRank, NPG_EachRank, Recv_Disp,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD );


//    free memory
      delete [] LBIdx0_ThisRank;
      delete [] Load_ThisRank;

      if ( MPI_Rank == 0 )
      {
         delete [] NPG_EachRank;
         delete [] Recv_Disp;
      }
   } // if ( InputLBIdx0AndLoad ) ... else ...


   if ( MPI_Rank == 0 )
   {
      double *Load_Record = ( OPT__VERBOSE ) ? new double [MPI_NRank] : NULL;
      double  Load_Ave;

//    3. sort LB_Idx
//    --> after sorting, we must use IdxTable to access the Load_AllRank[] array
      Mis_Heapsort( NPG_Total, LBIdx0_AllRank, IdxTable );


//    4. set the cut points
      for (int t=0; t<MPI_NRank+1; t++)   CutPoint[t] = -1;

//    4-1. take care of the case with no patches at all
      if ( NPG_Total == 0 )
      {
         Load_Ave = 0.0;

         if ( OPT__VERBOSE )
            for (int r=0; r<MPI_NRank; r++)  Load_Record[r] = 0.0;
      }

      else
      {
//       4-2. get the average workload for each rank
         Load_Ave = 0.0;
         for (int t=0; t<NPG_Total; t++)  Load_Ave += Load_AllRank[t];
         Load_Ave /= (double)MPI_NRank;

//       4-3. set the min and max cut points
         const long LBIdx0_Min = LBIdx0_AllRank[             0 ];
         const long LBIdx0_Max = LBIdx0_AllRank[ NPG_Total - 1 ];
         CutPoint[        0] = LBIdx0_Min;
         CutPoint[MPI_NRank] = LBIdx0_Max + 8;  // +8 since the maximum LBIdx in all patches is LBIdx0_Max + 7

//       4-4. find the LBIdx with an accumulated workload (LoadAcc) closest to the average workload of each rank (LoadTarget)
         int    CutIdx     = 1;                 // target array index for CutPoint[]
                                                // --> note that CutPoint[CutIdx] is the **exclusive** upper bound of rank "CutIdx-1"
         double LoadAcc    = 0.0;               // accumulated workload
         double LoadTarget = CutIdx*Load_Ave;   // target accumulated workload for the rank "CutIdx-1"
         double LoadThisPG;                     // workload of the target patch group

         for (int PG=0; PG<NPG_Total; PG++)
         {
//          nothing to do if all cut points have been set already
            if ( CutIdx == MPI_NRank )    break;

//          remember to use IdxTable to access Load_AllRank
            LoadThisPG = Load_AllRank[ IdxTable[PG] ];

//          check if adding a new patch group will exceed the target accumulated workload
            if ( LoadAcc+LoadThisPG >= LoadTarget )
            {
//             determine the cut point with an accumulated workload **closest** to the target accumulated workload
//             (a) if adding a new patch group will exceed the target accumulated workload too much
//                 --> exclude this patch group from the rank "CutIdx-1"
//             note that both "LoadAcc > LoadTarget" and "LoadAcc <= LoadTaget" can happen
               if ( fabs(LoadAcc-LoadTarget) < LoadAcc+LoadThisPG-LoadTarget )
               {
                  CutPoint[CutIdx] = LBIdx0_AllRank[PG];

                  PG --;   // because this patch group has been **excluded** from this cut point
               }

//             (b) if adding a new patch group will NOT exceed the target accumulated workload too much
//                 --> include this patch group in the rank "CutIdx-1"
               else
               {
//                be careful about the special case "PG == NPG_Total-1"
                  CutPoint[CutIdx] = ( PG == NPG_Total-1 ) ? CutPoint[MPI_NRank] : LBIdx0_AllRank[PG+1];

                  LoadAcc += LoadThisPG;
               }

//             record the **accumulated** workload of each rank
               if ( OPT__VERBOSE )  Load_Record[ CutIdx - 1 ] = LoadAcc;

               CutIdx ++;
               LoadTarget = CutIdx*Load_Ave;
            } // if ( LoadAcc+LoadThisPG >= LoadTarget )

            else
            {
               LoadAcc += LoadThisPG;
            } // if ( LoadAcc+LoadThisPG >= LoadTarget ) ... else ...

         } // for (int PG=0; PG<NPG_Total; PG++)

//       4.5 take care of the special case where the last several ranks have no patches at all
         for (int t=CutIdx; t<MPI_NRank; t++)
         {
            CutPoint[t] = CutPoint[MPI_NRank];

            if ( OPT__VERBOSE )  Load_Record[ t - 1 ] = Load_Ave*MPI_NRank;
         }

//       4.6 check
#        ifdef GAMER_DEBUG
//       all cut points must be set properly
         for (int t=0; t<MPI_NRank+1; t++)
            if ( CutPoint[t] == -1 )
               Aux_Error( ERROR_INFO, "lv %d, CutPoint[%d] == -1 !!\n", lv, t );

//       monotonicity
         for (int t=0; t<MPI_NRank; t++)
            if ( CutPoint[t+1] < CutPoint[t] )
               Aux_Error( ERROR_INFO, "lv %d, CutPoint[%d] (%ld) < CutPoint[%d] (%ld) !!\n",
                          lv, t+1, CutPoint[t+1], t, CutPoint[t] );
#        endif
      } // if ( NPG_Total == 0 ) ... else ...


//    5. output the cut points and workload of each MPI rank
      if ( OPT__VERBOSE )
      {
//       convert the accumulated workload to the actual workload of each rank
         Load_Record[ MPI_NRank - 1 ] = Load_Ave*MPI_NRank;
         for (int r=MPI_NRank-1; r>=1; r--)  Load_Record[r] -= Load_Record[r-1];

         double Load_Max = -1.0;

         for (int r=0; r<MPI_NRank; r++)
         {
            Aux_Message( stdout, "         Lv %2d: Rank %4d, Cut %15ld -> %15ld, Load_Weighted %9.3e\n",
                         lv, r, CutPoint[r], CutPoint[r+1], Load_Record[r] );

            if ( Load_Record[r] > Load_Max )    Load_Max = Load_Record[r];
         }

         Aux_Message( stdout, "         Load_Ave %9.3e, Load_Max %9.3e --> Load_Imbalance = %6.2f%%\n",
                      Load_Ave, Load_Max, (NPG_Total == 0) ? 0.0 : 100.0*(Load_Max-Load_Ave)/Load_Ave );
         Aux_Message( stdout, "         =============================================================================\n" );

         delete [] Load_Record;
      }
   } // if ( MPI_Rank == 0 )


// 6. broadcast the cut points
   MPI_Bcast( CutPoint, MPI_NRank+1, MPI_LONG, 0, MPI_COMM_WORLD );


// free memory
   if ( MPI_Rank == 0 )
   {
      delete [] IdxTable;

      if ( !InputLBIdx0AndLoad )
      {
         delete [] LBIdx0_AllRank;
         delete [] Load_AllRank;
      }
   }


   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "      %s at Lv %2d ... done\n", __FUNCTION__, lv );

} // FUNCTION : LB_SetCutPoint



#endif // #ifdef LOAD_BALANCE
