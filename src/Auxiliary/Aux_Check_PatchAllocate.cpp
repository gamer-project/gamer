#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_PatchAllocate
// Description :  Verify several things:
//                1. The corner coordinates of all patches (real+buffer) in all ranks lie inside
//                   the simulation box (padded with two buffer patches on each side in each direction for the periodic B.C.)
//                2. The corner coordinates of all level lv patches (real+buffer) in all ranks are multiples
//                   of the patch scale at level lv
//                3. All patches (real+buffer) in the same rank have different corner coordinates
//                4. All real patches in all ranks have different corner coordinates
//                5. The corner coordinates of any buffer patch do match the corner coordinates of one real patch
//                   (unless they lie outside the simulation box for the periodic B.C.)
//                6. Check if sibling patches have the correct coordinates
//                7. Work for both periodic and non-periodic BC's
//
// Note           This check will not detect the useless buffer patches (which do map to real patches but
//                do not provide any useful information during simulations)
//
// Parameter   :  lv      : Target refinement level
//                comment : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_PatchAllocate( const int lv, const char *comment )
{

// check
   if ( lv < 0  || lv >= NLEVEL )
      Aux_Error( ERROR_INFO, "lv %d lies outside the permitted range [0 ... %d] !!\n", lv, NLEVEL-1 );


   const int PScale             = PATCH_SIZE*amr->scale[lv];
   const int PGScale            = 2*PScale;
   const int BoxScale_Padded[3] = { amr->BoxScale[0] + 2*PGScale,
                                    amr->BoxScale[1] + 2*PGScale,
                                    amr->BoxScale[2] + 2*PGScale };
   int    NReal = amr->NPatchComma[lv][ 1];
   int    NTot  = amr->NPatchComma[lv][27];
   int    NBuff = NTot - NReal;
   int    Pass  = true;
   int    Cr_Padded[3], NReal_Tot=-1, NBuff_Tot=-1; 
   int   *NReal_Each=NULL, *NReal_Disp=NULL, *NBuff_Each=NULL, *NBuff_Disp=NULL; 
   ulong *Cr1D_Real=NULL, *Cr1D_Buff=NULL;
   char  *Match=NULL;
   ulong *Cr1D      = new ulong [ NReal+NBuff ];
   ulong *Cr1D_Sort = new ulong [ NReal+NBuff ];


// avoid over-flow
   if (  log2( (double)BoxScale_Padded[0]*(double)BoxScale_Padded[1]*(double)BoxScale_Padded[2] )
         > 8.0*sizeof(long)-1.0  )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : the check \"%s\" is disabled due to long integer overflow !!\n",
                      __FUNCTION__ );

      OPT__CK_PATCH_ALLOCATE = false;
      return;
   }


   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<NTot; PID++)
         {
            for (int d=0; d<3; d++)    
            {
               Cr_Padded[d] = amr->patch[0][lv][PID]->corner[d] + PGScale;

//             check 1 
               if ( PID < NReal  ||  OPT__BC_FLU[2*d] != BC_FLU_PERIODIC )    // real patches or non-periodic B.C
               {
                  if ( Cr_Padded[d] < PGScale  ||  Cr_Padded[d] > BoxScale_Padded[d]-PGScale-PScale )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Pass = false;
                     }

                     Aux_Message( stderr, "Check 1, Rank %4d, PID %8d (real patch/non-periodic B.C.): Cr_Padded[%d] = %7d\n",
                                  MPI_Rank, PID, d, Cr_Padded[d] );
                  }
               }
               else // buffer patches in periodic B.C.
               {
                  if ( Cr_Padded[d] < 0  ||  Cr_Padded[d] > BoxScale_Padded[d]-PScale )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Pass = false;
                     }

                     Aux_Message( stderr, "Check 1, Rank %4d, PID %8d (buff patch): Cr_Padded[%d] = %7d\n",
                                  MPI_Rank, PID, d, Cr_Padded[d] );
                  }
               }


//             check 2
               if ( Cr_Padded[d] % PScale != 0 )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Pass = false;
                  }

                  Aux_Message( stderr, "Check 2, Rank %4d, PID %8d: Cr_Padded[%d] = %7d\n",
                               MPI_Rank, PID, d, Cr_Padded[d] );
               }
            } // for (int d=0; d<3; d++)

//          check 6
            for (int s=0; s<26; s++)
            {
               const int SibPID = amr->patch[0][lv][PID]->sibling[s];

               if ( SibPID >= 0 )
               {
                  for (int d=0; d<3; d++)
                  {
                     if (    amr->patch[0][lv][   PID]->corner[d] + TABLE_01( s, 'x'+d, -PScale, 0, PScale )
                          != amr->patch[0][lv][SibPID]->corner[d]  )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, lv, Time[lv], Step );
                           Pass = false;
                        }

                        Aux_Message( stderr, "Check 6, Rank %4d, sibling %2d, PID %8d, SibPID %8d, dim %d: corner %8d + %8d != %8d\n",
                                     MPI_Rank, s, PID, SibPID, d, amr->patch[0][lv][PID]->corner[d],
                                     TABLE_01( s, 'x'+d, -PScale, 0, PScale ), amr->patch[0][lv][SibPID]->corner[d] );
                     }
                  }
               } // if ( SibPID >= 0 )
            } // for (int s=0; s<26; s++)


            Cr1D[PID] = Mis_Idx3D2Idx1D( BoxScale_Padded, Cr_Padded );

         } // for (int PID=0; PID<NTot; PID++)


//       check 3         
         memcpy( Cr1D_Sort, Cr1D, NTot*sizeof(ulong) );
         Mis_Heapsort( NTot, Cr1D_Sort, NULL );
         for (int t=0; t<NTot-1; t++)
         {
            if ( Cr1D_Sort[t] == Cr1D_Sort[t+1] )
            {
               if ( Pass )
               {
                  Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                               comment, __FUNCTION__, lv, Time[lv], Step );
                  Pass = false;
               }

               Aux_Message( stderr, "Check 3, Rank %4d: Cr1D_Sort[%d] %8lu\n", MPI_Rank, t, Cr1D_Sort[t] );
            }
         }
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


// gather data from all MPI ranks for the checks 4 and 5
   if ( MPI_Rank == 0 )
   {
      NReal_Each = new int [MPI_NRank];
      NReal_Disp = new int [MPI_NRank];
      NBuff_Each = new int [MPI_NRank];
      NBuff_Disp = new int [MPI_NRank];

      NReal_Disp[0] = 0;
      NBuff_Disp[0] = 0;
   }

#  ifndef SERIAL
   MPI_Gather( &NReal, 1, MPI_INT, NReal_Each, 1, MPI_INT, 0, MPI_COMM_WORLD );
   MPI_Gather( &NBuff, 1, MPI_INT, NBuff_Each, 1, MPI_INT, 0, MPI_COMM_WORLD );
#  else
   NReal_Each[0] = NReal;
   NBuff_Each[0] = NBuff;
#  endif

   if ( MPI_Rank == 0 )
   {
      for (int r=1; r<MPI_NRank; r++)
      {
         NReal_Disp[r] = NReal_Disp[r-1] + NReal_Each[r-1];
         NBuff_Disp[r] = NBuff_Disp[r-1] + NBuff_Each[r-1];
      }

      NReal_Tot = NReal_Disp[MPI_NRank-1] + NReal_Each[MPI_NRank-1];
      NBuff_Tot = NBuff_Disp[MPI_NRank-1] + NBuff_Each[MPI_NRank-1];

      Cr1D_Real = new ulong [NReal_Tot];
      Cr1D_Buff = new ulong [NBuff_Tot];
      Match     = new  char [NBuff_Tot];
   }
   
#  ifndef SERIAL
   MPI_Gatherv( Cr1D,       NReal, MPI_UNSIGNED_LONG, Cr1D_Real, NReal_Each, NReal_Disp, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
   MPI_Gatherv( Cr1D+NReal, NBuff, MPI_UNSIGNED_LONG, Cr1D_Buff, NBuff_Each, NBuff_Disp, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD );
#  else
   for (int t=0; t<NReal; t++)   Cr1D_Real[t] = Cr1D[t];
   for (int t=0; t<NBuff; t++)   Cr1D_Buff[t] = Cr1D[t+NReal];
#  endif


   if ( MPI_Rank == 0 )
   {
//    apply periodic boundary condition for the external buffer patches     
      if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_FLU[2] == BC_FLU_PERIODIC  ||  OPT__BC_FLU[4] == BC_FLU_PERIODIC )
      for (int t=0; t<NBuff_Tot; t++)
      {
         Mis_Idx1D2Idx3D( BoxScale_Padded, Cr1D_Buff[t], Cr_Padded );

         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[2*d] == BC_FLU_PERIODIC )
            {
               if      ( Cr_Padded[d] < PGScale )                             Cr_Padded[d] += amr->BoxScale[d];
               else if ( Cr_Padded[d] > BoxScale_Padded[d]-PGScale-PScale )   Cr_Padded[d] -= amr->BoxScale[d];
            }
         }

         Cr1D_Buff[t] = Mis_Idx3D2Idx1D( BoxScale_Padded, Cr_Padded );
      }


//    sort all real and buffer patches   
      Mis_Heapsort( NReal_Tot, Cr1D_Real, NULL );
      Mis_Heapsort( NBuff_Tot, Cr1D_Buff, NULL );


//    check 4      
      for (int t=0; t<NReal_Tot-1; t++)
      {
         if ( Cr1D_Real[t] == Cr1D_Real[t+1] )
         {
            if ( Pass )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                            comment, __FUNCTION__, lv, Time[lv], Step );
               Pass = false;
            }

            Aux_Message( stderr, "Check 4, Rank %4d: Cr1D_Real[%d] %8lu\n", MPI_Rank, t, Cr1D_Real[t] );
         }
      }


//    check 5      
      Mis_Matching_char( NReal_Tot, Cr1D_Real, NBuff_Tot, Cr1D_Buff, Match );

      for (int t=0; t<NBuff_Tot; t++)
      {
         if ( Match[t] != 1 )
         {
            if ( Pass )
            {
               Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                            comment, __FUNCTION__, lv, Time[lv], Step );
               Pass = false;
            }

            Aux_Message( stderr, "Check 5, Rank %4d: Cr1D_Buff[%d] %8lu\n", MPI_Rank, t, Cr1D_Buff[t] );
         }
      }
   } // if ( MPI_Rank == 0 )


   if ( Pass )
   {
      if ( MPI_Rank == 0 )    
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n", 
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }
   else
      MPI_Exit();


   delete [] Cr1D;
   delete [] Cr1D_Sort;
   if ( MPI_Rank == 0 )
   {
      delete [] NReal_Each;
      delete [] NReal_Disp;
      delete [] NBuff_Each;
      delete [] NBuff_Disp;
      delete [] Cr1D_Real;
      delete [] Cr1D_Buff;
      delete [] Match;
   }

} // FUNCTION : Aux_Check_PatchAllocate


