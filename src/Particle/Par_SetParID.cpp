#include "GAMER.h"



#ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    : Par_SetParID
// Description :
//
// Note        : 1. This function should be done before initializaion
//               2. Invoked by Init_GAMER()
//-------------------------------------------------------------------------------------------------------
void Par_SetParID( const bool init )
{

   // if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID ...\n" );

   long NPar_ThisRank, NPar_AllRank;
   NPar_ThisRank = amr->Par->NPar_AcPlusInac;
   int NPar_ThisRank_int = NPar_ThisRank;
   MPI_Allreduce( &NPar_ThisRank, &NPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];

   if ( init )
   {
      MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

      long_par *ParPIdx_AllRank  = new long_par [NPar_AllRank];
      for (long p=0; p<NPar_AllRank; p++)   ParPIdx_AllRank[p] = (long_par)(p+1);

      if ( MPI_Rank == 0 )
      {
         SendDisp[0] = 0;
         for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
      } // if ( MPI_Rank == 0 )

      MPI_Scatterv( ParPIdx_AllRank,  NSend, SendDisp, MPI_GAMER_LONG_PAR,
                    amr->Par->PIdx, NPar_ThisRank,   MPI_GAMER_LONG_PAR,
                    0, MPI_COMM_WORLD );

      delete [] ParPIdx_AllRank;
      amr->Par->NextUniqueIdx = (long_par)(NPar_AllRank+1);
   }
   else // if ( init )
   {

      long NNewPar_ThisRank = 0L, NNewPar_AllRank, pnew_idx;

      for (long p=0; p<NPar_ThisRank; p++) if ( amr->Par->PIdx[p] == (long_par)-1 ) NNewPar_ThisRank += 1L;
      int NNewPar_ThisRank_int = NNewPar_ThisRank;
      MPI_Allreduce( &NNewPar_ThisRank, &NNewPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

      if ( NNewPar_AllRank == 0L )   return;
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID ... new\n" );

      MPI_Allgather( &NNewPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, MPI_COMM_WORLD );

      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];

      long_par *NewParPIdx_AllRank   = new long_par [NNewPar_AllRank];
      long_par *NewParPIdx_ThisRank  = new long_par [NNewPar_ThisRank];
      real_par *NewParPos_AllRank[3];
      real_par *NewParPos_ThisRank[3];

      for (int d=0; d<3; d++)
      {
         NewParPos_AllRank[d]  = new real_par [NNewPar_AllRank];
         NewParPos_ThisRank[d] = new real_par [NNewPar_ThisRank];
      }

      pnew_idx = 0L;
      for (long p=0; p<NPar_ThisRank; p++)
      {
         if ( amr->Par->PIdx[p] != (long_par)-1 )   continue;
         NewParPos_ThisRank[0][pnew_idx] = amr->Par->PosX[p];
         NewParPos_ThisRank[1][pnew_idx] = amr->Par->PosY[p];
         NewParPos_ThisRank[2][pnew_idx] = amr->Par->PosZ[p];
         pnew_idx += 1L;
      }

      for (int d=0; d<3; d++)
      {
         MPI_Gatherv( NewParPos_ThisRank[d], NNewPar_ThisRank, MPI_GAMER_REAL_PAR,
                      NewParPos_AllRank[d], NSend, SendDisp, MPI_GAMER_REAL_PAR,
                      0, MPI_COMM_WORLD );
      }


      // get sorted position and set ID
      if ( MPI_Rank == 0 )
      {
         long *Sort_IdxTable = new long [NNewPar_AllRank];
         const int Sort_Order[3] = { 0, 1, 2 };
         Mis_SortByRows( NewParPos_AllRank, Sort_IdxTable, NNewPar_AllRank, Sort_Order, 3 );
         for (long p=0; p<NNewPar_AllRank; p++) NewParPIdx_AllRank[p] = (long_par)Sort_IdxTable[p] + amr->Par->NextUniqueIdx;
         // for (long p=0; p<NNewPar_AllRank; p++) printf( "Rank: %d, p: %ld, pos: (%+.14e, %+.14e, %+.14e), ID: %ld\n", MPI_Rank, p, NewParPos_AllRank[0][p], NewParPos_AllRank[1][p], NewParPos_AllRank[2][p], (long)NewParPIdx_AllRank[p] );
         delete [] Sort_IdxTable;
      }

      MPI_Barrier( MPI_COMM_WORLD );

      MPI_Scatterv( NewParPIdx_AllRank,  NSend, SendDisp,  MPI_GAMER_LONG_PAR,
                    NewParPIdx_ThisRank, NNewPar_ThisRank, MPI_GAMER_LONG_PAR,
                    0, MPI_COMM_WORLD );
      MPI_Barrier( MPI_COMM_WORLD );

      pnew_idx = 0L;
      for (long p=0; p<NPar_ThisRank; p++)
      {
         if ( amr->Par->PIdx[p] != (long_par)-1 )   continue;
         amr->Par->PIdx[p] = NewParPIdx_ThisRank[pnew_idx];
         // amr->Par->PIdx[p] = (long_par)10;
         // printf( "Rank: %d, p: %ld, ID: %ld\n", MPI_Rank, pnew_idx, (long)NewParPIdx_ThisRank[pnew_idx] );
         pnew_idx += 1L;
      }

      delete [] NewParPIdx_AllRank;
      delete [] NewParPIdx_ThisRank;
      for (int d=0; d<3; d++)
      {
         delete [] NewParPos_AllRank[d];
         delete [] NewParPos_ThisRank[d];
      }

      amr->Par->NextUniqueIdx += (long_par)NNewPar_AllRank;
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID ... done\n" );
   } // if ( init ) ... else ...


   MPI_Barrier( MPI_COMM_WORLD );

} // FUNCTION : Par_SetParID
#endif // #ifdef PARTICLE
