#include "GAMER.h"



#ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SetParID
// Description :  Set the particle unique UID
//
// Note        :  1. This function should be done before initializaion
//                2. Invoked by Init_GAMER() and SF_CreateStar()
//                3. Assuming all the particles is active when initializing
//
// Paramter    :  init : Initialization stage or not
//-------------------------------------------------------------------------------------------------------
void Par_SetParID( const bool init )
{

   long NPar_ThisRank = amr->Par->NPar_AcPlusInac;
   int  NSend[MPI_NRank], SendDisp[MPI_NRank];

   if ( init )
   {
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID (init) ...\n" );

      long NPar_AllRank;
      MPI_Allreduce( &NPar_ThisRank, &NPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

//    assign UID for all ranks
      long_par *ParPUid_AllRank = new long_par [NPar_AllRank];
      for (long p=0; p<NPar_AllRank; p++)   ParPUid_AllRank[p] = (long_par)(p+1);

//    calculate displacement
      int NPar_ThisRank_int = NPar_ThisRank;
      MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );
      if ( MPI_Rank == 0 )
      {
         SendDisp[0] = 0;
         for (int r=1; r<MPI_NRank; r++)   SendDisp[r] = SendDisp[r-1] + NSend[r-1];
      } // if ( MPI_Rank == 0 )

//    scatter particle UID to all ranks
      MPI_Scatterv( ParPUid_AllRank, NSend, SendDisp, MPI_GAMER_LONG_PAR,
                    amr->Par->PUid,  NPar_ThisRank,   MPI_GAMER_LONG_PAR,
                    0, MPI_COMM_WORLD );

//    update the next UID to all ranks
      amr->Par->NextUID = (long_par)(NPar_AllRank+1);

      delete [] ParPUid_AllRank;

      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID (init) ... done\n" );
   }
   else // if ( init )
   {
      long NNewPar_ThisRank = 0L, NNewPar_AllRank;

//    calculate number of new particles
      long *NewParIDList = new long [NPar_ThisRank];
      for (long p=0; p<NPar_ThisRank; p++)
      {
         if ( amr->Par->PUid[p] != (long_par)-1 )   continue;
         NewParIDList[NNewPar_ThisRank] = p;
         NNewPar_ThisRank += 1L;
      }
      int NNewPar_ThisRank_int = NNewPar_ThisRank;
      MPI_Allreduce( &NNewPar_ThisRank, &NNewPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

      if ( NNewPar_AllRank == 0L ) { delete [] NewParIDList; return; }

//    calculate displacement
      MPI_Allgather( &NNewPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, MPI_COMM_WORLD );
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];

      long_par *NewParPUid_AllRank  = new long_par [NNewPar_AllRank];
      long_par *NewParPUid_ThisRank = new long_par [NNewPar_ThisRank];
      real_par *NewParPos_AllRank [3];
      real_par *NewParPos_ThisRank[3];

      for (int d=0; d<3; d++)
      {
         NewParPos_AllRank[d]  = new real_par [NNewPar_AllRank];
         NewParPos_ThisRank[d] = new real_par [NNewPar_ThisRank];
      }

      for (long p=0; p<NNewPar_ThisRank; p++)
      {
         NewParPos_ThisRank[0][p] = amr->Par->PosX[NewParIDList[p]];
         NewParPos_ThisRank[1][p] = amr->Par->PosY[NewParIDList[p]];
         NewParPos_ThisRank[2][p] = amr->Par->PosZ[NewParIDList[p]];
      }

      for (int d=0; d<3; d++)
      {
         MPI_Gatherv( NewParPos_ThisRank[d], NNewPar_ThisRank, MPI_GAMER_REAL_PAR,
                      NewParPos_AllRank [d], NSend, SendDisp,  MPI_GAMER_REAL_PAR,
                      0, MPI_COMM_WORLD );
      }

//    get sorted position and assign ID
      if ( MPI_Rank == 0 )
      {
         long *Sort_IdxTable = new long [NNewPar_AllRank];
         const int Sort_Order[3] = { 0, 1, 2 };
         Mis_SortByRows( NewParPos_AllRank, Sort_IdxTable, NNewPar_AllRank, Sort_Order, 3 );
         for (long p=0; p<NNewPar_AllRank; p++) NewParPUid_AllRank[p] = (long_par)Sort_IdxTable[p] + amr->Par->NextUID;
         delete [] Sort_IdxTable;
      }

//    scatter particle UID to all ranks
      MPI_Scatterv( NewParPUid_AllRank,  NSend, SendDisp,  MPI_GAMER_LONG_PAR,
                    NewParPUid_ThisRank, NNewPar_ThisRank, MPI_GAMER_LONG_PAR,
                    0, MPI_COMM_WORLD );

      for (long p=0; p<NNewPar_ThisRank; p++)
         amr->Par->PUid[NewParIDList[p]] = NewParPUid_ThisRank[p];

      delete [] NewParIDList;
      delete [] NewParPUid_AllRank;
      delete [] NewParPUid_ThisRank;
      for (int d=0; d<3; d++)
      {
         delete [] NewParPos_AllRank[d];
         delete [] NewParPos_ThisRank[d];
      }

//    update the next UID to all ranks
      amr->Par->NextUID += (long_par)NNewPar_AllRank;
   } // if ( init ) ... else ...

   MPI_Barrier( MPI_COMM_WORLD );

} // FUNCTION : Par_SetParID



#endif // #ifdef PARTICLE
