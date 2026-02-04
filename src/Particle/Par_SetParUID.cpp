#include "GAMER.h"



#ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_SetParUID
// Description :  Set the particle unique UID
//
// Note        :  1. Invoked by Init_GAMER() and SF_CreateStar()
//                2. The new particle UID should be initialized as PPUID_TBA
//                3. Currently, the new particles UID are given by the sorted position
//
// Paramter    :  init : Initialization stage or not
//-------------------------------------------------------------------------------------------------------
void Par_SetParUID( const bool init )
{

   long NPar_ThisRank = amr->Par->NPar_AcPlusInac;
   int  NSend[MPI_NRank], SendDisp[MPI_NRank];

   if ( init  &&  MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID (init) ...\n" );

   if ( amr->Par->NextUID < 1L )
      Aux_Error( ERROR_INFO, "amr->Par->NextUID(%ld) is invalid !!\n", amr->Par->NextUID );

   long NNewPar_ThisRank = 0L, NNewPar_AllRank;

// calculate number of new particles
   long *NewParIDList = new long [NPar_ThisRank];
   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( amr->Par->PUid[p] > (long_par)0  &&  amr->Par->PUid[p] < amr->Par->NextUID )   continue;  // exclude particles that already have valid UID

      if ( amr->Par->PUid[p] != PPUID_TBA )
         Aux_Error( ERROR_INFO, "New particle before ParUID assignment has an invalid PUid[%ld] = %ld != %ld !!\n", p, (long)amr->Par->PUid[p], (long)PPUID_TBA );

      if ( amr->Par->Mass[p] < (real_par)0.0 )   continue;   // exclude inactive particles

      NewParIDList[NNewPar_ThisRank] = p;
      NNewPar_ThisRank += 1L;
   }
   MPI_Allreduce( &NNewPar_ThisRank, &NNewPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

// check the number of new particles
   if ( NNewPar_AllRank == 0L )
   {
      delete [] NewParIDList;
      if ( init  &&  MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID (init) ... done\n" );
      return;
   }
   if ( NNewPar_AllRank <  0L )
      Aux_Error( ERROR_INFO, "NNewPar_AllRank(%ld) is negative !!\n", NNewPar_AllRank );

// check whether the next UID will overflow
   if ( amr->Par->NextUID >= ( (sizeof(long_par) == sizeof(int)) ? __INT_MAX__ : __LONG_MAX__ ) - NNewPar_AllRank )
      Aux_Error( ERROR_INFO, "amr->Par->PUid = amr->Par->NextUID(%ld) + NNewPar_AllRank(%ld) will overflow !!\n", amr->Par->NextUID, NNewPar_AllRank );

// safety check for the type conversion from long to int
   int NNewPar_ThisRank_int = NNewPar_ThisRank;
   if ( NNewPar_ThisRank_int != NNewPar_ThisRank )
      Aux_Error( ERROR_INFO, "NNewPar_ThisRank_int(%d) != NNewPar_ThisRank(%ld) !!\n", NNewPar_ThisRank_int, NNewPar_ThisRank );

// calculate displacement
   MPI_Allgather( &NNewPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, MPI_COMM_WORLD );
   SendDisp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];

// allocate memory
   long_par *NewParPUid_AllRank  = new long_par [NNewPar_AllRank];
   long_par *NewParPUid_ThisRank = new long_par [NNewPar_ThisRank];
   real_par *NewParPos_AllRank [3];
   real_par *NewParPos_ThisRank[3];

   for (int d=0; d<3; d++)
   {
      NewParPos_AllRank[d]  = new real_par [NNewPar_AllRank];
      NewParPos_ThisRank[d] = new real_par [NNewPar_ThisRank];
   }

// collect new particle position in this rank
   for (long p=0; p<NNewPar_ThisRank; p++)
   {
      NewParPos_ThisRank[0][p] = amr->Par->PosX[NewParIDList[p]];
      NewParPos_ThisRank[1][p] = amr->Par->PosY[NewParIDList[p]];
      NewParPos_ThisRank[2][p] = amr->Par->PosZ[NewParIDList[p]];
   }

// gather particle position from all ranks
   for (int d=0; d<3; d++)
   {
      MPI_Gatherv( NewParPos_ThisRank[d], NNewPar_ThisRank_int, MPI_GAMER_REAL_PAR,
                   NewParPos_AllRank [d], NSend, SendDisp,      MPI_GAMER_REAL_PAR,
                   0, MPI_COMM_WORLD );
   }

// get sorted position and assign UID
   if ( MPI_Rank == 0 )
   {
      long *Sort_IdxTable = new long [NNewPar_AllRank];
      const int Sort_Order[3] = { 0, 1, 2 };
      Mis_SortByRows( NewParPos_AllRank, Sort_IdxTable, NNewPar_AllRank, Sort_Order, 3 );
      for (long p=0; p<NNewPar_AllRank; p++) NewParPUid_AllRank[ Sort_IdxTable[p] ] = (long_par)( p + amr->Par->NextUID );
      delete [] Sort_IdxTable;
   }

// scatter particle UID to all ranks
   MPI_Scatterv( NewParPUid_AllRank,  NSend, SendDisp,      MPI_GAMER_LONG_PAR,
                 NewParPUid_ThisRank, NNewPar_ThisRank_int, MPI_GAMER_LONG_PAR,
                 0, MPI_COMM_WORLD );

// set new particle UID in this rank
   for (long p=0; p<NNewPar_ThisRank; p++)
      amr->Par->PUid[NewParIDList[p]] = NewParPUid_ThisRank[p];

// free memory
   delete [] NewParIDList;
   delete [] NewParPUid_AllRank;
   delete [] NewParPUid_ThisRank;
   for (int d=0; d<3; d++)
   {
      delete [] NewParPos_AllRank[d];
      delete [] NewParPos_ThisRank[d];
   }

// update the next UID to all ranks
   amr->Par->NextUID += NNewPar_AllRank;

   if ( init  &&  MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID (init) ... done\n" );

   MPI_Barrier( MPI_COMM_WORLD );

} // FUNCTION : Par_SetParUID



#endif // #ifdef PARTICLE
