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

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID ...\n" );

   long NPar_ThisRank, NPar_AllRank;
   NPar_ThisRank = amr->Par->NPar_AcPlusInac;
   int NPar_ThisRank_int = NPar_ThisRank;
   MPI_Allreduce( &NPar_ThisRank, &NPar_AllRank, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];
   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( init )
   {
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

      delete ParPIdx_AllRank;
   }
   else // if ( init )
   {
      Aux_Error( ERROR_INFO, "%s Not supported yet", __FUNCTION__ );
   } // if ( init ) ... else ...

   if ( MPI_Rank == 0 )   Aux_Message( stdout, "Par_SetParID ... done\n" );

   MPI_Barrier( MPI_COMM_WORLD );


} // FUNCTION : Par_SetParID
#endif // #ifdef PARTICLE
