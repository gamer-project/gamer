#ifndef __SERIAL_H__
#define __SERIAL_H__



#include "Global.h"


// define the alternatives to all the MPI and buffer functions for the serial code
#define MPI_Barrier( MPI_COMM ) {}
#define MPI_Bcast( Buf, Count, Type, Root, MPI_COMM ) {}
#define MPI_Finalize() {}

#define MPI_Allreduce( SBuf, RBuf, Count, Type, Op, MPI_COMM ) \
        {   for (int t=0; t<(Count); t++)    (RBuf)[t] = (SBuf)[t];  }
#define MPI_Reduce( SBuf, RBuf, Count, Type, Op, Root, MPI_COMM ) \
        {   for (int t=0; t<(Count); t++)    (RBuf)[t] = (SBuf)[t];  }
#define MPI_Gather( SBuf, SCount, SType, RBuf, RCount, RType, Root, MPI_COMM ) \
        {   for (int t=0; t<(SCount); t++)   (RBuf)[t] = (SBuf)[t];  }
#define MPI_Gatherv( SBuf, SCount, SType, RBuf, RCount, RDisp, RType, Root, MPI_COMM ) \
        {   for (int t=0; t<(SCount); t++)   (RBuf)[t] = (SBuf)[t];  }
#define MPI_Alltoall( SBuf, SCount, SType, RBuf, RCount, RType, MPI_COMM ) \
        {   for (int t=0; t<(SCount); t++)   (RBuf)[t] = (SBuf)[t];  }
#define MPI_Alltoallv( SBuf, SCount, SDisp, SType, RBuf, RCount, RDisp, RType, MPI_COMM ) \
        {   for (int t=0; t<((SCount)[0]); t++)    (RBuf)[t] = (SBuf)[t];  }
#define MPI_Allgather( SBuf, SCount, SType, RBuf, RCount, RType, MPI_COMM ) \
        {   for (int t=0; t<(SCount); t++)   (RBuf)[t] = (SBuf)[t];  }
#define MPI_Allgatherv( SBuf, SCount, SType, RBuf, RCount, RDisp, RType, MPI_COMM ) \
        {   for (int t=0; t<(SCount); t++)   (RBuf)[t] = (SBuf)[t];  }

#define Buf_AllocateBufferPatch( lv ) {}
#define Buf_GetBufferData( lv, Package, ParaBuffer ) {}
#define Buf_RecordExchangeDataPatchID( lv ) {}
#define Buf_RecordExchangeFluxPatchID( lv ) {}
#define Buf_RecordBoundaryFlag( lv ) {}
#define Buf_RecordBoundaryPatch( lv ) {}
#define Buf_RecordBoundaryPatch_Base() {}
#define Buf_SortBoundaryPatch( NPatch, IDList, PosList ) {}
#define Buf_ResetBufferFlux( lv ) {}

#define MPI_ExchangeBoundaryFlag( lv ) {}
#define MPI_ExchangeBufferPosition( NSend, NRecv, Send_PosList, Recv_PosList ) {}
#define MPI_ExchangeData( TargetRank, SendSize, RecvSize, SendBuffer, RecvBuffer ) {}
#define Init_MPI( argc, argv )   { MyRank = 0; }
#define MPI_Exit() \
   {  cout << flush; fprintf( stderr, "\nProgram termination ...... rank %d\n\n", MyRank ); exit(-1);   }



#endif // #ifndef __SERIAL_H__
