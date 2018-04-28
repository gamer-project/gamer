#ifndef __SERIAL_H__
#define __SERIAL_H__



#include "Global.h"


// define the alternatives to all the MPI and buffer functions for the serial code
#define MPI_Barrier( MPI_COMM ) {}
#define MPI_Bcast( Buf, Count, Type, Root, MPI_COMM ) {}
#define MPI_Finalize() {}

// make sure to use a unique variable name for the loop counter (i.e., ssss below)
#define MPI_Allreduce( SBuf, RBuf, Count, Type, Op, MPI_COMM ) \
        {   for (int ssss=0; ssss<(Count); ssss++)    (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Reduce( SBuf, RBuf, Count, Type, Op, Root, MPI_COMM ) \
        {   for (int ssss=0; ssss<(Count); ssss++)    (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Gather( SBuf, SCount, SType, RBuf, RCount, RType, Root, MPI_COMM ) \
        {   for (int ssss=0; ssss<(SCount); ssss++)   (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Gatherv( SBuf, SCount, SType, RBuf, RCount, RDisp, RType, Root, MPI_COMM ) \
        {   for (int ssss=0; ssss<(SCount); ssss++)   (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Alltoall( SBuf, SCount, SType, RBuf, RCount, RType, MPI_COMM ) \
        {   for (int ssss=0; ssss<(SCount); ssss++)   (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Alltoallv( SBuf, SCount, SDisp, SType, RBuf, RCount, RDisp, RType, MPI_COMM ) \
        {   for (int ssss=0; ssss<((SCount)[0]); ssss++)    (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Allgather( SBuf, SCount, SType, RBuf, RCount, RType, MPI_COMM ) \
        {   for (int ssss=0; ssss<(SCount); ssss++)   (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Allgatherv( SBuf, SCount, SType, RBuf, RCount, RDisp, RType, MPI_COMM ) \
        {   for (int ssss=0; ssss<(SCount); ssss++)   (RBuf)[ssss] = (SBuf)[ssss];  }
#define MPI_Scatterv( SBuf, SCount, SDisp, SType, RBuf, RCount, RType, Root, MPI_COMM ) \
        {   for (int ssss=0; ssss<((SCount)[0]); ssss++)    (RBuf)[ssss] = (SBuf)[ssss];  }


#define Buf_AllocateBufferPatch( Tpatch, lv ) {}
#define Buf_GetBufferData( lv, FluSg, PotSg, GetBufMode, TVar, ParaBuf, UseLBFunc ) {}
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
#define MPI_CubeSlice( Dir, SendBuf, RecvBuf ) {}
#define Init_MPI( argc, argv )   { MPI_Rank = 0; }
#define MPI_Exit() \
   {  fflush( stdout ); fflush( stderr ); fprintf( stderr, "\nProgram termination ...... rank %d\n\n", MPI_Rank ); exit( 1 );   }
#define Flag_Buffer( lv ) {}
#define Refine_Buffer( lv, SonTable, GrandTable ) {}
#define Flu_AllocateFluxArray_Buffer( lv ) {}



#endif // #ifndef __SERIAL_H__
