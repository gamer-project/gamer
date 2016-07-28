#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_SendParticle
// Description :  Send particles to other MPI ranks
//
// Note        :  1. SendBuf_XXX must be preallocated and will NOT be deallocated in this function
//                2. RecvBuf_XXX will be allocated in this function (using call by reference) and must be
//                   deallocated manually after calling this function
//                3. SendBuf_ParDataEachPatch format: [ParID][ParAttribute] instead of [ParAttribute][ParID]
//
// Parameter   :  NParVar                  : Number of particle attributes to be sent
//                SendBuf_NPatchEachRank   : MPI send buffer --> number of patches sent to each rank
//                SendBuf_NParEachPatch    : MPI send buffer --> number of particles in each patch to be sent
//                SendBuf_LBIdxEachPatch   : MPI send buffer --> load-balance index of each patch to be sent
//                SendBuf_ParDataEachPatch : MPI send buffer --> particle data in each patch to be sent
//                RecvBuf_XXX              : MPI recv buffer
//                NRecvPatchTotal          : Total number of patches   received from all ranks
//                NRecvParTotal            : Total number of particles received from all ranks
//                Exchange_NPatchEachRank  : true  : Exchange SendBuf_NPatchEachRank to get RecvBuf_NPatchEachRank
//                                           false : Assuming RecvBuf_NPatchEachRank has already been set properly
//                                                   --> But one still needs to provide SendBuf_NPatchEachRank properly
//                                                   --> RecvBuf_NPatchEachRank will NOT be reallocated
//                                                   --> Useful in Par_LB_CollectParticleFromRealPatch.cpp
//                Exchange_LBIdxEachRank   : true  : Exchange SendBuf_LBIdxEachRank to get RecvBuf_LBIdxEachRank
//                                           false : Does NOT exchange SendBuf_LBIdxEachPatch at all
//                                                   --> One does NOT need to provide SendBuf_LBIdxEachPatch
//                                                   --> RecvBuf_LBIdxEachPatch will NOT be allocated
//                                                   --> Useful in Par_LB_CollectParticleFromRealPatch.cpp
//
// Return      :  RecvBuf_NPatchEachRank (if Exchange_NPatchEachRank == true), RecvBuf_NParEachPatch,
//                RecvBuf_LBIdxEachPatch (if Exchange_LBIdxEachRank == true), RecvBuf_ParDataEachPatch,
//                NRecvPatchTotal, NRecvPatchTotal
//-------------------------------------------------------------------------------------------------------
void Par_LB_SendParticle( const int NParVar, const int *SendBuf_NPatchEachRank,
                          const int *SendBuf_NParEachPatch, const long *SendBuf_LBIdxEachPatch, const real *SendBuf_ParDataEachPatch,
                          int *&RecvBuf_NPatchEachRank, int *&RecvBuf_NParEachPatch, long *&RecvBuf_LBIdxEachPatch,
                          real *&RecvBuf_ParDataEachPatch, int &NRecvPatchTotal, int &NRecvParTotal,
                          const bool Exchange_NPatchEachRank, const bool Exchange_LBIdxEachRank )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( NParVar < 0 )                        Aux_Error( ERROR_INFO, "NParVar = %d < 0 !!\n", NParVar );
   if ( SendBuf_NPatchEachRank   == NULL )   Aux_Error( ERROR_INFO, "SendBuf_NPatchEachRank == NULL !!\n" );
   if ( SendBuf_NParEachPatch    == NULL )   Aux_Error( ERROR_INFO, "SendBuf_NParEachPatch == NULL !!\n" );
   if ( SendBuf_ParDataEachPatch == NULL )   Aux_Error( ERROR_INFO, "SendBuf_ParDataEachPatch == NULL !!\n" );
   if ( Exchange_LBIdxEachRank  &&
        SendBuf_LBIdxEachPatch   == NULL )   Aux_Error( ERROR_INFO, "SendBuf_LBIdxEachPatch == NULL !!\n" );
   if ( !Exchange_NPatchEachRank  &&
        RecvBuf_NPatchEachRank   == NULL )   Aux_Error( ERROR_INFO, "RecvBuf_NParEachPatch == NULL !!\n" );
#  endif


// 1. get the number of patches received from each rank
   if ( Exchange_NPatchEachRank )
   {
      RecvBuf_NPatchEachRank = new int [MPI_NRank];

      MPI_Alltoall( SendBuf_NPatchEachRank, 1, MPI_INT, RecvBuf_NPatchEachRank, 1, MPI_INT, MPI_COMM_WORLD );
   }

   NRecvPatchTotal = 0;
   for (int r=0; r<MPI_NRank; r++)  NRecvPatchTotal += RecvBuf_NPatchEachRank[r];


// 2. get the number of particles received from each rank
   int *SendCount_NParEachPatch = new int [MPI_NRank];
   int *RecvCount_NParEachPatch = new int [MPI_NRank];
   int *SendDisp_NParEachPatch  = new int [MPI_NRank];
   int *RecvDisp_NParEachPatch  = new int [MPI_NRank];

   RecvBuf_NParEachPatch = new int [NRecvPatchTotal];

// send/recv count
   for (int r=0; r<MPI_NRank; r++)
   {
      SendCount_NParEachPatch[r] = SendBuf_NPatchEachRank[r];
      RecvCount_NParEachPatch[r] = RecvBuf_NPatchEachRank[r];
   }

// send/recv displacement
   SendDisp_NParEachPatch[0] = 0;
   RecvDisp_NParEachPatch[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      SendDisp_NParEachPatch[r] = SendDisp_NParEachPatch[r-1] + SendCount_NParEachPatch[r-1];
      RecvDisp_NParEachPatch[r] = RecvDisp_NParEachPatch[r-1] + RecvCount_NParEachPatch[r-1];
   }

// exchange data
   MPI_Alltoallv( SendBuf_NParEachPatch, SendCount_NParEachPatch, SendDisp_NParEachPatch, MPI_INT,
                  RecvBuf_NParEachPatch, RecvCount_NParEachPatch, RecvDisp_NParEachPatch, MPI_INT, MPI_COMM_WORLD );


// 3. collect LBIdx from all ranks
   if ( Exchange_LBIdxEachRank )
   {
      RecvBuf_LBIdxEachPatch = new long [NRecvPatchTotal];

      MPI_Alltoallv( SendBuf_LBIdxEachPatch, SendCount_NParEachPatch, SendDisp_NParEachPatch, MPI_LONG,
                     RecvBuf_LBIdxEachPatch, RecvCount_NParEachPatch, RecvDisp_NParEachPatch, MPI_LONG, MPI_COMM_WORLD );
   }


// 4. collect particle attributes from all ranks
   int *SendCount_ParDataEachPatch = new int [MPI_NRank];
   int *RecvCount_ParDataEachPatch = new int [MPI_NRank];
   int *SendDisp_ParDataEachPatch  = new int [MPI_NRank];
   int *RecvDisp_ParDataEachPatch  = new int [MPI_NRank];

// send/recv count
   const int *SendPtr = NULL, *RecvPtr = NULL;
   NRecvParTotal = 0;

   for (int r=0; r<MPI_NRank; r++)
   {
      SendCount_ParDataEachPatch[r] = 0;
      RecvCount_ParDataEachPatch[r] = 0;

      SendPtr = SendBuf_NParEachPatch + SendDisp_NParEachPatch[r];
      RecvPtr = RecvBuf_NParEachPatch + RecvDisp_NParEachPatch[r];

      for (int p=0; p<SendBuf_NPatchEachRank[r]; p++)    SendCount_ParDataEachPatch[r] += SendPtr[p];
      for (int p=0; p<RecvBuf_NPatchEachRank[r]; p++)    RecvCount_ParDataEachPatch[r] += RecvPtr[p];

      NRecvParTotal += RecvCount_ParDataEachPatch[r];

      SendCount_ParDataEachPatch[r] *= NParVar;
      RecvCount_ParDataEachPatch[r] *= NParVar;
   }

// send/recv displacement
   SendDisp_ParDataEachPatch[0] = 0;
   RecvDisp_ParDataEachPatch[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      SendDisp_ParDataEachPatch[r] = SendDisp_ParDataEachPatch[r-1] + SendCount_ParDataEachPatch[r-1];
      RecvDisp_ParDataEachPatch[r] = RecvDisp_ParDataEachPatch[r-1] + RecvCount_ParDataEachPatch[r-1];
   }

// exchange data
   RecvBuf_ParDataEachPatch = new real [ NRecvParTotal*NParVar ];

#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_ParDataEachPatch, SendCount_ParDataEachPatch, SendDisp_ParDataEachPatch, MPI_DOUBLE,
                  RecvBuf_ParDataEachPatch, RecvCount_ParDataEachPatch, RecvDisp_ParDataEachPatch, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_ParDataEachPatch, SendCount_ParDataEachPatch, SendDisp_ParDataEachPatch, MPI_FLOAT,
                  RecvBuf_ParDataEachPatch, RecvCount_ParDataEachPatch, RecvDisp_ParDataEachPatch, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif


// 5. free memory
   delete [] SendCount_NParEachPatch;
   delete [] SendDisp_NParEachPatch;
   delete [] RecvCount_NParEachPatch;
   delete [] RecvDisp_NParEachPatch;

   delete [] SendCount_ParDataEachPatch;
   delete [] SendDisp_ParDataEachPatch;
   delete [] RecvCount_ParDataEachPatch;
   delete [] RecvDisp_ParDataEachPatch;

} // FUNCTION : Par_LB_SendParticle



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
