#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_SendParticleData
// Description :  Exchange particles between different MPI ranks
//
// Note        :  1. SendBuf_XXX must be preallocated and will NOT be deallocated in this function
//                2. RecvBuf_XXX will be allocated in this function (using call by reference) and must be
//                   deallocated manually after calling this function
//                   --> Except for RecvBuf_ParFlt/IntDataEachPatch, which are just pointers to the MPI recv buffer
//                       declared in LB_GetBufferData
//                3. SendBuf_ParFlt/IntDataEachPatch format: [ParID][ParAttribute] instead of [ParAttribute][ParID]
//                4. Called by Par_LB_CollectParticleFromRealPatch(), Par_LB_CollectParticle2OneLevel(), and
//                   Par_LB_ExchangeParticleBetweenPatch()
//                   --> Par_LB_ExchangeParticleBetweenPatch() is called by
//                       Par_PassParticle2Sibling() and Par_PassParticle2Son_MultiPatch()
//
// Parameter   :  NParAttFlt                  : Number of particle floating-point attributes to be sent
//                NParAttInt                  : Number of particle integer        attributes to be sent
//                SendBuf_NPatchEachRank      : MPI send buffer --> number of patches sent to each rank
//                SendBuf_NParEachPatch       : MPI send buffer --> number of particles in each patch to be sent
//                SendBuf_LBIdxEachPatch      : MPI send buffer --> load-balance index of each patch to be sent
//                SendBuf_ParFltDataEachPatch : MPI send buffer --> particle floating-point data in each patch to be sent
//                SendBuf_ParIntDataEachPatch : MPI send buffer --> particle integer        data in each patch to be sent
//                NSendParTotal               : Total number of particles sent to all ranks (used by OPT__TIMING_MPI only)
//                RecvBuf_XXX                 : MPI recv buffer
//                NRecvPatchTotal             : Total number of patches   received from all ranks
//                NRecvParTotal               : Total number of particles received from all ranks
//                Exchange_NPatchEachRank     : true  : Exchange SendBuf_NPatchEachRank to get RecvBuf_NPatchEachRank
//                                              false : Assuming RecvBuf_NPatchEachRank has already been set properly
//                                                      --> But one still needs to provide SendBuf_NPatchEachRank properly
//                                                      --> RecvBuf_NPatchEachRank will NOT be reallocated
//                                                      --> Useful in Par_LB_CollectParticleFromRealPatch.cpp
//                Exchange_LBIdxEachRank      : true  : Exchange SendBuf_LBIdxEachRank to get RecvBuf_LBIdxEachRank
//                                              false : Does NOT exchange SendBuf_LBIdxEachPatch at all
//                                                      --> One does NOT need to provide SendBuf_LBIdxEachPatch
//                                                      --> RecvBuf_LBIdxEachPatch will NOT be allocated
//                                                      --> Useful in Par_LB_CollectParticleFromRealPatch.cpp
//                Exchange_ParDataEachRank    : true  : Exchange SendBuf_ParFltDataEachPatch to get RecvBuf_ParFltDataEachPatch
//                                                      Exchange SendBuf_ParIntDataEachPatch to get RecvBuf_ParIntDataEachPatch
//                Timer                       : Timer used by the options "TIMING" and "OPT__TIMING_MPI"
//                                              --> Do nothing if Timer == NULL
//                Timer_Comment               : String used by "OPT__TIMING_MPI"
//
// Return      :  RecvBuf_NPatchEachRank (if Exchange_NPatchEachRank == true), RecvBuf_NParEachPatch,
//                RecvBuf_LBIdxEachPatch (if Exchange_LBIdxEachRank == true),
//                RecvBuf_ParFltDataEachPatch (if Exchange_ParDataEachRank == true),
//                RecvBuf_ParIntDataEachPatch (if Exchange_ParDataEachRank == true),
//                NRecvPatchTotal, NRecvPatchTotal
//-------------------------------------------------------------------------------------------------------
void Par_LB_SendParticleData( const int NParAttFlt, const int NParAttInt, int *SendBuf_NPatchEachRank, int *SendBuf_NParEachPatch,
                              long *SendBuf_LBIdxEachPatch, real_par *SendBuf_ParFltDataEachPatch, long_par *SendBuf_ParIntDataEachPatch,
                              const long NSendParTotal, int *&RecvBuf_NPatchEachRank, int *&RecvBuf_NParEachPatch, long *&RecvBuf_LBIdxEachPatch,
                              real_par *&RecvBuf_ParFltDataEachPatch, long_par *&RecvBuf_ParIntDataEachPatch, int &NRecvPatchTotal,
                              long &NRecvParTotal, const bool Exchange_NPatchEachRank, const bool Exchange_LBIdxEachRank,
                              const bool Exchange_ParDataEachRank, Timer_t *Timer, const char *Timer_Comment )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( NParAttFlt < 0 )                        Aux_Error( ERROR_INFO, "NParAttFlt = %d < 0 !!\n", NParAttFlt );
   if ( NParAttInt < 0 )                        Aux_Error( ERROR_INFO, "NParAttInt = %d < 0 !!\n", NParAttInt );
   if ( SendBuf_NPatchEachRank   == NULL )      Aux_Error( ERROR_INFO, "SendBuf_NPatchEachRank == NULL !!\n" );
   if ( SendBuf_NParEachPatch    == NULL )      Aux_Error( ERROR_INFO, "SendBuf_NParEachPatch == NULL !!\n" );
   if ( Exchange_ParDataEachRank  &&
        SendBuf_ParFltDataEachPatch == NULL )   Aux_Error( ERROR_INFO, "SendBuf_ParFltDataEachPatch == NULL !!\n" );
   if ( Exchange_ParDataEachRank  &&
        SendBuf_ParIntDataEachPatch == NULL )   Aux_Error( ERROR_INFO, "SendBuf_ParIntDataEachPatch == NULL !!\n" );
   if ( Exchange_LBIdxEachRank  &&
        SendBuf_LBIdxEachPatch   == NULL )      Aux_Error( ERROR_INFO, "SendBuf_LBIdxEachPatch == NULL !!\n" );
   if ( !Exchange_NPatchEachRank  &&
        RecvBuf_NPatchEachRank   == NULL )      Aux_Error( ERROR_INFO, "RecvBuf_NParEachPatch == NULL !!\n" );
#  ifdef TIMING
   if ( Timer != NULL  &&  OPT__TIMING_MPI  &&  Timer_Comment == NULL )
      Aux_Error( ERROR_INFO, "Timer_Comment == NULL !!\n" );
#  endif
#  endif


// start timing
#  ifdef TIMING
   double time0, dtime;

   if ( Timer != NULL )
   {
//    it's better to add barrier before timing transferring data through MPI
//    --> so that the timing results (i.e., the MPI bandwidth reported by OPT__TIMING_MPI ) does NOT include
//        the time waiting for other ranks to reach here
//    --> make the MPI bandwidth measured here more accurate
      if ( OPT__TIMING_BARRIER )    MPI_Barrier( MPI_COMM_WORLD );

      if ( OPT__TIMING_MPI )  time0 = Timer->GetValue();

      Timer->Start();
   }
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
   if ( Exchange_ParDataEachRank )
   {
      long *SendCount_ParFltDataEachPatch = new long [MPI_NRank];
      long *RecvCount_ParFltDataEachPatch = new long [MPI_NRank];
      long *SendDisp_ParFltDataEachPatch  = new long [MPI_NRank];
      long *RecvDisp_ParFltDataEachPatch  = new long [MPI_NRank];
      long *SendCount_ParIntDataEachPatch = new long [MPI_NRank];
      long *RecvCount_ParIntDataEachPatch = new long [MPI_NRank];
      long *SendDisp_ParIntDataEachPatch  = new long [MPI_NRank];
      long *RecvDisp_ParIntDataEachPatch  = new long [MPI_NRank];

//    send/recv count
      const int *SendPtr = NULL, *RecvPtr = NULL;
      NRecvParTotal = 0L;

      for (int r=0; r<MPI_NRank; r++)
      {
         SendCount_ParFltDataEachPatch[r] = 0L;
         RecvCount_ParFltDataEachPatch[r] = 0L;
         SendCount_ParIntDataEachPatch[r] = 0L;
         RecvCount_ParIntDataEachPatch[r] = 0L;

         SendPtr = SendBuf_NParEachPatch + SendDisp_NParEachPatch[r];
         RecvPtr = RecvBuf_NParEachPatch + RecvDisp_NParEachPatch[r];

         for (int p=0; p<SendBuf_NPatchEachRank[r]; p++)
         {
            SendCount_ParFltDataEachPatch[r] += (long)SendPtr[p];
            SendCount_ParIntDataEachPatch[r] += (long)SendPtr[p];
         }

         for (int p=0; p<RecvBuf_NPatchEachRank[r]; p++)
         {
            RecvCount_ParFltDataEachPatch[r] += (long)RecvPtr[p];
            RecvCount_ParIntDataEachPatch[r] += (long)RecvPtr[p];
         }

         NRecvParTotal += RecvCount_ParFltDataEachPatch[r];

         SendCount_ParFltDataEachPatch[r] *= (long)NParAttFlt;
         RecvCount_ParFltDataEachPatch[r] *= (long)NParAttFlt;
         SendCount_ParIntDataEachPatch[r] *= (long)NParAttInt;
         RecvCount_ParIntDataEachPatch[r] *= (long)NParAttInt;
      }

//    send/recv displacement
      SendDisp_ParFltDataEachPatch[0] = 0L;
      RecvDisp_ParFltDataEachPatch[0] = 0L;
      SendDisp_ParIntDataEachPatch[0] = 0L;
      RecvDisp_ParIntDataEachPatch[0] = 0L;
      for (int r=1; r<MPI_NRank; r++)
      {
         SendDisp_ParFltDataEachPatch[r] = SendDisp_ParFltDataEachPatch[r-1] + SendCount_ParFltDataEachPatch[r-1];
         RecvDisp_ParFltDataEachPatch[r] = RecvDisp_ParFltDataEachPatch[r-1] + RecvCount_ParFltDataEachPatch[r-1];
         SendDisp_ParIntDataEachPatch[r] = SendDisp_ParIntDataEachPatch[r-1] + SendCount_ParIntDataEachPatch[r-1];
         RecvDisp_ParIntDataEachPatch[r] = RecvDisp_ParIntDataEachPatch[r-1] + RecvCount_ParIntDataEachPatch[r-1];
      }

//    reuse the MPI recv buffer declared in LB_GetBufferData for better MPI performance
      const long ParAllAttSize = NRecvParTotal * ( (long)NParAttFlt*sizeof(real_par) + (long)NParAttInt*sizeof(long_par) );
      RecvBuf_ParFltDataEachPatch = (real_par *)LB_GetBufferData_MemAllocate_Recv( ParAllAttSize );
      RecvBuf_ParIntDataEachPatch = (long_par *)( RecvBuf_ParFltDataEachPatch + NRecvParTotal*NParAttFlt );


//    exchange data
      MPI_Alltoallv_GAMER( SendBuf_ParFltDataEachPatch, SendCount_ParFltDataEachPatch, SendDisp_ParFltDataEachPatch, MPI_GAMER_REAL_PAR,
                           RecvBuf_ParFltDataEachPatch, RecvCount_ParFltDataEachPatch, RecvDisp_ParFltDataEachPatch, MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
      MPI_Alltoallv_GAMER( SendBuf_ParIntDataEachPatch, SendCount_ParIntDataEachPatch, SendDisp_ParIntDataEachPatch, MPI_GAMER_LONG_PAR,
                           RecvBuf_ParIntDataEachPatch, RecvCount_ParIntDataEachPatch, RecvDisp_ParIntDataEachPatch, MPI_GAMER_LONG_PAR, MPI_COMM_WORLD );

//    free memory
      delete [] SendCount_ParFltDataEachPatch;
      delete [] SendDisp_ParFltDataEachPatch;
      delete [] RecvCount_ParFltDataEachPatch;
      delete [] RecvDisp_ParFltDataEachPatch;
      delete [] SendCount_ParIntDataEachPatch;
      delete [] SendDisp_ParIntDataEachPatch;
      delete [] RecvCount_ParIntDataEachPatch;
      delete [] RecvDisp_ParIntDataEachPatch;
   } // if ( Exchange_ParDataEachRank )


// 5. free memory
   delete [] SendCount_NParEachPatch;
   delete [] SendDisp_NParEachPatch;
   delete [] RecvCount_NParEachPatch;
   delete [] RecvDisp_NParEachPatch;


// stop timing
#  ifdef TIMING
   if ( Timer != NULL )
   {
      Timer->Stop();

      if ( OPT__TIMING_MPI )
      {
         dtime = Timer->GetValue() - time0;

//       output to the same log file as LB_GetBufferData
         char FileName[2*MAX_STRING];
         sprintf( FileName, "%s/Record__TimingMPI_Rank%05d", OUTPUT_DIR, MPI_Rank );

         FILE *File = fopen( FileName, "a" );

         const double SendMB = (double)NSendParTotal*NParAttFlt*sizeof(real_par)*1.0e-6 + (double)NSendParTotal*NParAttInt*sizeof(long_par)*1.0e-6;
         const double RecvMB = (double)NRecvParTotal*NParAttFlt*sizeof(real_par)*1.0e-6 + (double)NRecvParTotal*NParAttInt*sizeof(long_par)*1.0e-6;

         fprintf( File, "%19s %2d+%1d %4s %10s %10s %10.5f %8.3f %8.3f %10.3f %10.3f\n",
                  Timer_Comment, NParAttFlt, NParAttInt, "X", "X", "X", dtime, SendMB, RecvMB, SendMB/dtime, RecvMB/dtime );

         fclose( File );
      } // if ( OPT__TIMING_MPI )
   } // if ( Timer != NULL )
#  endif // #ifdef TIMING

} // FUNCTION : Par_LB_SendParticleData



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
