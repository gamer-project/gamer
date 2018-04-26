#include "GAMER.h"

#if ( defined PARTICLE  &&  !defined SERIAL )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_Init_RedistributeByRectangular
// Description :  Redistribute particles assuming the rectangular domain decomposition
//
// Note        :  1. Usefully only during the initialization when calling Init_GAMER()
//                2. Inactive particles will NOT be redistributed
//                3. Redistribute the Time, Mass, Pos, Vel, and all passive variables of particles
//                   --> Acc is not redistributed since it has not been initialized
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass,PosX/Y/Z,VelX/Y/Z,Passive
//-------------------------------------------------------------------------------------------------------
void Par_LB_Init_RedistributeByRectangular()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   for (int d=0; d<3; d++)
      if ( MPI_NRank_X[d] < 0 )  Aux_Error( ERROR_INFO, "MPI_NRank_X[%d] = %d < 0 !!\n", d, MPI_NRank_X[d] );


   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };

   int Send_Count[MPI_NRank], Recv_Count[MPI_NRank], Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];
   int Send_Count_Sum, Recv_Count_Sum;

   for (int r=0; r<MPI_NRank; r++)  Send_Count[r] = 0;


// 1. get the target MPI rank of each particle
   const int NRankStride[3] = { 1, MPI_NRank_X[0], MPI_NRank_X[0]*MPI_NRank_X[1] };
   double (*SubDomain_EdgeR_AllRank)[3];
   int    TRank_3D[3];
   int   *TRank = new int [amr->Par->NPar_AcPlusInac];

   SubDomain_EdgeR_AllRank = new double [MPI_NRank][3];

   MPI_Allgather( amr->ParaVar->SubDomain_EdgeR, 3, MPI_DOUBLE, SubDomain_EdgeR_AllRank, 3, MPI_DOUBLE, MPI_COMM_WORLD );

   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    TRank is set to -1 for inactive particles
      if ( Mass[ParID] < (real)0.0 )
      {
         TRank[ParID] = -1;
         continue;
      }

      for (int d=0; d<3; d++)
      {
         TRank_3D[d] = -1;

         for (int r=0; r<MPI_NRank_X[d]; r++)
         {
            if ( Pos[d][ParID] < SubDomain_EdgeR_AllRank[ r*NRankStride[d] ][d] )
            {
               TRank_3D[d] = r;
               break;
            }
         }

         if ( TRank_3D[d] == -1 )
            Aux_Error( ERROR_INFO, "cannot determine target rank for ParID %ld with Pos = (%20.14e, %20.14e, %20.14e) !!\n",
                       ParID, Pos[0][ParID], Pos[1][ParID], Pos[2][ParID] );

         if ( Pos[d][ParID] < (real)0.0 )
            Aux_Error( ERROR_INFO, "Pos[%d][%ld] = %21.14e < 0.0 !!\n", d, ParID, Pos[d][ParID] );
      } // for (int d=0; d<3; d++)

      TRank[ParID] = (int)Mis_Idx3D2Idx1D( MPI_NRank_X, TRank_3D );

      Send_Count[ TRank[ParID] ] ++;
   } // for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)


// 2. construct the MPI send and recv data list
   MPI_Alltoall( Send_Count, 1, MPI_INT, Recv_Count, 1, MPI_INT, MPI_COMM_WORLD );

   Send_Disp[0] = 0;
   Recv_Disp[0] = 0;

   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + Send_Count[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + Recv_Count[r-1];
   }

   Send_Count_Sum = Send_Disp[ MPI_NRank-1 ] + Send_Count[ MPI_NRank-1 ];
   Recv_Count_Sum = Recv_Disp[ MPI_NRank-1 ] + Recv_Count[ MPI_NRank-1 ];

   if ( Send_Count_Sum != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "Send_Count_Sum (%d) != NPar_Active (%ld) !!\n", Send_Count_Sum, amr->Par->NPar_Active );


// 3. redistribute particle attributes (one attribute at a time to save memory)
   const int NAtt = 8 + PAR_NPASSIVE;  // 8 = mass, pos*3, vel*3, time --> do not transfer acceleration
   int Offset[MPI_NRank];

   real *SendBuf = new real [Send_Count_Sum];
   real *RecvBuf = NULL;

// 3-1. record attribute pointers
   real *AttPtr[NAtt], **AttPtrPtr[NAtt];

   AttPtr[0] = amr->Par->ParVar[PAR_MASS];
   AttPtr[1] = amr->Par->ParVar[PAR_POSX];
   AttPtr[2] = amr->Par->ParVar[PAR_POSY];
   AttPtr[3] = amr->Par->ParVar[PAR_POSZ];
   AttPtr[4] = amr->Par->ParVar[PAR_VELX];
   AttPtr[5] = amr->Par->ParVar[PAR_VELY];
   AttPtr[6] = amr->Par->ParVar[PAR_VELZ];
   AttPtr[7] = amr->Par->ParVar[PAR_TIME];

   for (int v=8; v<NAtt; v++)    AttPtr[v] = amr->Par->Passive[v-8];

   AttPtrPtr[0] = &amr->Par->ParVar[PAR_MASS];
   AttPtrPtr[1] = &amr->Par->ParVar[PAR_POSX];
   AttPtrPtr[2] = &amr->Par->ParVar[PAR_POSY];
   AttPtrPtr[3] = &amr->Par->ParVar[PAR_POSZ];
   AttPtrPtr[4] = &amr->Par->ParVar[PAR_VELX];
   AttPtrPtr[5] = &amr->Par->ParVar[PAR_VELY];
   AttPtrPtr[6] = &amr->Par->ParVar[PAR_VELZ];
   AttPtrPtr[7] = &amr->Par->ParVar[PAR_TIME];

   for (int v=8; v<NAtt; v++)    AttPtrPtr[v] = &amr->Par->Passive[v-8];

   for (int v=0; v<NAtt; v++)
   {
//    3-2. initialize the array offsets of all target ranks
      for (int r=0; r<MPI_NRank; r++)  Offset[r] = Send_Disp[r];

//    3-3. prepare send buffer
      for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
         if ( TRank[ParID] != -1 )  SendBuf[ Offset[TRank[ParID]] ++ ] = AttPtr[v][ParID];

//    3-4. free the old particle array
      free( AttPtr[v] );

//    3-5. allocate the new particle array and set the recv buffer
      *(AttPtrPtr[v]) = (real*)malloc( Recv_Count_Sum*sizeof(real) );
      RecvBuf         = *(AttPtrPtr[v]);

//    3-6. redistribute data
#     ifdef FLOAT8
      MPI_Alltoallv( SendBuf, Send_Count, Send_Disp, MPI_DOUBLE, RecvBuf, Recv_Count, Recv_Disp, MPI_DOUBLE, MPI_COMM_WORLD );
#     else
      MPI_Alltoallv( SendBuf, Send_Count, Send_Disp, MPI_FLOAT,  RecvBuf, Recv_Count, Recv_Disp, MPI_FLOAT,  MPI_COMM_WORLD );
#     endif
   } // for (int t=0; t<NAtt; t++)


// 4. reset particle parameters
   amr->Par->NPar_AcPlusInac     = (long)Recv_Count_Sum;
   amr->Par->NPar_Active         = (long)Recv_Count_Sum;
   amr->Par->NPar_Inactive       = 0;                       // since we don't redistribute inactive particles
   amr->Par->ParListSize         = (long)Recv_Count_Sum;
   amr->Par->InactiveParListSize = (long)MAX( 1, Recv_Count_Sum/100 );

   if ( amr->Par->InactiveParList != NULL )
   {
      free( amr->Par->InactiveParList );

      amr->Par->InactiveParList = (long*)malloc( amr->Par->InactiveParListSize*sizeof(long) );
   }


// 5. reset attribute pointers and reallocate the acceleration arrays
   amr->Par->Mass = amr->Par->ParVar[PAR_MASS];
   amr->Par->PosX = amr->Par->ParVar[PAR_POSX];
   amr->Par->PosY = amr->Par->ParVar[PAR_POSY];
   amr->Par->PosZ = amr->Par->ParVar[PAR_POSZ];
   amr->Par->VelX = amr->Par->ParVar[PAR_VELX];
   amr->Par->VelY = amr->Par->ParVar[PAR_VELY];
   amr->Par->VelZ = amr->Par->ParVar[PAR_VELZ];
   amr->Par->Time = amr->Par->ParVar[PAR_TIME];

#  ifdef STORE_PAR_ACC
   free( amr->Par->ParVar[PAR_ACCX] );
   free( amr->Par->ParVar[PAR_ACCY] );
   free( amr->Par->ParVar[PAR_ACCZ] );

   amr->Par->ParVar[PAR_ACCX] = (real*)malloc( amr->Par->ParListSize*sizeof(real) );
   amr->Par->ParVar[PAR_ACCY] = (real*)malloc( amr->Par->ParListSize*sizeof(real) );
   amr->Par->ParVar[PAR_ACCZ] = (real*)malloc( amr->Par->ParListSize*sizeof(real) );

   amr->Par->AccX = amr->Par->ParVar[PAR_ACCX];
   amr->Par->AccY = amr->Par->ParVar[PAR_ACCY];
   amr->Par->AccZ = amr->Par->ParVar[PAR_ACCZ];
#  endif // #ifdef STORE_PAR_ACC


// 6. check
   Mass   = amr->Par->Mass;
   Pos[0] = amr->Par->PosX;
   Pos[1] = amr->Par->PosY;
   Pos[2] = amr->Par->PosZ;

   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    there should be no inactive particles
      if ( Mass[ParID] < 0.0 )   Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, Mass[ParID] );

      for (int d=0; d<3; d++)
         if ( Pos[d][ParID] < amr->ParaVar->SubDomain_EdgeL[d]  ||  Pos[d][ParID] >= amr->ParaVar->SubDomain_EdgeR[d] )
            Aux_Error( ERROR_INFO, "Pos[%d][%ld] = %20.14e lies outside the sub-domain [%20.14e ... %20.14e] !!\n",
                       d, ParID, Pos[d][ParID], amr->ParaVar->SubDomain_EdgeL[d], amr->ParaVar->SubDomain_EdgeR[d] );
   }


// free memory
   delete [] TRank;
   delete [] SendBuf;
   delete [] SubDomain_EdgeR_AllRank;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_LB_Init_RedistributeByRectangular



#endif // #if ( defined PARTICLE  &&  !defined SERIAL )
