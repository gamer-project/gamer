#include "GAMER.h"

#ifdef PARTICLE

static long ParPos2LBIdx( const int lv, const real ParPos[] );
static void SendParticle2HomeRank( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_FindHomePatch_UniformGrid
// Description :  Find the home patches of all particles on a fully occupied level
//
// Note        :  1. The target level must be fully occupied by patches
//                   --> So either "lv=0" or "lv-1 is fully refined"
//                2. Invoked by Init_UniformGrid() and Init_ByFile()
//                3. After calling this function, the amr->Par structure will be reconstructed and
//                   all particles will be associated with their home patches on lv
//
// Parameter   :  lv : Target level
//
// Return      :  1. amr->Par
//                2. NPar, ParListSize and ParList[] of all real patches on lv
//-------------------------------------------------------------------------------------------------------
void Par_FindHomePatch_UniformGrid( const int lv )
{

// check
#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv > TOP_LEVEL )   Aux_Error( ERROR_INFO, "incorrect lv = %d !!\n", lv );
#  endif


// 1. redistribute all particles to their home ranks
   SendParticle2HomeRank( lv );


// 2. find the home patch
//    --> note that SendParticle2HomeRank() should already remove all inactive particles (i.e., those with Mass<0.0)
//    --> so in the following we can assume that all particles are active and their home patches exist
#  ifdef DEBUG_PARTICLE
   if ( amr->Par->NPar_Active != amr->Par->NPar_AcPlusInac )
      Aux_Error( ERROR_INFO, "NPar_Active (%ld) != NPar_AcPlusInac (%ld) !!\n",
                 amr->Par->NPar_Active, amr->Par->NPar_AcPlusInac );
#  endif

   const real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   const int   NReal  = amr->NPatchComma[lv][1];

   real TParPos[3];

   long *HomeLBIdx          = new long [amr->Par->NPar_Active];
   int  *HomeLBIdx_IdxTable = new int  [amr->Par->NPar_Active];
   int  *HomePID            = new int  [amr->Par->NPar_Active];
   int  *MatchIdx           = new int  [amr->Par->NPar_Active];

// get the load-balance index of the particle's home patch
   for (long ParID=0; ParID<amr->Par->NPar_Active; ParID++)
   {
      for (int d=0; d<3; d++)    TParPos[d] = Pos[d][ParID];
      HomeLBIdx[ParID] = ParPos2LBIdx( lv, TParPos );
   }

// sort the LBIdx list
   Mis_Heapsort( amr->Par->NPar_Active, HomeLBIdx, HomeLBIdx_IdxTable );

// construct and sort the LBIdx list of all real patches
// --> do not use amr->LB->IdxList_Real[] since it may not be constructed yet
   long *RealPatchLBIdx          = new long [NReal];
   int  *RealPatchLBIdx_IdxTable = new int  [NReal];

   for (int PID=0; PID<NReal; PID++)   RealPatchLBIdx[PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NReal, RealPatchLBIdx, RealPatchLBIdx_IdxTable );

// LBIdx --> home (real) patch indices
   Mis_Matching_int( NReal, RealPatchLBIdx, amr->Par->NPar_Active, HomeLBIdx, MatchIdx );

// check: every particle must have a home patch
#  ifdef DEBUG_PARTICLE
   for (long t=0; t<amr->Par->NPar_Active; t++)
   {
      if ( MatchIdx[t] == -1 )
      {
         const long ParID = HomeLBIdx_IdxTable[t];

         Aux_Error( ERROR_INFO, "lv %d, ParID %ld, ParPos (%14.7e, %14.7e, %14.7e) --> found no home patch !!\n",
                    lv, ParID, Pos[0][ParID], Pos[1][ParID], Pos[2][ParID] );
      }
   }
#  endif

// record the home patch indices
   for (long t=0; t<amr->Par->NPar_Active; t++)
   {
      const long ParID = HomeLBIdx_IdxTable[t];
      const int  PID   = RealPatchLBIdx_IdxTable[ MatchIdx[t] ];

      HomePID[ParID] = PID;
   }


// 3. count the number of particles in each patch and allocate the particle list
   for (int PID=0; PID<NReal; PID++)   amr->patch[0][lv][PID]->ParListSize = 0;

   for (long ParID=0; ParID<amr->Par->NPar_Active; ParID++)
      amr->patch[0][lv][ HomePID[ParID] ]->ParListSize ++;

   for (int PID=0; PID<NReal; PID++)
   {
      if ( amr->patch[0][lv][PID]->ParList != NULL )
      {
         free( amr->patch[0][lv][PID]->ParList );
         amr->patch[0][lv][PID]->ParList = NULL;
      }

      if ( amr->patch[0][lv][PID]->ParListSize > 0 )
      {
         amr->patch[0][lv][PID]->ParList = (long*)malloc( amr->patch[0][lv][PID]->ParListSize*sizeof(long) );
      }
   }


// 4. associate particles with their home patches
   amr->Par->NPar_Lv[lv] = 0;

// no OpenMP since AddParticle() will modify amr->Par->NPar_Lv[]
   for (long ParID=0; ParID<amr->Par->NPar_Active; ParID++)
   {
      const int PID = HomePID[ParID];

#     ifdef DEBUG_PARTICLE
      amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv],
                                           Pos, amr->Par->NPar_Active, __FUNCTION__ );
#     else
      amr->patch[0][lv][PID]->AddParticle( 1, &ParID, &amr->Par->NPar_Lv[lv] );
#     endif
   }


   delete [] HomeLBIdx;
   delete [] HomeLBIdx_IdxTable;
   delete [] HomePID;
   delete [] MatchIdx;
   delete [] RealPatchLBIdx;
   delete [] RealPatchLBIdx_IdxTable;

} // FUNCTION : Par_FindHomePatch_UniformGrid



//-------------------------------------------------------------------------------------------------------
// Function    :  SendParticle2HomeRank
// Description :  Redistribute all particles to their home ranks
//
// Note        :  1. amr->LB->CutPoint[lv][] must be set in advance (e.g., by calling LB_SetCutPoint())
//                2. Invoked by Par_FindHomePatch_UniformGrid()
//                3. Inactive particles will NOT be redistributed
//                4. Currently only redistribute the Time, Mass, Pos, Vel, and all passive variables
//                   of particles
//                   --> Acc is not redistributed since it has not been initialized when calling
//                       Par_FindHomePatch_UniformGrid()
//
// Parameter   :  lv : Target level
//
// Return      :  1. amr->Par
//                2. NPar, ParListSize and ParList[] of all real patches on lv
//-------------------------------------------------------------------------------------------------------
void SendParticle2HomeRank( const int lv )
{

   real *Mass   =   amr->Par->Mass;
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };

   int Send_Count[MPI_NRank], Recv_Count[MPI_NRank], Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];
   int Send_Count_Sum, Recv_Count_Sum;

   for (int r=0; r<MPI_NRank; r++)  Send_Count[r] = 0;


// 1. get the target MPI rank of each particle
   int *TRank = new int [amr->Par->NPar_AcPlusInac];
   real TParPos[3];
   long LBIdx;

   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    TRank is set to -1 for inactive particles
      if ( Mass[ParID] < (real)0.0 )
      {
         TRank[ParID] = -1;
         continue;
      }

//    get the load-balance index of the particle's home patch
      for (int d=0; d<3; d++)    TParPos[d] = Pos[d][ParID];
      LBIdx = ParPos2LBIdx( lv, TParPos );

//    record the home patch
#     ifdef SERIAL
      TRank[ParID] = 0;
#     else
      TRank[ParID] = LB_Index2Rank( lv, LBIdx, CHECK_ON );
#     endif
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
   } // for (int v=0; v<NAtt; v++)


// 4. reset particle parameters
   if ( amr->Par->InactiveParList != NULL )  free( amr->Par->InactiveParList );

   amr->Par->NPar_AcPlusInac     = (long)Recv_Count_Sum;
   amr->Par->NPar_Active         = (long)Recv_Count_Sum;
   amr->Par->NPar_Inactive       = 0;                       // since we don't redistribute inactive particles
   amr->Par->ParListSize         = (long)Recv_Count_Sum;
   amr->Par->InactiveParListSize = (long)MAX( 1, Recv_Count_Sum/100 );
   amr->Par->InactiveParList     = (long*)malloc( amr->Par->InactiveParListSize*sizeof(long) );


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
#  ifdef DEBUG_PARTICLE
   for (long ParID=0; ParID<amr->Par->NPar_AcPlusInac; ParID++)
   {
//    there should be no inactive particles
      if ( amr->Par->Mass[ParID] < 0.0 )
         Aux_Error( ERROR_INFO, "ParMass[%ld] = %14.7e < 0.0 !!\n", ParID, amr->Par->Mass[ParID] );
   }
#  endif


// free memory
   delete [] TRank;
   delete [] SendBuf;

} // FUNCTION : SendParticle2HomeRank



//-------------------------------------------------------------------------------------------------------
// Function    :  ParPos2LBIdx
// Description :  Convert the input particle position to the load-balance index (LBIdx) of its home patch 
//
// Note        :  1. Input particle position must lie within the simulation box
//                   --> Does not check periodicity
//                2. Invoked by SendParticle2HomeRank() and Par_FindHomePatch_UniformGrid()
//
// Parameter   :  lv     : Target level
//                ParPos : Particle position
//
// Return      :  Load-balance index (LBIdx)
//-------------------------------------------------------------------------------------------------------
long ParPos2LBIdx( const int lv, const real ParPos[] )
{

// check
#  ifdef DEBUG_PARTICLE
   for (int d=0; d<3; d++)
   {
      if ( ParPos[d] < amr->BoxEdgeL[d]  ||  ParPos[d] >= amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "ParPos[%d] = %14.7e lies outside the simulation box (BoxEdgeL/R = %14.7e/%14.7e) !!\n",
                    d, ParPos[d], amr->BoxEdgeL[d], amr->BoxEdgeR[d] );
   }
#  endif


// calculate the home patch corner
   const double dh_min        = amr->dh[TOP_LEVEL];
   const double _PatchPhySize = 1.0/( PS1*amr->dh[lv] );
   const int    PatchScale    = PS1*amr->scale[lv];

   double PatchEdgeL, PatchEdgeR;
   int    Cr[3];

   for (int d=0; d<3; d++)
   {
      Cr[d] = (int)floor(  ( (double)ParPos[d] - amr->BoxEdgeL[d] )*_PatchPhySize  )*PatchScale;

//    check the home patch corner carefully to prevent from any issue resulting from round-off errors
//    --> make sure to adopt the same procedure of calculating the patch left/right edges as Patch.h
      PatchEdgeL = (double)( Cr[d]              )*dh_min;
      PatchEdgeR = (double)( Cr[d] + PatchScale )*dh_min;

      if      ( ParPos[d] <  PatchEdgeL )    Cr[d] -= PatchScale;
      else if ( ParPos[d] >= PatchEdgeR )    Cr[d] += PatchScale;

//    check if the target home patch does enclose the given particle position
#     ifdef DEBUG_PARTICLE
      PatchEdgeL = (double)( Cr[d]              )*dh_min;
      PatchEdgeR = (double)( Cr[d] + PatchScale )*dh_min;

      if ( ParPos[d] < PatchEdgeL  ||  ParPos[d] >= PatchEdgeR )
         Aux_Error( ERROR_INFO, "incorrect home patch (dim %d, ParPos %14.7e, PatchEdgeL/R %14.7e/%14.7e) !!\n",
                    d, ParPos[d], PatchEdgeL, PatchEdgeR );

      if ( PatchEdgeL < amr->BoxEdgeL[d]  ||  PatchEdgeR > amr->BoxEdgeR[d] )
         Aux_Error( ERROR_INFO, "incorrect home patch (dim %d, PatchEdgeL/R %14.7e/%14.7e, BoxEdgeL/R %14.7e/%14.7e) !!\n",
                    d, PatchEdgeL, PatchEdgeR, amr->BoxEdgeL[d], amr->BoxEdgeR[d] );
#     endif
   } // for (int d=0; d<3; d++)


// corner --> LBIdx and return
   const long LBIdx = LB_Corner2Index( lv, Cr, CHECK_ON );

   return LBIdx;

} // FUNCTION : ParPos2LBIdx



#endif // #ifdef PARTICLE
