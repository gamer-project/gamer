#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

extern char    ClusterMerger_File_Par1[1000];
extern char    ClusterMerger_File_Par2[1000];
extern bool    ClusterMerger_Coll;
extern double  ClusterMerger_Coll_D;
extern double  ClusterMerger_Coll_B;
extern double  ClusterMerger_Coll_BulkVel1;
extern double  ClusterMerger_Coll_BulkVel2;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize all particle attributes
//
// Note        :  1. Invoked by "Init_GAMER"
//                2. Assuming all particles are active initially
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass/Pos(X/Y/Z)/Vel(X/Y/Z)/Time
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int NParVar = 7;  // mass, pos*3, vel*3
   double (*ParData_AllRank)[NParVar] = NULL;

// load data --> for simplicity, currently only the root rank will load data from disk
   if ( MPI_Rank == 0 )
   {
//    check file existence
      if ( !Aux_CheckFileExist(ClusterMerger_File_Par1) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", ClusterMerger_File_Par1 );

      if ( ClusterMerger_Coll  &&  !Aux_CheckFileExist(ClusterMerger_File_Par2) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", ClusterMerger_File_Par2 );


//    check file size
      FILE *File=NULL;
      long  NPar1=0, NPar2=0, NParAll;

      File  = fopen( ClusterMerger_File_Par1, "rb" );
      fseek( File, 0, SEEK_END );
      NPar1 = ftell( File ) / ( long(NParVar)*sizeof(double) );
      fclose( File );
      Aux_Message( stdout, "   Number of particles in cluster 1 = %ld\n", NPar1 );

      if ( ClusterMerger_Coll ) {
      File  = fopen( ClusterMerger_File_Par2, "rb" );
      fseek( File, 0, SEEK_END );
      NPar2 = ftell( File ) / ( long(NParVar)*sizeof(double) );
      fclose( File );
      Aux_Message( stdout, "   Number of particles in cluster 2 = %ld\n", NPar2 );
      }

      NParAll = NPar1 + NPar2;

      if ( NParAll != amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                    NParAll, amr->Par->NPar_Active_AllRank );


//    load data
      Aux_Message( stdout, "   Load particles from disk ... " );

      ParData_AllRank = new double [NParAll][NParVar];

//    note that fread may fail for large files if sizeof(size_t) == 4 instead of 8
      File = fopen( ClusterMerger_File_Par1, "rb" );
      fread( ParData_AllRank,       sizeof(double), long(NParVar)*NPar1, File );
      fclose( File );

      if ( ClusterMerger_Coll ) {
      File = fopen( ClusterMerger_File_Par2, "rb" );
      fread( ParData_AllRank+NPar1, sizeof(double), long(NParVar)*NPar2, File );
      fclose( File );
      }

      Aux_Message( stdout, "done\n" );


//    convert to code units (assuming loaded data are in cgs)
//    --> assuming particle data are stored in the order [ID][Mass/PosX/Y/Z/VelX/Y/Z]
      Aux_Message( stdout, "   Convert particle attributes to code units ... " );

      for (long p=0; p<NParAll; p++)
      {
         ParData_AllRank[p][0] /= UNIT_M;

         ParData_AllRank[p][1] /= UNIT_L;
         ParData_AllRank[p][2] /= UNIT_L;
         ParData_AllRank[p][3] /= UNIT_L;

         ParData_AllRank[p][4] /= UNIT_V;
         ParData_AllRank[p][5] /= UNIT_V;
         ParData_AllRank[p][6] /= UNIT_V;
      }

      Aux_Message( stdout, "done\n" );


//    shift center (assuming the center of loaded particles = [0,0,0])
      Aux_Message( stdout, "   Shift particle center and add bulk velocity ... " );
      const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

      if ( ClusterMerger_Coll )
      {
         const double ClusterCenter1[3]
            = { BoxCenter[0]-0.5*ClusterMerger_Coll_D, BoxCenter[1]-0.5*ClusterMerger_Coll_B, BoxCenter[2] };
         const double ClusterCenter2[3]
            = { BoxCenter[0]+0.5*ClusterMerger_Coll_D, BoxCenter[1]+0.5*ClusterMerger_Coll_B, BoxCenter[2] };

         for (long p=0; p<NPar1; p++)
         for (int d=1; d<=3; d++)
            ParData_AllRank[p][d] += ClusterCenter1[d-1];

         for (long p=NPar1; p<NParAll; p++)
         for (int d=1; d<=3; d++)
            ParData_AllRank[p][d] += ClusterCenter2[d-1];
      }

      else
      {
         for (long p=0; p<NParAll; p++)
         for (int d=1; d<=3; d++)
            ParData_AllRank[p][d] += BoxCenter[d-1];
      }


//    add the bulk velocity (to velocity-x only)
      if ( ClusterMerger_Coll )
      {
         for (long p=0;     p<NPar1;   p++)  ParData_AllRank[p][4] += ClusterMerger_Coll_BulkVel1;
         for (long p=NPar1; p<NParAll; p++)  ParData_AllRank[p][4] += ClusterMerger_Coll_BulkVel2;
      }

      Aux_Message( stdout, "done\n" );
   } // if ( MPI_Rank == 0 )

   MPI_Barrier( MPI_COMM_WORLD );


// get the number of particles in each rank and set the corresponding offsets
   if ( NParVar*amr->Par->NPar_Active_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "Total number of particle attributes to be sent (%ld) exceeds the maximum integer (%ld) !!\n",
                 (long)NParVar*amr->Par->NPar_Active_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank], NPar_MyRank=(int)amr->Par->NPar_AcPlusInac;

   MPI_Gather( &NPar_MyRank, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)  NSend[r] *= NParVar;

      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particles from the root rank to all ranks
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Send particles to all ranks ... " );
   double (*ParData_MyRank)[NParVar] = new double [NPar_MyRank][NParVar];

   MPI_Scatterv( ParData_AllRank[0], NSend, SendDisp, MPI_DOUBLE, ParData_MyRank[0], NPar_MyRank*NParVar, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   delete [] ParData_AllRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Store particles to the repository ... " );

   for (long p=0; p<NPar_MyRank; p++)
   {
      amr->Par->Mass[p] = ParData_MyRank[p][0];
      amr->Par->PosX[p] = ParData_MyRank[p][1];
      amr->Par->PosY[p] = ParData_MyRank[p][2];
      amr->Par->PosZ[p] = ParData_MyRank[p][3];
      amr->Par->VelX[p] = ParData_MyRank[p][4];
      amr->Par->VelY[p] = ParData_MyRank[p][5];
      amr->Par->VelZ[p] = ParData_MyRank[p][6];

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_MyRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
