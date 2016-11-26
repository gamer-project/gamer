#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE

extern bool ClusterMerger_Coll;
extern char ClusterMerger_File_Par[1000];




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction
// Description :  Initialize the particle position and velocity
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


// check
   if ( !Aux_CheckFileExist(ClusterMerger_File_Par) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", ClusterMerger_File_Par );

   FILE *FileTemp = fopen( ClusterMerger_File_Par, "rb" );

   fseek( FileTemp, 0, SEEK_END );

// assuming that particle data are stored in double precision
   const int NParVar     = 7; // mass, pos*3, vel*3
   const long ExpectSize = long(NParVar)*amr->Par->NPar_Active_AllRank*sizeof(double);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 ClusterMerger_File_Par, FileSize, ExpectSize );

   fclose( FileTemp );

   MPI_Barrier( MPI_COMM_WORLD );


// set the file offset for this rank
   long NPar_EachRank[MPI_NRank], NPar_Check=0, FileOffset=0;

   MPI_Allgather( &amr->Par->NPar_AcPlusInac, 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

// check if the total number of particles is correct
   for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_EachRank[r];
   if ( NPar_Check != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found (%ld) != expect (%ld) !!\n",
                 NPar_Check, amr->Par->NPar_Active_AllRank );

   for (int r=0; r<MPI_Rank; r++)   FileOffset += long(NParVar)*NPar_EachRank[r]*sizeof(double);


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading particles ... " );

   double (*ParData_ThisRank)[NParVar] = new double [amr->Par->NPar_AcPlusInac][NParVar];

// note that fread may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( ClusterMerger_File_Par, "rb" );

   fseek( File, FileOffset, SEEK_SET );
   fread( ParData_ThisRank, sizeof(double), long(NParVar)*amr->Par->NPar_AcPlusInac, File );

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// convert to code units (assuming loaded data are in cgs)
// --> assuming particle data are stored in the order [ID][Mass/PosX/Y/Z/VelX/Y/Z]
   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      ParData_ThisRank[p][0] /= UNIT_M;

      ParData_ThisRank[p][1] /= UNIT_L;
      ParData_ThisRank[p][2] /= UNIT_L;
      ParData_ThisRank[p][3] /= UNIT_L;

      ParData_ThisRank[p][4] /= UNIT_V;
      ParData_ThisRank[p][5] /= UNIT_V;
      ParData_ThisRank[p][6] /= UNIT_V;
   }


// shift center (assuming the center of loaded particles = [0,0,0])
   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   if ( ClusterMerger_Coll )
   {
      Aux_Error( ERROR_INFO, "NOT SUPPORETD YET !!\n" );
   }

   else
   {
      for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
      for (int d=1; d<=3; d++)
         ParData_ThisRank[p][d] += BoxCenter[d-1];
   }


// add the relative bulk velocity
   if ( ClusterMerger_Coll )
   {
      Aux_Error( ERROR_INFO, "NOT SUPPORETD YET !!\n" );
   }


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
      amr->Par->Mass[p] = ParData_ThisRank[p][0];
      amr->Par->PosX[p] = ParData_ThisRank[p][1];
      amr->Par->PosY[p] = ParData_ThisRank[p][2];
      amr->Par->PosZ[p] = ParData_ThisRank[p][3];
      amr->Par->VelX[p] = ParData_ThisRank[p][4];
      amr->Par->VelY[p] = ParData_ThisRank[p][5];
      amr->Par->VelZ[p] = ParData_ThisRank[p][6];

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_ThisRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
