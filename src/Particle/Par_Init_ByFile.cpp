#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFile
// Description :  Initialize particle attributes from a file
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box even for the periodic BC
//                3. Particles lying outside the active region will be removed by Par_Aux_InitCheck()
//                   if non-periodic B.C. is adopted
//                4. Particles loaded here are only temporarily stored in this rank
//                   --> They will be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                5. Currently the target file name is fixed to "PAR_IC"
//                6. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
//                7. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by "Init_Parallelization"
//
// Parameter   :  None
//
// Return      :  amr->Par->Time,Mass,PosX/Y/Z,VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFile()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char FileName[] = "PAR_IC";
   const int  NParVar    = 7;             // mass, pos*3, vel*3


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist for PAR_INIT == PAR_INIT_BY_FILE !!\n", FileName );

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = long(NParVar)*amr->Par->NPar_Active_AllRank*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 FileName, FileSize, ExpectSize );

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

   for (int r=0; r<MPI_Rank; r++)   FileOffset = FileOffset + long(NParVar)*NPar_EachRank[r]*sizeof(real);


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ... " );

   real (*ParData_ThisRank)[NParVar] = new real [amr->Par->NPar_AcPlusInac][NParVar];

// note that fread may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( FileName, "rb" );

   fseek( File, FileOffset, SEEK_SET );
   fread( ParData_ThisRank, sizeof(real), long(NParVar)*amr->Par->NPar_AcPlusInac, File );

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
   {
//    we have assumed that particle data in the file are stored in the order [ID][Mass/PosX/Y/Z/VelX/Y/Z]
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

} // FUNCTION : Par_Init_ByFile



#endif // #ifdef PARTICLE
