#include "GAMER.h"

#ifdef PARTICLE

// floating-point type in the input particle file
typedef double real_par_in;
//typedef float  real_par_in;

extern char    Merger_File_Par1[1000];
extern char    Merger_File_Par2[1000];
extern bool    Merger_Coll;
extern double  Merger_Coll_D;
extern double  Merger_Coll_B;
extern double  Merger_Coll_BulkVel1;
extern double  Merger_Coll_BulkVel2;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Merger
// Description :  Initialize particle attributes for the merging cluster test
//                --> Modified from "Par_Init_ByFile.cpp"
//
// Note        :  1. Invoked by "Init_GAMER" using the function pointer "Par_Init_ByFunction_Ptr"
//                2. Particles lying outside the active region will be removed by "Par_Aux_InitCheck"
//                   if non-periodic B.C. is adopted
//                3. Particles set here are only temporarily stored in this rank
//                   --> They will be redistributed when calling "Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()"
//                4. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by "Init_Parallelization()" and stored in "amr->Par->NPar_AcPlusInac"
//                5. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
//
// Parameter   :  None
//
// Return      :  amr->Par->Time,Mass,PosX/Y/Z,VelX/Y/Z
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Merger()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const long NParVar = 7;   // mass, pos*3, vel*3

// check file existence
   if ( !Aux_CheckFileExist(Merger_File_Par1) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par1 );

   if ( Merger_Coll  &&  !Aux_CheckFileExist(Merger_File_Par2) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par2 );


// check file size
   FILE *File=NULL;
   long  NPar[2]={0,0}, NParAll;

   File  = fopen( Merger_File_Par1, "rb" );
   fseek( File, 0, SEEK_END );
   NPar[0] = ftell( File ) / ( NParVar*sizeof(real_par_in) );
   fclose( File );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Number of particles in cluster 1 = %ld\n", NPar[0] );

   if ( Merger_Coll ) {
   File  = fopen( Merger_File_Par2, "rb" );
   fseek( File, 0, SEEK_END );
   NPar[1] = ftell( File ) / ( NParVar*sizeof(real_par_in) );
   fclose( File );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Number of particles in cluster 2 = %ld\n", NPar[1] );
   }

   NParAll = NPar[0] + NPar[1];

   if ( NParAll != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                 NParAll, amr->Par->NPar_Active_AllRank );

   MPI_Barrier( MPI_COMM_WORLD );


// prepare to load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Preparing to load data ... " );

   const int NCluster = ( Merger_Coll ) ? 2 : 1;
   long NPar_ThisRank[2]={0,0}, FileOffset[2];  // [0/1] --> cluster 1/2

   for (int c=0; c<NCluster; c++)
   {
//    get the number of particles loaded by each rank for each cluster
      long NPar_EachRank[MPI_NRank];

      if ( c == 0 )  NPar_ThisRank[0] = NPar[0] / MPI_NRank + ( (MPI_Rank<NPar[0]%MPI_NRank)?1:0 );
      else           NPar_ThisRank[1] = amr->Par->NPar_AcPlusInac - NPar_ThisRank[0];

      MPI_Allgather( &NPar_ThisRank[c], 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

//    check if the total number of particles is correct
      long NPar_Check = 0;
      for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_EachRank[r];
      if ( NPar_Check != NPar[c] )
         Aux_Error( ERROR_INFO, "total number of particles in cluster %d: found (%ld) != expect (%ld) !!\n",
                    c, NPar_Check, NPar[c] );

//    set the file offset for this rank
      FileOffset[c] = 0;
      for (int r=0; r<MPI_Rank; r++)   FileOffset[c] = FileOffset[c] + NPar_EachRank[r]*NParVar*sizeof(real_par_in);
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// load data to the particle repository
   real_par_in (*ParData_ThisRank)[NParVar] = new real_par_in [ MAX(NPar_ThisRank[0], NPar_ThisRank[1]) ][NParVar];

   for (int c=0; c<NCluster; c++)
   {
//    load data
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading cluster %d ... ", c );

//    note that fread may fail for large files if sizeof(size_t) == 4 instead of 8
      FILE *File = fopen( (c==0)?Merger_File_Par1:Merger_File_Par2, "rb" );

      fseek( File, FileOffset[c], SEEK_SET );
      fread( ParData_ThisRank, sizeof(real_par_in), NPar_ThisRank[c]*NParVar, File );
      fclose( File );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


//    store data to the particle repository
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing cluster %d to the particle repository ... ", c );

      for (long p=0; p<NPar_ThisRank[c]; p++)
      {
//       particle index offset
         const long pp = p + c*NPar_ThisRank[0];

//       we have assumed that particle data in the file are stored in the order [ID][Mass/PosX/Y/Z/VelX/Y/Z]
//       --> convert to code unit before storing to the particle repository to avoid floating-point overflow
//       --> we have assumed that the loaded data are in cgs
         amr->Par->Mass[pp] = real( ParData_ThisRank[p][0] / UNIT_M );

         amr->Par->PosX[pp] = real( ParData_ThisRank[p][1] / UNIT_L );
         amr->Par->PosY[pp] = real( ParData_ThisRank[p][2] / UNIT_L );
         amr->Par->PosZ[pp] = real( ParData_ThisRank[p][3] / UNIT_L );

         amr->Par->VelX[pp] = real( ParData_ThisRank[p][4] / UNIT_V );
         amr->Par->VelY[pp] = real( ParData_ThisRank[p][5] / UNIT_V );
         amr->Par->VelZ[pp] = real( ParData_ThisRank[p][6] / UNIT_V );

//       synchronize all particles to the physical time at the base level
         amr->Par->Time[pp] = Time[0];
      }

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int c=0; c<NCluster; c++)

   delete [] ParData_ThisRank;


// shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center and adding bulk velocity ... " );

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   real *Pos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };

   if ( Merger_Coll )
   {
      const double ClusterCenter1[3]
         = { BoxCenter[0]-0.5*Merger_Coll_D, BoxCenter[1]-0.5*Merger_Coll_B, BoxCenter[2] };
      const double ClusterCenter2[3]
         = { BoxCenter[0]+0.5*Merger_Coll_D, BoxCenter[1]+0.5*Merger_Coll_B, BoxCenter[2] };

      for (long p=0; p<NPar_ThisRank[0]; p++)
      for (int d=0; d<3; d++)
         Pos[d][p] += ClusterCenter1[d];

      for (long p=NPar_ThisRank[0]; p<amr->Par->NPar_AcPlusInac; p++)
      for (int d=0; d<3; d++)
         Pos[d][p] += ClusterCenter2[d];
   }

   else
   {
      for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
      for (int d=0; d<3; d++)
         Pos[d][p] += BoxCenter[d];
   }


// add the bulk velocity (to velocity-x only)
   if ( Merger_Coll )
   {
      for (long p=0;                p<NPar_ThisRank[0];          p++)  amr->Par->VelX[p] += Merger_Coll_BulkVel1;
      for (long p=NPar_ThisRank[0]; p<amr->Par->NPar_AcPlusInac; p++)  amr->Par->VelX[p] += Merger_Coll_BulkVel2;
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Merger



#endif // #ifdef PARTICLE





