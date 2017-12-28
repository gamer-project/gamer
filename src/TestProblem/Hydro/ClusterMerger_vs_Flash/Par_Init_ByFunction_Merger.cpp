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
// Description :  Initialize all particle attributes for the merging cluster test
//                --> Modified from "Par_Init_ByFile.cpp"
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//                4. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
//
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParPassive    : Particle passive attributes pointer array with the size [PAR_NPASSIVE][NPar_ThisRank]
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParPassive
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Merger( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                 real *ParPassive[PAR_NPASSIVE] )
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
   long  NPar_EachCluster[2]={0,0}, NPar_AllCluster;

   File  = fopen( Merger_File_Par1, "rb" );
   fseek( File, 0, SEEK_END );
   NPar_EachCluster[0] = ftell( File ) / ( NParVar*sizeof(real_par_in) );
   fclose( File );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Number of particles in cluster 1 = %ld\n", NPar_EachCluster[0] );

   if ( Merger_Coll ) {
   File  = fopen( Merger_File_Par2, "rb" );
   fseek( File, 0, SEEK_END );
   NPar_EachCluster[1] = ftell( File ) / ( NParVar*sizeof(real_par_in) );
   fclose( File );
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Number of particles in cluster 2 = %ld\n", NPar_EachCluster[1] );
   }

   NPar_AllCluster = NPar_EachCluster[0] + NPar_EachCluster[1];

   if ( NPar_AllCluster != NPar_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                 NPar_AllCluster, NPar_AllRank );

   MPI_Barrier( MPI_COMM_WORLD );


// prepare to load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Preparing to load data ... " );

   const int NCluster = ( Merger_Coll ) ? 2 : 1;
   long NPar_ThisRank_EachCluster[2]={0,0}, FileOffset[2];   // [0/1] --> cluster 1/2

   for (int c=0; c<NCluster; c++)
   {
//    get the number of particles loaded by each rank for each cluster
      long NPar_ThisCluster_EachRank[MPI_NRank];

      if ( c == 0 )  NPar_ThisRank_EachCluster[0] = NPar_EachCluster[0] / MPI_NRank + ( (MPI_Rank<NPar_EachCluster[0]%MPI_NRank)?1:0 );
      else           NPar_ThisRank_EachCluster[1] = NPar_ThisRank - NPar_ThisRank_EachCluster[0];

      MPI_Allgather( &NPar_ThisRank_EachCluster[c], 1, MPI_LONG, NPar_ThisCluster_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

//    check if the total number of particles is correct
      long NPar_Check = 0;
      for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_ThisCluster_EachRank[r];
      if ( NPar_Check != NPar_EachCluster[c] )
         Aux_Error( ERROR_INFO, "total number of particles in cluster %d: found (%ld) != expect (%ld) !!\n",
                    c, NPar_Check, NPar_EachCluster[c] );

//    set the file offset for this rank
      FileOffset[c] = 0;
      for (int r=0; r<MPI_Rank; r++)   FileOffset[c] = FileOffset[c] + NPar_ThisCluster_EachRank[r]*NParVar*sizeof(real_par_in);
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// load data to the particle repository
   real_par_in (*ParData_ThisRank)[NParVar] = new real_par_in [ MAX(NPar_ThisRank_EachCluster[0], NPar_ThisRank_EachCluster[1]) ][NParVar];

   for (int c=0; c<NCluster; c++)
   {
//    load data
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading cluster %d ... ", c );

//    note that fread may fail for large files if sizeof(size_t) == 4 instead of 8
      FILE *File = fopen( (c==0)?Merger_File_Par1:Merger_File_Par2, "rb" );

      fseek( File, FileOffset[c], SEEK_SET );
      fread( ParData_ThisRank, sizeof(real_par_in), NPar_ThisRank_EachCluster[c]*NParVar, File );
      fclose( File );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


//    store data to the particle repository
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing cluster %d to the particle repository ... ", c );

      for (long p=0; p<NPar_ThisRank_EachCluster[c]; p++)
      {
//       particle index offset
         const long pp = p + c*NPar_ThisRank_EachCluster[0];

//       we have assumed that particle data in the file are stored in the order [ID][Mass/PosX/Y/Z/VelX/Y/Z]
//       --> convert to code unit before storing to the particle repository to avoid floating-point overflow
//       --> we have assumed that the loaded data are in cgs
         ParMass[pp] = real( ParData_ThisRank[p][0] / UNIT_M );

         ParPosX[pp] = real( ParData_ThisRank[p][1] / UNIT_L );
         ParPosY[pp] = real( ParData_ThisRank[p][2] / UNIT_L );
         ParPosZ[pp] = real( ParData_ThisRank[p][3] / UNIT_L );

         ParVelX[pp] = real( ParData_ThisRank[p][4] / UNIT_V );
         ParVelY[pp] = real( ParData_ThisRank[p][5] / UNIT_V );
         ParVelZ[pp] = real( ParData_ThisRank[p][6] / UNIT_V );

//       synchronize all particles to the physical time at the base level
         ParTime[pp] = Time[0];
      }

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int c=0; c<NCluster; c++)

   delete [] ParData_ThisRank;


// shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center and adding bulk velocity ... " );

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };

   if ( Merger_Coll )
   {
      const double ClusterCenter1[3]
         = { BoxCenter[0]-0.5*Merger_Coll_D, BoxCenter[1]-0.5*Merger_Coll_B, BoxCenter[2] };
      const double ClusterCenter2[3]
         = { BoxCenter[0]+0.5*Merger_Coll_D, BoxCenter[1]+0.5*Merger_Coll_B, BoxCenter[2] };

      for (long p=0; p<NPar_ThisRank_EachCluster[0]; p++)
      for (int d=0; d<3; d++)
         ParPos[d][p] += ClusterCenter1[d];

      for (long p=NPar_ThisRank_EachCluster[0]; p<NPar_ThisRank; p++)
      for (int d=0; d<3; d++)
         ParPos[d][p] += ClusterCenter2[d];
   }

   else
   {
      for (long p=0; p<NPar_ThisRank; p++)
      for (int d=0; d<3; d++)
         ParPos[d][p] += BoxCenter[d];
   }


// add the bulk velocity (to velocity-x only)
   if ( Merger_Coll )
   {
      for (long p=0;                            p<NPar_ThisRank_EachCluster[0]; p++)   ParVelX[p] += Merger_Coll_BulkVel1;
      for (long p=NPar_ThisRank_EachCluster[0]; p<NPar_ThisRank;                p++)   ParVelX[p] += Merger_Coll_BulkVel2;
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Merger



#endif // #ifdef PARTICLE





