#include "GAMER.h"

#ifdef PARTICLE

// floating-point type in the input particle file
typedef double real_par_in;
//typedef float  real_par_in;

extern char    Merger_File_Par1[1000];
extern char    Merger_File_Par2[1000];
extern bool    Merger_Coll;
extern double  Merger_Coll_PosX1;
extern double  Merger_Coll_PosY1;
extern double  Merger_Coll_PosX2;
extern double  Merger_Coll_PosY2;
extern double  Merger_Coll_VelX1;
extern double  Merger_Coll_VelY1;
extern double  Merger_Coll_VelX2;
extern double  Merger_Coll_VelY2;

long Read_Particle_Number(std::string filename);
void Read_Particles(std::string filename, long offset, long num,
                    real_par_in mass[], real_par_in xpos[], real_par_in ypos[], 
                    real_par_in zpos[], real_par_in xvel[], real_par_in yvel[], 
                    real_par_in zvel[]);

//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_ClusterMerger
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
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
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
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_ClusterMerger( const long NPar_ThisRank, const long NPar_AllRank,
                                        real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                        real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                        real *AllAttribute[PAR_NATT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   const long NParAtt = 7;   // mass, pos*3, vel*3
// check file existence
   if ( !Aux_CheckFileExist(Merger_File_Par1) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par1 );

   if ( Merger_Coll  &&  !Aux_CheckFileExist(Merger_File_Par2) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par2 );

   const string filename1(Merger_File_Par1);
   const string filename2(Merger_File_Par2);

// check file size
   long NPar_EachCluster[2]={0,0}, NPar_AllCluster;

   if ( MPI_Rank == 0 ) 

      NPar_EachCluster[0] = Read_Particle_Number(filename1);

      Aux_Message( stdout, "   Number of particles in cluster 1 = %ld\n", 
                  NPar_EachCluster[0] );

      if ( Merger_Coll ) {
         NPar_EachCluster[1] = Read_Particle_Number(filename2);
         Aux_Message( stdout, "   Number of particles in cluster 2 = %ld\n", NPar_EachCluster[1] );
      }

   }

#ifndef SERIAL
   MPI_Bcast(NPar_EachCluster, 2, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

   NPar_AllCluster = NPar_EachCluster[0] + NPar_EachCluster[1];

   if ( NPar_AllCluster != NPar_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found in cluster [%ld] != expect [%ld] !!\n",
                 NPar_AllCluster, NPar_AllRank );

// prepare to load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Preparing to load data ... " );

   const int NCluster = ( Merger_Coll ) ? 2 : 1;
   long NPar_ThisRank_EachCluster[2]={0,0}, Offset[2];   // [0/1] --> cluster 1/2

   for (int c=0; c<NCluster; c++)
   {
//    get the number of particles loaded by each rank for each cluster
      long NPar_ThisCluster_EachRank[MPI_NRank];

      if ( c == 0 )  
         NPar_ThisRank_EachCluster[0] = NPar_EachCluster[0] / MPI_NRank + ( (MPI_Rank<NPar_EachCluster[0]%MPI_NRank)?1:0 );
      else           
         NPar_ThisRank_EachCluster[1] = NPar_ThisRank - NPar_ThisRank_EachCluster[0];

      MPI_Allgather( &NPar_ThisRank_EachCluster[c], 1, MPI_LONG, NPar_ThisCluster_EachRank, 1, MPI_LONG,                         MPI_COMM_WORLD );

//    check if the total number of particles is correct
      long NPar_Check = 0;
      for (int r=0; r<MPI_NRank; r++)  
         NPar_Check += NPar_ThisCluster_EachRank[r];
      if ( NPar_Check != NPar_EachCluster[c] )
         Aux_Error( ERROR_INFO, "total number of particles in cluster %d: found (%ld) != expect (%ld) !!\n",
                    c, NPar_Check, NPar_EachCluster[c] );

//    set the file offset for this rank
      Offset[c] = 0;
      for (int r=0; r<MPI_Rank; r++) 
         Offset[c] = Offset[c] + NPar_ThisCluster_EachRank[r];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// load data to the particle repository

   for (int c=0; c<NCluster; c++)
   {
//    load data
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading cluster %d ... ", c );

      real_par_in *mass = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *xpos = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *ypos = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *zpos = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *xvel = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *yvel = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *zvel = new real_par_in [NPar_ThisRank_EachCluster[c]];

      const string filename((c==0)? Merger_File_Par1:Merger_File_Par2);

      Read_Particles(filename, Offset[c], NPar_ThisRank_EachCluster[c],
                     mass, xpos, ypos, zpos, xvel, yvel, zvel);

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "done\n" );

//    store data to the particle repository
      if ( MPI_Rank == 0 )    
         Aux_Message( stdout, "   Storing cluster %d to the particle repository ... ", c );

      for (long p=0; p<NPar_ThisRank_EachCluster[c]; p++)
      {
//       particle index offset
         const long pp = p + c*NPar_ThisRank_EachCluster[0];

//       --> convert to code unit before storing to the particle repository to avoid floating-point overflow
//       --> we have assumed that the loaded data are in cgs
         ParMass[pp] = real( mass[p] / UNIT_M );

         ParPosX[pp] = real( xpos[p] / UNIT_L );
         ParPosY[pp] = real( ypos[p] / UNIT_L );
         ParPosZ[pp] = real( zpos[p] / UNIT_L );

         ParVelX[pp] = real( xvel[p] / UNIT_V );
         ParVelY[pp] = real( yvel[p] / UNIT_V );
         ParVelZ[pp] = real( zvel[p] / UNIT_V );

//       synchronize all particles to the physical time at the base level
         ParTime[pp] = Time[0];
      }

      delete [] mass;
      delete [] xpos;
      delete [] ypos;
      delete [] zpos;
      delete [] xvel;
      delete [] yvel;
      delete [] zvel;

   } // for (int c=0; c<NCluster; c++)

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

// shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center and adding bulk velocity ... " );

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };

   if ( Merger_Coll )
   {
      const double ClusterCenter1[3]
         = { Merger_Coll_PosX1, Merger_Coll_PosY1, BoxCenter[2] };
      const double ClusterCenter2[3]
         = { Merger_Coll_PosX2, Merger_Coll_PosY2, BoxCenter[2] };

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


// add the bulk velocity
   if ( Merger_Coll )
   {
      for (long p=0; p<NPar_ThisRank_EachCluster[0]; p++) {  
         ParVelX[p] += Merger_Coll_VelX1;
         ParVelY[p] += Merger_Coll_VelY1;
      }
      for (long p=NPar_ThisRank_EachCluster[0]; p<NPar_ThisRank; p++) {  
         ParVelX[p] += Merger_Coll_VelX2;
         ParVelY[p] += Merger_Coll_VelY2;
      }
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_ClusterMerger



#endif // #ifdef PARTICLE





