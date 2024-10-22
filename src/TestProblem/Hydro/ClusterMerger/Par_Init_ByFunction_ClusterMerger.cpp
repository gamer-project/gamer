#include "GAMER.h"

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

#include <string>


// floating-point type in the input particle file
typedef double real_par_in;
//typedef float  real_par_in;

extern char    Merger_File_Par1[1000];
extern char    Merger_File_Par2[1000];
extern char    Merger_File_Par3[1000];
extern int     Merger_Coll_NumHalos;
extern double  Merger_Coll_PosX1;
extern double  Merger_Coll_PosY1;
extern double  Merger_Coll_PosX2;
extern double  Merger_Coll_PosY2;
extern double  Merger_Coll_PosX3;
extern double  Merger_Coll_PosY3;
extern double  Merger_Coll_VelX1;
extern double  Merger_Coll_VelY1;
extern double  Merger_Coll_VelX2;
extern double  Merger_Coll_VelY2;
extern double  Merger_Coll_VelX3;
extern double  Merger_Coll_VelY3;
extern bool    Merger_Coll_LabelCenter;
extern long    NPar_EachCluster[3];
extern long    NPar_AllCluster;


// variables that need to be record in Record__Center
// =======================================================================================
extern double  Bondi_MassBH1;
extern double  Bondi_MassBH2;
extern double  Bondi_MassBH3;
extern double  Mdot_BH1;
extern double  Mdot_BH2;
extern double  Mdot_BH3;

extern double CM_Bondi_SinkMass[3];
extern double CM_Bondi_SinkMomX[3];
extern double CM_Bondi_SinkMomY[3];
extern double CM_Bondi_SinkMomZ[3];
extern double CM_Bondi_SinkMomXAbs[3];
extern double CM_Bondi_SinkMomYAbs[3];
extern double CM_Bondi_SinkMomZAbs[3];
extern double CM_Bondi_SinkE[3];
extern double CM_Bondi_SinkEk[3];
extern double CM_Bondi_SinkEt[3];
extern int    CM_Bondi_SinkNCell[3];

extern double GasVel[3][3];  // gas velocity
extern double SoundSpeed[3];
extern double GasDens[3];
extern double RelativeVel[3];
extern double ColdGasMass[3];
extern double Jet_Vec[3][3]; // jet direction  
extern double Mdot[3]; // the feedback injection rate
extern double Pdot[3];
extern double Edot[3];
extern double E_inj_exp[3];
extern double M_inj_exp[3];
extern double dt_base;
double E_power_inj[3];  // the injection power
extern double ClusterCen[3][3];
extern double BH_Vel[3][3];
extern double R_acc;
extern bool   fixBH;
int num_par_sum[3] = {0, 0, 0};   // total number of particles inside the target region of each cluster
// =======================================================================================

#ifdef MASSIVE_PARTICLES
void Read_Particles_ClusterMerger(std::string filename, long offset, long num,
                                  real_par_in xpos[], real_par_in ypos[],
                                  real_par_in zpos[], real_par_in xvel[],
                                  real_par_in yvel[], real_par_in zvel[],
                                  real_par_in mass[], real_par_in ptype[]);
void GetClusterCenter( int lv, bool AdjustPos, bool AdjustVel, double Cen_old[][3], double Cen_new[][3], double Cen_Vel[][3] );

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
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------

void Par_Init_ByFunction_ClusterMerger( const long NPar_ThisRank, const long NPar_AllRank,
                                        real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                        real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                        real_par *ParType, real_par *AllAttribute[PAR_NATT_TOTAL] )
{

#ifdef SUPPORT_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   // prepare to load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Preparing to load data ... " );

   const int NCluster = Merger_Coll_NumHalos;
   long NPar_ThisRank_EachCluster[3]={0,0,0}, Offset[3];   // [0/1/2] --> cluster 1/2/3

   for (int c=0; c<NCluster; c++)
   {
      // get the number of particles loaded by each rank for each cluster
      long NPar_ThisCluster_EachRank[MPI_NRank];

      switch (c) {
      case 0:
         NPar_ThisRank_EachCluster[0] = NPar_EachCluster[0] / MPI_NRank + ( (MPI_Rank<NPar_EachCluster[0]%MPI_NRank)?1:0 );
         break;
      case 1:
         if (NCluster == 2)
            NPar_ThisRank_EachCluster[1] = NPar_ThisRank - NPar_ThisRank_EachCluster[0];
         else
            NPar_ThisRank_EachCluster[1] = NPar_EachCluster[1] / MPI_NRank + ( (MPI_Rank<NPar_EachCluster[1]%MPI_NRank)?1:0 );
         break;
      case 2:
         NPar_ThisRank_EachCluster[2] = NPar_ThisRank - NPar_ThisRank_EachCluster[0] - NPar_ThisRank_EachCluster[1];
         break;
      }

      MPI_Allgather( &NPar_ThisRank_EachCluster[c], 1, MPI_LONG, NPar_ThisCluster_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

      // check if the total number of particles is correct
      long NPar_Check = 0;
      for (int r=0; r<MPI_NRank; r++)
         NPar_Check += NPar_ThisCluster_EachRank[r];
      if ( NPar_Check != NPar_EachCluster[c] )
         Aux_Error( ERROR_INFO, "total number of particles in cluster %d: found (%ld) != expect (%ld) !!\n",
                    c, NPar_Check, NPar_EachCluster[c] );

      // set the file offset for this rank
      Offset[c] = 0;
      for (int r=0; r<MPI_Rank; r++)
         Offset[c] = Offset[c] + NPar_ThisCluster_EachRank[r];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

   // load data to the particle repository

   const std::string filenames[3] = { Merger_File_Par1, Merger_File_Par2, Merger_File_Par3 };

   for ( int c=0; c<NCluster; c++ )
   {
      // load data
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading cluster %d ... \n", c+1 );

      real_par_in *mass  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *xpos  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *ypos  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *zpos  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *xvel  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *yvel  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *zvel  = new real_par_in [NPar_ThisRank_EachCluster[c]];
      real_par_in *ptype = new real_par_in [NPar_ThisRank_EachCluster[c]];

      Read_Particles_ClusterMerger( filenames[c], Offset[c], NPar_ThisRank_EachCluster[c],
                                    xpos, ypos, zpos, xvel, yvel, zvel, mass, ptype );

#     ifndef TRACER
      for (long p=0; p<NPar_ThisRank_EachCluster[c]; p++) {
         if ( ptype[p] == PTYPE_TRACER )
            Aux_Error( ERROR_INFO,
"Tracer particles were found in the input data for cluster %d, but TRACER is not defined!\n",
                       c );
      }
#     endif

//    Reset the particle type for different clusters
      for (long p=0; p<NPar_ThisRank_EachCluster[c]; p++)   ptype[p] = PTYPE_CLUSTER + c;

      if ( MPI_Rank == 0 ) Aux_Message( stdout, "done\n" );

      // store data to the particle repository
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "   Storing cluster %d to the particle repository ... \n", c+1 );

      // Compute offsets for assigning particles

      double coffset;
      switch (c) {
      case 0:
         coffset = 0;
         break;
      case 1:
         coffset = NPar_ThisRank_EachCluster[0];
         break;
      case 2:
         coffset = NPar_ThisRank_EachCluster[0]+NPar_ThisRank_EachCluster[1];
         break;
      }

      for (long p=0; p<NPar_ThisRank_EachCluster[c]; p++)
      {
         // particle index offset
         const long pp = p + coffset;

         // set the particle type
         ParType[pp] = real_par( ptype[p] );

         // --> convert to code unit before storing to the particle repository to avoid floating-point overflow
         // --> we have assumed that the loaded data are in cgs

         ParPosX[pp] = real_par( xpos[p] / UNIT_L );
         ParPosY[pp] = real_par( ypos[p] / UNIT_L );
         ParPosZ[pp] = real_par( zpos[p] / UNIT_L );

         if ( ptype[p] == PTYPE_TRACER ) {
            // tracer particles have zero mass
            // and their velocities will be set by
            // the grid later
            ParMass[pp] = (real_par)0.0;
            ParVelX[pp] = (real_par)0.0;
            ParVelY[pp] = (real_par)0.0;
            ParVelX[pp] = (real_par)0.0;
         } else {
            // For massive particles get their mass
            // and velocity
            ParMass[pp] = real_par( mass[p] / UNIT_M );
            ParVelX[pp] = real_par( xvel[p] / UNIT_V );
            ParVelY[pp] = real_par( yvel[p] / UNIT_V );
            ParVelZ[pp] = real_par( zvel[p] / UNIT_V );
         }

         // synchronize all particles to the physical time at the base level
         ParTime[pp] = (real_par)Time[0];

      }

      delete [] mass;
      delete [] xpos;
      delete [] ypos;
      delete [] zpos;
      delete [] xvel;
      delete [] yvel;
      delete [] zvel;
      delete [] ptype;

   } // for (int c=0; c<NCluster; c++)

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

   // shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center and adding bulk velocity ... " );

   real_par *ParPos[3] = { ParPosX, ParPosY, ParPosZ };

   const double ClusterCenter1[3]
      = { Merger_Coll_PosX1, Merger_Coll_PosY1, amr->BoxCenter[2] };
   const double ClusterCenter2[3]
      = { Merger_Coll_PosX2, Merger_Coll_PosY2, amr->BoxCenter[2] };
   const double ClusterCenter3[3]
      = { Merger_Coll_PosX3, Merger_Coll_PosY3, amr->BoxCenter[2] };

   for (long p=0; p<NPar_ThisRank_EachCluster[0]; p++) {
      if ( (int)ParType[p] != PTYPE_TRACER ) {
         ParVelX[p] += Merger_Coll_VelX1;
         ParVelY[p] += Merger_Coll_VelY1;
      }
      for (int d=0; d<3; d++)
         ParPos[d][p] += ClusterCenter1[d];
   }

   for (long p=NPar_ThisRank_EachCluster[0]; p<NPar_ThisRank_EachCluster[0]+NPar_ThisRank_EachCluster[1]; p++) {
      if ( (int)ParType[p] != PTYPE_TRACER ) {
         ParVelX[p] += Merger_Coll_VelX2;
         ParVelY[p] += Merger_Coll_VelY2;
      }
      for (int d=0; d<3; d++)
         ParPos[d][p] += ClusterCenter2[d];
   }

   for (long p=NPar_ThisRank_EachCluster[0]+NPar_ThisRank_EachCluster[1]; p<NPar_ThisRank; p++) {
      if ( (int)ParType[p] != PTYPE_TRACER ) {
         ParVelX[p] += Merger_Coll_VelX3;
         ParVelY[p] += Merger_Coll_VelY3;
      }
      for (int d=0; d<3; d++)
         ParPos[d][p] += ClusterCenter3[d];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   // label cluster centers
   if ( Merger_Coll_LabelCenter ) {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Labeling cluster centers ... " );

      const double Centers[3][3] = {  { ClusterCenter1[0], ClusterCenter1[1], ClusterCenter1[2] },
                                      { ClusterCenter2[0], ClusterCenter2[1], ClusterCenter2[2] },
                                      { ClusterCenter3[0], ClusterCenter3[1], ClusterCenter3[2] }  };
      long pidx_offset = 0;

      for (int c=0; c<NCluster; c++) {
         long   min_pidx   = -1;
         real   min_pos[3] = { NULL_REAL, NULL_REAL, NULL_REAL };
         double min_r      = __DBL_MAX__;

         // get the particle in this rank closest to the cluster center
         for (long p=pidx_offset; p<pidx_offset+NPar_ThisRank_EachCluster[c]; p++) {
            const double r = SQR( ParPos[0][p] - Centers[c][0] ) +
                             SQR( ParPos[1][p] - Centers[c][1] ) +
                             SQR( ParPos[2][p] - Centers[c][2] );
            if ( r < min_r ) {
               min_pidx   = p;
               min_r      = r;
               min_pos[0] = ParPos[0][p];
               min_pos[1] = ParPos[1][p];
               min_pos[2] = ParPos[2][p];
            }
         }

         // collect data among all ranks
         double min_r_allrank;
         int    NFound_ThisRank=0, NFound_AllRank;
         MPI_Allreduce( &min_r, &min_r_allrank, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
         if ( min_r == min_r_allrank ) {
            ParType[min_pidx] = PTYPE_CEN + c;
            NFound_ThisRank = 1;
         }

         // check if one and only one particle is labeled
         MPI_Allreduce( &NFound_ThisRank, &NFound_AllRank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
         if ( NFound_AllRank != 1 )
            Aux_Error( ERROR_INFO, "NFound_AllRank (%d) != 1 for cluster %d !!\n", NFound_AllRank, c );

         // update the particle index offset for the next cluster
         pidx_offset += NPar_ThisRank_EachCluster[c];
      } // for (int c=0; c<NCluster; c++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // if ( Merger_Coll_LabelCenter )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

#endif // #ifdef SUPPORT_HDF5

} // FUNCTION : Par_Init_ByFunction_ClusterMerger

#ifdef SUPPORT_HDF5

void Read_Particles_ClusterMerger( std::string filename, long offset, long num,
                                   real_par_in xpos[], real_par_in ypos[],
                                   real_par_in zpos[], real_par_in xvel[],
                                   real_par_in yvel[], real_par_in zvel[],
                                   real_par_in mass[], real_par_in ptype[] )
{

   hid_t   file_id, dataset, dataspace, memspace;
   herr_t  status;

   hsize_t start[2], stride[2], count[2], dims[2], maxdims[2];
   hsize_t start1d[1], stride1d[1], count1d[1], dims1d[1], maxdims1d[1];
   hsize_t start0[1];

   int rank;

   stride[0] = 1;
   stride[1] = 1;
   start[0] = (hsize_t)offset;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset   = H5Dopen(file_id, "particle_position", H5P_DEFAULT);

   dataspace = H5Dget_space(dataset);

   rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   count[0] = (hsize_t)num;
   count[1] = 1;

   dims[0] = count[0];
   dims[1] = 1;

   count1d[0] = (hsize_t)num;
   dims1d[0] = count1d[0];
   stride1d[0] = 1;
   start1d[0] = 0;
   start[1] = 0;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, xpos);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle x-position!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   dataspace = H5Dget_space(dataset);

   start[1] = 1;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, ypos);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle y-position!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   dataspace = H5Dget_space(dataset);

   start[1] = 2;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, zpos);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle z-position!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   dataset   = H5Dopen(file_id, "particle_velocity", H5P_DEFAULT);

   dataspace = H5Dget_space(dataset);

   start[1] = 0;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, xvel);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle x-velocity!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);

   dataspace = H5Dget_space(dataset);

   start[1] = 1;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, yvel);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle y-velocity!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);

   dataspace = H5Dget_space(dataset);

   start[1] = 2;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start,
                                 stride, count, NULL);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, zvel);

   if (status < 0) {
      Aux_Message(stderr, "Could not read particle z-velocity!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   dataset   = H5Dopen(file_id, "particle_mass", H5P_DEFAULT);

   dataspace = H5Dget_space(dataset);

   start1d[0] = (hsize_t)offset;
   start0[0] = 0;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   rank      = H5Sget_simple_extent_dims(dataspace, dims1d, maxdims1d);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start0,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, mass);

   if (status < 0) {
      Aux_Message( stderr, "Could not read particle mass!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   dataset   = H5Dopen(file_id, "particle_type", H5P_DEFAULT);

   dataspace = H5Dget_space(dataset);

   start1d[0] = (hsize_t)offset;
   start0[0] = 0;

   status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start1d,
                                 stride1d, count1d, NULL);
   rank      = H5Sget_simple_extent_dims(dataspace, dims1d, maxdims1d);
   memspace = H5Screate_simple(1, dims1d, NULL);
   status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start0,
                                 stride1d, count1d, NULL);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
                     H5P_DEFAULT, ptype);

   if (status < 0) {
      Aux_Message( stderr, "Could not read particle type!!\n");
   }

   H5Sclose(memspace);
   H5Sclose(dataspace);
   H5Dclose(dataset);

   H5Fclose(file_id);

   return;

} // FUNCTION : Read_Particles_ClusterMerger

#endif // #ifdef SUPPORT_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_ClusterMerger
// Description :  Record the cluster centers
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//                4. Must enable Merger_Coll_LabelCenter
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_ClusterMerger()
{

   const char FileName[] = "Record__Center";
   static bool FirstTime = true;

   // header
   if ( FirstTime )
   {
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%13s%14s",  "Time", "Step" );
         for (int c=0; c<Merger_Coll_NumHalos; c++){
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "x", c, "y", c, "z", c);
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "BHVel_x[km/s]", c, "BHVel_y", c, "BHVel_z", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "GasVel_x", c,"GasVel_y", c,"GasVel_z", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "RelativeVel", c, "SoundSpeed", c, "GasDens(cgs)", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "mass_BH[Msun]", c, "Mdot(cgs)", c, "NVoidCell", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "MomXInj(cgs)", c, "MomYInj", c,"MomZInj", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "MomXInjAbs", c, "MomYInjAbs", c, "MomZInjAbs", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "EInj_exp[erg]", c, "E_Inj[erg]", c, "E_Inj_err", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "Ek_Inj[erg]", c, "Et_Inj[erg]", c, "PowerInj(cgs)", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "MInjexp[Msun]", c, "MassInj[Msun]", c, "M_Inj_err", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "Mdot(cgs)", c, "Pdot(cgs)", c, "Edot(cgs)", c );
            fprintf( File_User, " %13s%1d %13s%1d %13s%1d", "Jet_Vec_x", c, "Jet_Vec_y", c, "Jet_Vec_z", c );
            fprintf( File_User, " %13s%1d %13s%1d", "num_par_sum", c, "ColdGasMass", c );
         }
         fprintf( File_User, "\n" );
         fclose( File_User );
      }

      FirstTime = false;
   } // if ( FirstTime )

   double Bondi_MassBH[3] = { Bondi_MassBH1, Bondi_MassBH2, Bondi_MassBH3 };
   double Mdot_BH[3] = { Mdot_BH1, Mdot_BH2, Mdot_BH3 };

// sum over the variables and convert units 
   int SinkNCell_Sum[3];
   double Mass_Sum[3], MomX_Sum[3], MomY_Sum[3], MomZ_Sum[3], MomXAbs_Sum[3], MomYAbs_Sum[3], MomZAbs_Sum[3], E_Sum[3], Ek_Sum[3], Et_Sum[3];

   for (int c=0; c<Merger_Coll_NumHalos; c++){
      MPI_Reduce( &CM_Bondi_SinkNCell[c], &SinkNCell_Sum[c], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMass[c],    &Mass_Sum[c],    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomX[c],    &MomX_Sum[c],    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomY[c],    &MomY_Sum[c],    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomZ[c],    &MomZ_Sum[c],    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomXAbs[c], &MomXAbs_Sum[c], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomYAbs[c], &MomYAbs_Sum[c], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkMomZAbs[c], &MomZAbs_Sum[c], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkE[c],       &E_Sum[c],       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkEk[c],      &Ek_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CM_Bondi_SinkEt[c],      &Et_Sum[c],      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   
      Mass_Sum[c]    *= UNIT_M/Const_Msun;
      MomX_Sum[c]    *= UNIT_M*UNIT_V;
      MomY_Sum[c]    *= UNIT_M*UNIT_V;
      MomZ_Sum[c]    *= UNIT_M*UNIT_V;
      MomXAbs_Sum[c] *= UNIT_M*UNIT_V;
      MomYAbs_Sum[c] *= UNIT_M*UNIT_V;
      MomZAbs_Sum[c] *= UNIT_M*UNIT_V;
      E_Sum[c]       *= UNIT_E;
      Ek_Sum[c]      *= UNIT_E;
      Et_Sum[c]      *= UNIT_E;
      E_inj_exp[c]   *= UNIT_E; 
      M_inj_exp[c]   *= UNIT_M/Const_Msun;
   }

   for (int c=0; c<Merger_Coll_NumHalos; c++)   E_power_inj[c] = E_Sum[c]/(dt_base*UNIT_T);

   // output the properties of the cluster centers 
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e%14ld", Time[0], Step );
      for (int c=0; c<Merger_Coll_NumHalos; c++){
         fprintf( File_User, " %14.7e %14.7e %14.7e", ClusterCen[c][0], ClusterCen[c][1], ClusterCen[c][2] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", BH_Vel[c][0]*UNIT_V/(Const_km/Const_s), BH_Vel[c][1]*UNIT_V/(Const_km/Const_s), BH_Vel[c][2]*UNIT_V/(Const_km/Const_s) );
         fprintf( File_User, " %14.7e %14.7e %14.7e", GasVel[c][0]*UNIT_V/(Const_km/Const_s), GasVel[c][1]*UNIT_V/(Const_km/Const_s), GasVel[c][2]*UNIT_V/(Const_km/Const_s) );
         fprintf( File_User, " %14.7e %14.7e %14.7e", RelativeVel[c]*UNIT_V/(Const_km/Const_s), SoundSpeed[c]*UNIT_V/(Const_km/Const_s), GasDens[c]*UNIT_D );
         fprintf( File_User, " %14.7e %14.7e %14d",   Bondi_MassBH[c]*UNIT_M/Const_Msun, Mdot_BH[c]*UNIT_M/UNIT_T, SinkNCell_Sum[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", MomX_Sum[c], MomY_Sum[c], MomZ_Sum[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", MomXAbs_Sum[c], MomYAbs_Sum[c], MomZAbs_Sum[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", E_inj_exp[c], E_Sum[c], (E_Sum[c]-E_inj_exp[c])/E_inj_exp[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", Ek_Sum[c], Et_Sum[c], E_power_inj[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", M_inj_exp[c], Mass_Sum[c], (Mass_Sum[c]-M_inj_exp[c])/M_inj_exp[c] );
         fprintf( File_User, " %14.7e %14.7e %14.7e", Mdot[c]*UNIT_M/UNIT_T, Pdot[c]*UNIT_M*UNIT_V/UNIT_T, Edot[c]*UNIT_E/UNIT_T );
         fprintf( File_User, " %14.7e %14.7e %14.7e", Jet_Vec[c][0], Jet_Vec[c][1], Jet_Vec[c][2] );
         fprintf( File_User, " %14d %14.7e",          num_par_sum[c], ColdGasMass[c]*UNIT_M/Const_Msun );
      }
      fprintf( File_User, "\n" );
      fclose( File_User );
   }

// reset the cumulative variables to zero
   for (int c=0; c<Merger_Coll_NumHalos; c++){
      CM_Bondi_SinkMass[c]    = 0.0;
      CM_Bondi_SinkMomX[c]    = 0.0;
      CM_Bondi_SinkMomY[c]    = 0.0;
      CM_Bondi_SinkMomZ[c]    = 0.0;
      CM_Bondi_SinkMomXAbs[c] = 0.0;
      CM_Bondi_SinkMomYAbs[c] = 0.0;
      CM_Bondi_SinkMomZAbs[c] = 0.0;
      CM_Bondi_SinkE[c]       = 0.0;
      CM_Bondi_SinkEk[c]      = 0.0;
      CM_Bondi_SinkEt[c]      = 0.0;
      E_inj_exp[c]            = 0.0;
      M_inj_exp[c]            = 0.0;
   }

} // FUNCTION : Aux_Record_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  GetClusterCenter
// Description :  Get the cluster centers
//
// Note        :  1. Must enable Merger_Coll_LabelCenter
//
// Parameter   :  Cen : Cluster centers
//
// Return      :  Cen[]
//-------------------------------------------------------------------------------------------------------
void GetClusterCenter( int lv, bool AdjustPos, bool AdjustVel, double Cen_old[][3], double Cen_new[][3], double Cen_Vel[][3] )
{

   if ( fixBH == false ){
   double min_pos[3][3], DM_Vel[3][3];   // The updated BH position / velocity 
   const bool CurrentMaxLv = (  NPatchTotal[lv] > 0  &&  ( lv == MAX_LEVEL || NPatchTotal[lv+1] == 0 )  );

// Initialize min_pos to be the old center
   for (int c=0; c<Merger_Coll_NumHalos; c++){
      for (int d=0; d<3; d++)   min_pos[c][d] = Cen_old[c][d];
   }

   if ( (CurrentMaxLv  &&  AdjustPos == true) || (CurrentMaxLv  &&  AdjustVel == true) ){

//    Do not support periodic BC
      for (int f=0; f<6; f++)                               
         if (OPT__BC_FLU[f] == BC_FLU_PERIODIC) Aux_Error( ERROR_INFO, "do not support periodic BC (OPT__BC_FLU* = 1)!\n" ); 

#     ifdef GRAVITY                                            
      if ( OPT__BC_POT == BC_POT_PERIODIC )  Aux_Error( ERROR_INFO, "do not support periodic BC (OPT__BC_POT = 1)!\n" );
#     endif


      double dis_exp = 1e-6;   // To check if the output BH positions of each calculaiton are close enough 
      bool   IfConverge = false;   // If the BH positions are close enough, then complete the calculation
      int    count = 0;   // How many times the calculation is performed (minimum: 2, maximum: 10)
      double Cen_new_pre[3][3];

      while ( IfConverge == false  &&  count <= 10 ){ 

         for (int c=0; c<Merger_Coll_NumHalos; c++){
            for (int d=0; d<3; d++)  Cen_new_pre[c][d] = min_pos[c][d];
         }

         int N = 10000;   // Maximum particle numbers (to allocate the array size)
         int num_par[3] = {0, 0, 0};   // (each rank) number of particles inside the target region of each cluster
         double **ParX = new double*[Merger_Coll_NumHalos];
         double **ParY = new double*[Merger_Coll_NumHalos];
         double **ParZ = new double*[Merger_Coll_NumHalos];
         double **ParM = new double*[Merger_Coll_NumHalos];
         double **VelX = new double*[Merger_Coll_NumHalos];
         double **VelY = new double*[Merger_Coll_NumHalos];
         double **VelZ = new double*[Merger_Coll_NumHalos];
         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            ParX[c] = new double[N];
            ParY[c] = new double[N];
            ParZ[c] = new double[N];
            ParM[c] = new double[N];
            VelX[c] = new double[N];
            VelY[c] = new double[N];
            VelZ[c] = new double[N];
         }
      
//       Find the particles within the arrection radius
         for (int c=0; c<Merger_Coll_NumHalos; c++) {   
            num_par_sum[c] = 0; 
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++) {
               const double *EdgeL = amr->patch[0][lv][PID]->EdgeL;
               const double *EdgeR = amr->patch[0][lv][PID]->EdgeR;
               const double patch_pos[3] = { (EdgeL[0]+EdgeR[0])*0.5, (EdgeL[1]+EdgeR[1])*0.5, (EdgeL[2]+EdgeR[2])*0.5 };
               const double patch_d = sqrt(SQR(EdgeL[0]-EdgeR[0])+SQR(EdgeL[1]-EdgeR[1])+SQR(EdgeL[2]-EdgeR[2]))*0.5;
               if (SQR(patch_pos[0]-Cen_new_pre[c][0])+SQR(patch_pos[1]-Cen_new_pre[c][1])+SQR(patch_pos[2]-Cen_new_pre[c][2]) <= SQR(20*R_acc+patch_d)){
                  for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++) {
                     const long ParID = amr->patch[0][lv][PID]->ParList[p];
                     const real ParX_tmp = amr->Par->PosX[ParID];
                     const real ParY_tmp = amr->Par->PosY[ParID];
                     const real ParZ_tmp = amr->Par->PosZ[ParID];
                     const real ParM_tmp = amr->Par->Mass[ParID];
                     const real VelX_tmp = amr->Par->VelX[ParID];
                     const real VelY_tmp = amr->Par->VelY[ParID];
                     const real VelZ_tmp = amr->Par->VelZ[ParID];
                     bool if_cluster = false;
                     if ( amr->Par->Type[ParID] == real(PTYPE_CLUSTER+c) || amr->Par->Type[ParID] == real(PTYPE_CEN+c) ){
                        if_cluster = true;
                     }
                     if ( if_cluster && SQR(ParX_tmp-Cen_new_pre[c][0])+SQR(ParY_tmp-Cen_new_pre[c][1])+SQR(ParZ_tmp-Cen_new_pre[c][2]) <= SQR(10*R_acc) ){
//                      Record the mass, position and velocity of this particle
                        ParX[c][num_par[c]] = ParX_tmp;
                        ParY[c][num_par[c]] = ParY_tmp;
                        ParZ[c][num_par[c]] = ParZ_tmp;
                        ParM[c][num_par[c]] = ParM_tmp;
                        VelX[c][num_par[c]] = VelX_tmp;
                        VelY[c][num_par[c]] = VelY_tmp;
                        VelZ[c][num_par[c]] = VelZ_tmp;
                        num_par[c] += 1; 
                     }
                     if (num_par[c] >= N) {
                         N = num_par[c] + 1;  // Increase the new maximum size if needed
                         ParX[c] = (double*)realloc(ParX[c], N * sizeof(double));
                         ParY[c] = (double*)realloc(ParY[c], N * sizeof(double));
                         ParZ[c] = (double*)realloc(ParZ[c], N * sizeof(double));
                         ParM[c] = (double*)realloc(ParM[c], N * sizeof(double));
                         VelX[c] = (double*)realloc(VelX[c], N * sizeof(double));
                         VelY[c] = (double*)realloc(VelY[c], N * sizeof(double));
                         VelZ[c] = (double*)realloc(VelZ[c], N * sizeof(double));
         }  }  }  }  }
      
//       Collect the number of target particles from each rank
         MPI_Allreduce( num_par, num_par_sum, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

         int num_par_eachRank[3][MPI_NRank];
         int displs[3][MPI_NRank];
         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            MPI_Allgather( &num_par[c], 1, MPI_INT, num_par_eachRank[c], 1, MPI_INT, MPI_COMM_WORLD );
            displs[c][0] = 0;
            for (int i=1; i<MPI_NRank; i++)  displs[c][i] = displs[c][i-1] + num_par_eachRank[c][i-1];
         }
      
//       Collect the mass, position and velocity of target particles to the root rank
         double **ParX_sum = new double*[Merger_Coll_NumHalos];
         double **ParY_sum = new double*[Merger_Coll_NumHalos];
         double **ParZ_sum = new double*[Merger_Coll_NumHalos];
         double **ParM_sum = new double*[Merger_Coll_NumHalos];
         double **VelX_sum = new double*[Merger_Coll_NumHalos];
         double **VelY_sum = new double*[Merger_Coll_NumHalos];
         double **VelZ_sum = new double*[Merger_Coll_NumHalos];
         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            ParX_sum[c] = new double[num_par_sum[c]];
            ParY_sum[c] = new double[num_par_sum[c]];
            ParZ_sum[c] = new double[num_par_sum[c]];
            ParM_sum[c] = new double[num_par_sum[c]];
            VelX_sum[c] = new double[num_par_sum[c]];
            VelY_sum[c] = new double[num_par_sum[c]];
            VelZ_sum[c] = new double[num_par_sum[c]];
         }
         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            MPI_Allgatherv( ParX[c], num_par[c], MPI_DOUBLE, ParX_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( ParY[c], num_par[c], MPI_DOUBLE, ParY_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( ParZ[c], num_par[c], MPI_DOUBLE, ParZ_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( ParM[c], num_par[c], MPI_DOUBLE, ParM_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( VelX[c], num_par[c], MPI_DOUBLE, VelX_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( VelY[c], num_par[c], MPI_DOUBLE, VelY_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
            MPI_Allgatherv( VelZ[c], num_par[c], MPI_DOUBLE, VelZ_sum[c], num_par_eachRank[c], displs[c], MPI_DOUBLE, MPI_COMM_WORLD );
         }

//       Compute potential and find the minimum position, and calculate the average DM velocity on the root rank
         if ( AdjustPos == true ){
            double soften = amr->dh[MAX_LEVEL];
            for (int c=0; c<Merger_Coll_NumHalos; c++) {
               double *pote = new double[num_par_sum[c]];

//             Distribute MPI jobs
               int par_per_rank = num_par_sum[c]/MPI_NRank;
               int remainder = num_par_sum[c]%MPI_NRank;
               int start = MPI_Rank*par_per_rank + MIN(MPI_Rank, remainder);
               int end = start + par_per_rank + (MPI_Rank < remainder ? 1 : 0);
               double *pote_local = new double[end-start];
    
#              pragma omp parallel for schedule(static)
               for (int i=start; i<end; i++){
                  pote_local[i-start] = 0.0;
                  for (int j=0; j<num_par_sum[c]; j++){
                     double rel_pos = sqrt(SQR(ParX_sum[c][i]-ParX_sum[c][j])+SQR(ParY_sum[c][i]-ParY_sum[c][j])+SQR(ParZ_sum[c][i]-ParZ_sum[c][j]));
                     if ( rel_pos > soften )   pote_local[i-start] += ParM_sum[c][j]/rel_pos; 
                     else if  ( rel_pos <= soften && i != j )   pote_local[i-start] += ParM_sum[c][j]/soften;
                  }
                  pote_local[i-start] *= -NEWTON_G;
               }
               int recvcounts[MPI_NRank], displs[MPI_NRank];
               for (int i=0; i<MPI_NRank; i++)  recvcounts[i] = (i < remainder ? par_per_rank+1 : par_per_rank);
               displs[0] = 0;
               for (int i=1; i<MPI_NRank; i++)  displs[i] = displs[i-1] + recvcounts[i-1]; 
               MPI_Allgatherv( pote_local, end-start, MPI_DOUBLE, pote, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

               double Pote_min = 0.0;
               for (int i=0; i<num_par_sum[c]; i++){
                  if ( pote[i] < Pote_min ){
                     Pote_min = pote[i];
                     min_pos[c][0] = ParX_sum[c][i];
                     min_pos[c][1] = ParY_sum[c][i];
                     min_pos[c][2] = ParZ_sum[c][i];
                  }
               }
               delete[] pote;
               delete[] pote_local;
            } // for (int c=0; c<Merger_Coll_NumHalos; c++)
         } // if ( AdjustPos == true )

//       Calculate the average DM velocity
         if ( AdjustVel == true ){
            for (int c=0; c<Merger_Coll_NumHalos; c++) {
               for (int d=0; d<3; d++)  DM_Vel[c][d] = 0.0;
               for (int i=0; i<num_par_sum[c]; i++){
                  DM_Vel[c][0] += VelX_sum[c][i]; 
                  DM_Vel[c][1] += VelY_sum[c][i]; 
                  DM_Vel[c][2] += VelZ_sum[c][i]; 
               }
               for (int d=0; d<3; d++)  DM_Vel[c][d] /= num_par_sum[c];
            }
         } // if ( AdjustVel == true )

//       Iterate the above calculation until the output BH positions become close enough
         count += 1;
         double dis[3] = {0.0, 0.0, 0.0};
         for (int c=0; c<Merger_Coll_NumHalos; c++){
            for (int d=0; d<3; d++)   dis[c] += SQR( min_pos[c][d] - Cen_new_pre[c][d] ); 
         }
         if ( count > 1  &&  sqrt(dis[0]) < dis_exp  &&  sqrt(dis[1]) < dis_exp )   IfConverge = true;

         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            delete[] ParX_sum[c];
            delete[] ParY_sum[c];
            delete[] ParZ_sum[c];
            delete[] ParM_sum[c];
            delete[] VelX_sum[c];
            delete[] VelY_sum[c];
            delete[] VelZ_sum[c];
         }
         delete[] ParX_sum;
         delete[] ParY_sum;
         delete[] ParZ_sum;
         delete[] ParM_sum;
         delete[] VelX_sum;
         delete[] VelY_sum;
         delete[] VelZ_sum;
         for (int c=0; c<Merger_Coll_NumHalos; c++) {
            delete[] ParX[c];
            delete[] ParY[c];
            delete[] ParZ[c];
            delete[] ParM[c];
            delete[] VelX[c];
            delete[] VelY[c];
            delete[] VelZ[c];
         }
         delete[] ParX;
         delete[] ParY;
         delete[] ParZ;
         delete[] ParM;
         delete[] VelX;
         delete[] VelY;
         delete[] VelZ;
      } // while ( IfConverge == false )
   } // if ( (CurrentMaxLv  &&  AdjustPos == true) || (CurrentMaxLv  &&  AdjustVel == true) ) 


// Find the BH particles and adjust their position and velocity
   for (int c=0; c<Merger_Coll_NumHalos; c++) {
      double Cen_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ };   // set to -inf
      double Vel_Tmp[3] = { -__FLT_MAX__, -__FLT_MAX__, -__FLT_MAX__ };
      for (long p=0; p<amr->Par->NPar_AcPlusInac; p++) {
         if ( amr->Par->Mass[p] >= (real)0.0  &&  amr->Par->Type[p] == real(PTYPE_CEN+c) ){
            if ( CurrentMaxLv  &&  AdjustPos == true ){
               amr->Par->PosX[p] = min_pos[c][0];
               amr->Par->PosY[p] = min_pos[c][1];
               amr->Par->PosZ[p] = min_pos[c][2];    
            }
            if ( CurrentMaxLv  &&  AdjustVel == true ){
               amr->Par->VelX[p] = DM_Vel[c][0];
               amr->Par->VelY[p] = DM_Vel[c][1];
               amr->Par->VelZ[p] = DM_Vel[c][2];        
            }   
            Cen_Tmp[0] = amr->Par->PosX[p];
            Cen_Tmp[1] = amr->Par->PosY[p];
            Cen_Tmp[2] = amr->Par->PosZ[p];
            Vel_Tmp[0] = amr->Par->VelX[p];
            Vel_Tmp[1] = amr->Par->VelY[p];
            Vel_Tmp[2] = amr->Par->VelZ[p];
            break;
         }
      }

//    use MPI_MAX since Cen_Tmp[] is initialized as -inf
      MPI_Allreduce( Cen_Tmp, Cen_new[c], 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
      MPI_Allreduce( Vel_Tmp, Cen_Vel[c], 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
   } // for (int c=0; c<Merger_Coll_NumHalos; c++)

   if ( CurrentMaxLv  &&  AdjustPos == true ){
      const bool TimingSendPar_No = false;
      Par_PassParticle2Sibling( lv, TimingSendPar_No );
      Par_PassParticle2Son_MultiPatch( lv, PAR_PASS2SON_EVOLVE, TimingSendPar_No, NULL_INT, NULL );
   }
   } // if ( fixBH == false )
   else {   // if ( fixBH == true )
      for (int d=0; d<3; d++)   Cen_new[0][d] = 7.5; 
      for (int d=0; d<3; d++)   Cen_Vel[0][d] = 0.0; 
   }

   for (int c=0; c<Merger_Coll_NumHalos; c++) {
      for (int d=0; d<3; d++)  Cen_old[c][d] = Cen_new[c][d];
   }
} // FUNCTION : GetClusterCenter

#endif // #ifdef MASSIVE_PARTICLES
