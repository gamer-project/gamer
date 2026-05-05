#include "GAMER.h"
#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_DiskHeating
// Description :  User-specified function to initialize particle attributes
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
//                4. The initialization of the PUID routine has been separated into amr->Par->InitRepo()
//                   --> If needed, you can still modify PUID through the AllAttributeInt array
//
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_DiskHeating( const long NPar_ThisRank, const long NPar_AllRank,
                                      real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                      real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                      long_par *ParType,
                                      real_par *AllAttribute[PAR_NATT_FLT_TOTAL],
                                      long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// NOTE: DiskHeatingParticleIC uses the floating-point type for particle type and assumes single precision
   const char FileName[]   = "DiskHeatingParticleIC";
   const long NParAllRank  = amr->Par->NPar_Active_AllRank;
         long NParThisRank = amr->Par->NPar_AcPlusInac;        // cannot be "const" due to MPI_Allgather()
   const int  NParAtt      = 8;                                // mass, pos*3, vel*3, type


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *FileTemp = fopen( FileName, "rb" );

   fseek( FileTemp, 0, SEEK_END );

   const long ExpectSize = long(NParAtt)*NParAllRank*sizeof(real);
   const long FileSize   = ftell( FileTemp );
   if ( FileSize != ExpectSize )
      Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                 FileName, FileSize, ExpectSize );

   fclose( FileTemp );

   MPI_Barrier( MPI_COMM_WORLD );


// set the file offset for this rank
   long NPar_EachRank[MPI_NRank], NPar_Check=0, FileOffset=0;

   MPI_Allgather( &NParThisRank, 1, MPI_LONG, NPar_EachRank, 1, MPI_LONG, MPI_COMM_WORLD );

// check if the total number of particles is correct
   for (int r=0; r<MPI_NRank; r++)  NPar_Check += NPar_EachRank[r];
   if ( NPar_Check != NParAllRank )
      Aux_Error( ERROR_INFO, "total number of particles found (%ld) != expect (%ld) !!\n", NPar_Check, NParAllRank );

   for (int r=0; r<MPI_Rank; r++)
   {
      FileOffset = FileOffset + NPar_EachRank[r]*sizeof(real);
   }


// load data
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Loading data ... " );

   real *ParData_ThisRank = new real [ NParThisRank*NParAtt ];

// note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
   FILE *File = fopen( FileName, "rb" );

   for (int v=0; v<NParAtt; v++)
   {
      fseek( File, FileOffset+v*NParAllRank*sizeof(real), SEEK_SET );
      fread( ParData_ThisRank+v*NParThisRank, sizeof(real), NParThisRank, File );
   }

   fclose( File );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data into the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing data into particle repository ... " );

   real *ParData1 = new real [NParAtt];

   for (long p=0; p<NParThisRank; p++)
   {
//    collect data for the target particle
//    [att][id]
      for (int v=0; v<NParAtt; v++)
         ParData1[v] = ParData_ThisRank[ v*NParThisRank + p ];

//    assuming that the orders of the particle attributes stored on the disk and in Par->Attribute[] are the same
      ParMass[p] = ParData1[0];
      ParPosX[p] = ParData1[1];
      ParPosY[p] = ParData1[2];
      ParPosZ[p] = ParData1[3];
      ParVelX[p] = ParData1[4];
      ParVelY[p] = ParData1[5];
      ParVelZ[p] = ParData1[6];
      ParType[p] = (long_par)ParData1[7]; // 1=CDM halo, 2=disk

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_ThisRank;
   delete [] ParData1;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_DiskHeating



#endif // #ifdef PARTICLE
