#include "GAMER.h"

#ifdef PARTICLE


extern char AGORA_HaloPar_Filename[1000];
extern char AGORA_DiskPar_Filename[1000];
extern char AGORA_BulgePar_Filename[1000];
extern bool AGORA_UseMetal;




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_AGORA
// Description :  Initialize all particle attributes for the AGORA isolated galaxy simulation
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
//                4. Input particle IC text file format: (x, y, z, vx, vy, vz, mass)
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
void Par_Init_ByFunction_AGORA( const long NPar_ThisRank, const long NPar_AllRank,
                                real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                real *ParPassive[PAR_NPASSIVE] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const int NParVar = 7;  // mass, pos*3, vel*3
   real *ParData_AllRank = NULL;

// load data --> for simplicity, currently only the root rank will load data from disk
   if ( MPI_Rank == 0 )
   {
//    check file existence
      if ( !Aux_CheckFileExist(AGORA_HaloPar_Filename) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", AGORA_HaloPar_Filename );

      if ( !Aux_CheckFileExist(AGORA_DiskPar_Filename) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", AGORA_DiskPar_Filename );

      if ( !Aux_CheckFileExist(AGORA_BulgePar_Filename) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", AGORA_BulgePar_Filename );


//    check file size
      long NPar_Halo, NPar_Disk, NPar_Bulge, NPar_Sum;

      NPar_Halo  = Aux_CountRow( AGORA_HaloPar_Filename );
      NPar_Disk  = Aux_CountRow( AGORA_DiskPar_Filename );
      NPar_Bulge = Aux_CountRow( AGORA_BulgePar_Filename );
      NPar_Sum   = NPar_Halo + NPar_Disk + NPar_Bulge;

      Aux_Message( stdout, "   Number of halo  particles = %ld\n", NPar_Halo );
      Aux_Message( stdout, "   Number of disk  particles = %ld\n", NPar_Disk );
      Aux_Message( stdout, "   Number of bulge particles = %ld\n", NPar_Bulge );
      Aux_Message( stdout, "   Sum    of all   particles = %ld\n", NPar_Sum );

      if ( NPar_Sum != NPar_AllRank )
         Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                    NPar_Sum, NPar_AllRank );


//    allocate memory to store all particles loaded from disk
      ParData_AllRank = new real [NPar_Sum*NParVar];


//    load data from the three particle tables
      const char *Filename[3]  = { AGORA_HaloPar_Filename, AGORA_DiskPar_Filename, AGORA_BulgePar_Filename };
      const bool  RowMajor_Yes = true;                   // load data into the row-major order
      const bool  AllocMem_No  = false;                  // do not allocate memory for ParData_AllRank
      const int   NCol         = 7;                      // total number of columns to load
      const int   Col[NCol]    = {0, 1, 2, 3, 4, 5, 6};  // target columns: (x, y, z, vx, vy, vz, mass)

      long NPar_Loaded = 0;

      for (int t=0; t<3; t++)
      {
         Aux_Message( stdout, "   Loading particles from the file \"%s\" ... ", Filename[t] );

//       must use a temporary pointer "tmp_ptr" for Aux_LoadTable() because of the call-by-reference approach
         real *tmp_ptr = ParData_AllRank + NPar_Loaded*NParVar;

         NPar_Loaded += Aux_LoadTable( tmp_ptr, Filename[t], NCol, Col, RowMajor_Yes, AllocMem_No );

         Aux_Message( stdout, "done\n" );
      }

//    check
      if ( NPar_Loaded != NPar_Sum )
         Aux_Error( ERROR_INFO, "total number of loaded particles (%ld) != expect (%ld) !!\n", NPar_Loaded, NPar_Sum );

   } // if ( MPI_Rank == 0 )

   MPI_Barrier( MPI_COMM_WORLD );


// get the number of particles in each rank and set the corresponding offsets
   if ( NParVar*NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "Total number of particle attributes to be sent (%ld) exceeds the maximum integer (%ld) !!\n",
                 (long)NParVar*NPar_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];

   MPI_Gather( &NPar_ThisRank, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)  NSend[r] *= NParVar;

      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particles from the root rank to all ranks
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Sending particles from the root rank to all ranks ... " );
   real (*ParData_MyRank)[NParVar] = new real [NPar_ThisRank][NParVar];

#  ifdef FLOAT8
   MPI_Scatterv( ParData_AllRank,   NSend, SendDisp,       MPI_DOUBLE,
                 ParData_MyRank[0], NPar_ThisRank*NParVar, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Scatterv( ParData_AllRank,   NSend, SendDisp,        MPI_FLOAT,
                 ParData_MyRank[0], NPar_ThisRank*NParVar,  MPI_FLOAT, 0, MPI_COMM_WORLD );
#  endif

   delete [] ParData_AllRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data in the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing particles in the particle repository ... " );

   for (int p=0; p<NPar_ThisRank; p++)
   {
//    note that the on-disk data format is (x, y, z, vx, vy, vz, mass)
      ParPosX[p] = ParData_MyRank[p][0];
      ParPosY[p] = ParData_MyRank[p][1];
      ParPosZ[p] = ParData_MyRank[p][2];
      ParVelX[p] = ParData_MyRank[p][3];
      ParVelY[p] = ParData_MyRank[p][4];
      ParVelZ[p] = ParData_MyRank[p][5];
      ParMass[p] = ParData_MyRank[p][6];

//    synchronize all particles to the physical time at the base level
      ParTime[p] = Time[0];
   }

   delete [] ParData_MyRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// convert to code units (assuming the loaded data have the units 1e9 Msun, kpc, and km/s)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Converting particle attributes to the code units ... " );

   const real UnitM = 1.0e9*Const_Msun/UNIT_M;
   const real UnitL = 1.0e0*Const_kpc /UNIT_L;
   const real UnitV = 1.0e0*Const_km  /UNIT_V;

   for (int p=0; p<NPar_ThisRank; p++)
   {
      ParMass[p] *= UnitM;
      ParPosX[p] *= UnitL;
      ParPosY[p] *= UnitL;
      ParPosZ[p] *= UnitL;
      ParVelX[p] *= UnitV;
      ParVelY[p] *= UnitV;
      ParVelZ[p] *= UnitV;
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center ... " );

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   for (int p=0; p<NPar_ThisRank; p++)
   {
      ParPosX[p] += BoxCenter[0];
      ParPosY[p] += BoxCenter[1];
      ParPosZ[p] += BoxCenter[2];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// initialize the passive particle attributes, specifically, PAR_CREATION_TIME and PAR_METAL_FRAC
// --> set to an arbitrary negative value to indicate that they are actually useless since they are only used
//     for star particles created during the evolution
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing passive particle attributes ... " );

   const real Useless = -1.0;

#  ifdef STAR_FORMATION
   for (int p=0; p<NPar_ThisRank; p++)    ParPassive[PAR_CREATION_TIME][p] = Useless;
#  endif

//###: HARD-CODED FIELDS
#  if ( PAR_NPASSIVE > 0 )
   if ( AGORA_UseMetal )
   for (int p=0; p<NPar_ThisRank; p++)    ParPassive[PAR_METAL_FRAC   ][p] = Useless;
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_AGORA



#endif // #ifdef PARTICLE
