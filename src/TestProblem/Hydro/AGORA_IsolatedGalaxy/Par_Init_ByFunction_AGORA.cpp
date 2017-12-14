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
// Note        :  1. Invoked by "Init_GAMER" using the function pointer "Par_Init_ByFunction_Ptr"
//                2. Particles lying outside the active region will be removed by "Par_Aux_InitCheck"
//                   if non-periodic B.C. is adopted
//                3. Particles set here are only temporarily stored in this rank
//                   --> They will be redistributed when calling "Par_LB_Init_RedistributeByRectangular()
//                       and LB_Init_LoadBalance()"
//                4. For LOAD_BALANCE, the number of particles in each rank must be set in advance
//                   --> Currently it's set by "Init_Parallelization()" and stored in "amr->Par->NPar_AcPlusInac"
//                5. Input text file format: (x, y, z, vx, vy, vz, mass)
//
// Parameter   :  None
//
// Return      :  amr->Par->Mass/Pos(X/Y/Z)/Vel(X/Y/Z)/Time
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_AGORA()
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
      long NPar_Halo, NPar_Disk, NPar_Bulge, NPar_All;

      NPar_Halo  = Aux_CountRow( AGORA_HaloPar_Filename );
      NPar_Disk  = Aux_CountRow( AGORA_DiskPar_Filename );
      NPar_Bulge = Aux_CountRow( AGORA_BulgePar_Filename );
      NPar_All   = NPar_Halo + NPar_Disk + NPar_Bulge;

      Aux_Message( stdout, "   Number of halo  particles = %ld\n", NPar_Halo );
      Aux_Message( stdout, "   Number of disk  particles = %ld\n", NPar_Disk );
      Aux_Message( stdout, "   Number of bulge particles = %ld\n", NPar_Bulge );
      Aux_Message( stdout, "   Number of all   particles = %ld\n", NPar_All );

      if ( NPar_All != amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                    NPar_All, amr->Par->NPar_Active_AllRank );


//    allocate memory to store all particles loaded from disk
      ParData_AllRank = new real [NPar_All*NParVar];


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
      if ( NPar_Loaded != NPar_All )
         Aux_Error( ERROR_INFO, "total number of loaded particles (%ld) != expect (%ld) !!\n", NPar_Loaded, NPar_All );

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
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Sending particles from the root rank to all ranks ... " );
   real (*ParData_MyRank)[NParVar] = new real [NPar_MyRank][NParVar];

#  ifdef FLOAT8
   MPI_Scatterv( ParData_AllRank, NSend, SendDisp, MPI_DOUBLE, ParData_MyRank[0], NPar_MyRank*NParVar, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Scatterv( ParData_AllRank, NSend, SendDisp, MPI_FLOAT,  ParData_MyRank[0], NPar_MyRank*NParVar, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   delete [] ParData_AllRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data in the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing particles in the particle repository ... " );

   for (int p=0; p<NPar_MyRank; p++)
   {
//    note that the on-disk data format is (x, y, z, vx, vy, vz, mass)
      amr->Par->PosX[p] = ParData_MyRank[p][0];
      amr->Par->PosY[p] = ParData_MyRank[p][1];
      amr->Par->PosZ[p] = ParData_MyRank[p][2];
      amr->Par->VelX[p] = ParData_MyRank[p][3];
      amr->Par->VelY[p] = ParData_MyRank[p][4];
      amr->Par->VelZ[p] = ParData_MyRank[p][5];
      amr->Par->Mass[p] = ParData_MyRank[p][6];

//    synchronize all particles to the physical time at the base level
      amr->Par->Time[p] = Time[0];
   }

   delete [] ParData_MyRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// convert to code units (assuming the loaded data have the units 1e9 Msun, kpc, and km/s)
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Converting particle attributes to the code units ... " );

   const real UnitM = 1.0e9*Const_Msun/UNIT_M;
   const real UnitL = 1.0e0*Const_kpc /UNIT_L;
   const real UnitV = 1.0e0*Const_km  /UNIT_V;

   for (int p=0; p<NPar_MyRank; p++)
   {
      amr->Par->Mass[p] *= UnitM;
      amr->Par->PosX[p] *= UnitL;
      amr->Par->PosY[p] *= UnitL;
      amr->Par->PosZ[p] *= UnitL;
      amr->Par->VelX[p] *= UnitV;
      amr->Par->VelY[p] *= UnitV;
      amr->Par->VelZ[p] *= UnitV;
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// shift center (assuming the center of loaded particles = [0,0,0])
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Shifting particle center ... " );

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   for (int p=0; p<NPar_MyRank; p++)
   {
      amr->Par->PosX[p] += BoxCenter[0];
      amr->Par->PosY[p] += BoxCenter[1];
      amr->Par->PosZ[p] += BoxCenter[2];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// initialize the passive particle attributes, specifically, PAR_CREATION_TIME and PAR_METAL_FRAC
// --> set to an arbitrary negative value to indicate that they are actually useless since they are only used
//     for star particles created during the evolution
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Initializing passive particle attributes ... " );

   const real Useless = -1.0;

#  ifdef STAR_FORMATION
   for (int p=0; p<NPar_MyRank; p++)   amr->Par->Passive[PAR_CREATION_TIME][p] = Useless;
#  endif

//###: HARD-CODED FIELDS
#  if ( PAR_NPASSIVE > 0 )
   if ( AGORA_UseMetal )
   for (int p=0; p<NPar_MyRank; p++)   amr->Par->Passive[PAR_METAL_FRAC   ][p] = Useless;
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_AGORA



#endif // #ifdef PARTICLE
