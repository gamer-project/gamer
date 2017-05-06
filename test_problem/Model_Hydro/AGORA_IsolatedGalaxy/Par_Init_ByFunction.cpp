#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE



extern bool CheckEmptyString( const char *InputString );
extern int  CountRow( const char *Filename );

extern char AGORA_HaloPar_Filename[1000];
extern char AGORA_DiskPar_Filename[1000];
extern char AGORA_BulgePar_Filename[1000];




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
   real (*ParData_AllRank)[NParVar] = NULL;

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

      NPar_Halo  = CountRow( AGORA_HaloPar_Filename );
      NPar_Disk  = CountRow( AGORA_DiskPar_Filename );
      NPar_Bulge = CountRow( AGORA_BulgePar_Filename );
      NPar_All   = NPar_Halo + NPar_Disk + NPar_Bulge;

      Aux_Message( stdout, "   Number of halo  particles = %ld\n", NPar_Halo );
      Aux_Message( stdout, "   Number of disk  particles = %ld\n", NPar_Disk );
      Aux_Message( stdout, "   Number of bulge particles = %ld\n", NPar_Bulge );
      Aux_Message( stdout, "   Number of all   particles = %ld\n", NPar_All );

      if ( NPar_All != amr->Par->NPar_Active_AllRank )
         Aux_Error( ERROR_INFO, "total number of particles found in disk [%ld] != expect [%ld] !!\n",
                    NPar_All, amr->Par->NPar_Active_AllRank );


//    allocate memory to store all particles loaded from disk
      ParData_AllRank = new real [NPar_All][NParVar];


//    load data
      const char *Filename[3] = { AGORA_HaloPar_Filename, AGORA_DiskPar_Filename, AGORA_BulgePar_Filename };
      const int   MaxLine     = 1024;  // maximum number of characters per line

      char  *Line      = new char [MaxLine];
      char  *FirstChar = NULL;
      FILE  *File      = NULL;
      long   NLoad     = 0;

//    loop over all three particle tables
      for (int t=0; t<3; t++)
      {
         Aux_Message( stdout, "   Loading particles from the file \"%s\" ... ", Filename[t] );

         File = fopen( Filename[t], "r" );

         while ( fgets( Line, MaxLine, File ) != NULL )
         {
//          skip empty lines
            if (  !CheckEmptyString( Line )  )
            {
               FirstChar = Line;

//             find the first non-empty character
               while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//             skip lines starting with "#"
               if ( *FirstChar != '#' )
               {
                  if ( NLoad >= NPar_All )   Aux_Error( ERROR_INFO, "NLoad (%ld) >= NPar_All (%ld) !!\n", NLoad, NPar_All );

//                note that the on-disk data format is (x, y, z, vx, vy, vz, mass),
//                while here we want to convert it to  (mass, x, y, z, vx, vy, vz)
#                 ifdef FLOAT8
                  sscanf( Line, "%lf%lf%lf%lf%lf%lf%lf",
                          ParData_AllRank[NLoad]+1, ParData_AllRank[NLoad]+2, ParData_AllRank[NLoad]+3,
                          ParData_AllRank[NLoad]+4, ParData_AllRank[NLoad]+5, ParData_AllRank[NLoad]+6,
                          ParData_AllRank[NLoad]+0 );
#                 else
                  sscanf( Line, "%f%f%f%f%f%f%f",
                          ParData_AllRank[NLoad]+1, ParData_AllRank[NLoad]+2, ParData_AllRank[NLoad]+3,
                          ParData_AllRank[NLoad]+4, ParData_AllRank[NLoad]+5, ParData_AllRank[NLoad]+6,
                          ParData_AllRank[NLoad]+0 );
#                 endif

                  NLoad ++;
               }
            } // if (  !CheckEmptyString( Line )  )
         } // while ( fgets( Line, MaxLine, File ) != NULL )

         fclose( File );

         Aux_Message( stdout, "done\n" );

      } //for (int t=0; t<3; t++)

      delete [] Line;

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
   MPI_Scatterv( ParData_AllRank[0], NSend, SendDisp, MPI_DOUBLE, ParData_MyRank[0], NPar_MyRank*NParVar, MPI_DOUBLE, 0, MPI_COMM_WORLD );
#  else
   MPI_Scatterv( ParData_AllRank[0], NSend, SendDisp, MPI_FLOAT,  ParData_MyRank[0], NPar_MyRank*NParVar, MPI_FLOAT,  0, MPI_COMM_WORLD );
#  endif

   delete [] ParData_AllRank;

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


// store data in the particle repository
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Storing particles in the particle repository ... " );

   for (int p=0; p<NPar_MyRank; p++)
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

   const real BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };

   for (int p=0; p<NPar_MyRank; p++)
   {
      amr->Par->PosX[p] += BoxCenter[0];
      amr->Par->PosY[p] += BoxCenter[1];
      amr->Par->PosZ[p] += BoxCenter[2];
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction



#endif // #ifdef PARTICLE
