#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_Record_ParticleCount
// Description :  Count the number of active particles at each level in each rank
//
// Note        :  1. OPT__PARTICLE_COUNT = 1/2 --> count the number of particles every step/sub-step
//                2. Output filename is "Record__ParticleCount"
//-------------------------------------------------------------------------------------------------------
void Par_Aux_Record_ParticleCount()
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__ParticleCount", OUTPUT_DIR );

   if ( MPI_Rank == 0  &&  FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      FirstTime = false;
   }


// 1. gather information from all ranks
   const long NPar_Tot = amr->Par->NPar_Active_AllRank;

   long   NPar_AllRank[NLEVEL];
   double Frac_AllRank[NLEVEL];

   long   (*NPar_EachRank)[NLEVEL] = new long   [MPI_NRank][NLEVEL];
   double (*Frac_EachRank)[NLEVEL] = new double [MPI_NRank][NLEVEL];

   MPI_Gather( amr->Par->NPar_Lv, NLEVEL, MPI_LONG, NPar_EachRank[0], NLEVEL, MPI_LONG, 0, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )
   {
//    2. calculate the fraction of particles at each level in each rank
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       get the total number of particles at each level in all ranks
         NPar_AllRank[lv] = 0;
         for (int r=0; r<MPI_NRank; r++)  NPar_AllRank[lv] += NPar_EachRank[r][lv];

//       get the fraction of particles at each level in each rank
         for (int r=0; r<MPI_NRank; r++)  Frac_EachRank[r][lv] = 100.0*NPar_EachRank[r][lv]/NPar_Tot;

//       get the fraction of particles at each level in all ranks
         Frac_AllRank[lv] = 100.0*NPar_AllRank[lv]/NPar_Tot;
      }


//    3. write to the file "Record__ParticleCount"
      FILE *File = fopen( FileName, "a" );

      fprintf( File, "Time = %13.7e,  Step = %7ld,  NPar = %10ld\n\n", Time[0], Step, NPar_Tot );

      fprintf( File, "%4s", "Rank" );
      for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%15s %-2d", "Level", lv );
      fprintf( File, "\n" );

      for (int r=0; r<MPI_NRank; r++)
      {
         fprintf( File, "%4d", r );
         for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%9ld(%6.2lf%%)", NPar_EachRank[r][lv], Frac_EachRank[r][lv] );
         fprintf( File, "\n" );
      }

      fprintf( File, "-------------------------------------------------------------------------------------" );
      fprintf( File, "-------------------------------------------------------------------------------------\n" );
      fprintf( File, "%4s", "Sum:" );
      for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, "%9ld(%6.2lf%%)", NPar_AllRank[lv], Frac_AllRank[lv] );
      fprintf( File, "\n" );
      fprintf( File, "%4s", "Ave:" );
      for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, "%12.2f%6s", (double)NPar_AllRank[lv]/MPI_NRank, "" );
      fprintf( File, "\n" );


//    4. get the load and load-imbalance factor at each level
      long   NPar_Max[NLEVEL];
      double NPar_Ave[NLEVEL], NPar_Imb[NLEVEL];   // Imb : Imbalance
      double WLoad_Ave=0.0, WLoad_Max=0.0, WLI;    // WLoad : weighted load

//    4-1. estimate the load-imbalance factor at each level
      for (int lv=0; lv<NLEVEL; lv++)
      {
         NPar_Max[lv] = -1;
         NPar_Ave[lv] = (double)NPar_AllRank[lv] / MPI_NRank;

         for (int r=0; r<MPI_NRank; r++)
            if ( NPar_EachRank[r][lv] > NPar_Max[lv] )   NPar_Max[lv] = NPar_EachRank[r][lv];

         NPar_Imb[lv] = (NPar_Max[lv]==0) ? 0.0 : ( NPar_Max[lv] - NPar_Ave[lv] ) / NPar_Ave[lv];
      }

//    4-2. estimate the weighted load-imbalance factor of particles at all levels
      for (int lv=0; lv<NLEVEL; lv++)
      {
         WLoad_Max += (double)amr->NUpdateLv[lv]*NPar_Max[lv];
         WLoad_Ave += (double)amr->NUpdateLv[lv]*NPar_Ave[lv];
      }

      WLI = ( WLoad_Max - WLoad_Ave ) / WLoad_Ave;

//    4-3. record the load-imbalance factors
      fprintf( File, "%4s", "Imb:" );
      for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%9ld(%6.2f%%)", NPar_Max[lv], 100.0*NPar_Imb[lv] );
      fprintf( File, "\n" );
      fprintf( File, "Weighted load-imbalance factor = %6.2f%%\n", 100.0*WLI );

      fprintf( File, "-------------------------------------------------------------------------------------" );
      fprintf( File, "-------------------------------------------------------------------------------------\n" );
      fprintf( File, "\n\n" );

      fclose( File );
   } // if ( MPI_Rank == 0 )


   delete [] NPar_EachRank;
   delete [] Frac_EachRank;

} // FUNCTION : Par_Aux_Record_ParticleCount



#endif // #ifdef PARTICLE
