#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_EstimateLoadImbalance
// Description :  Estimate the current load-imbalance factor
//
// Note        :  1. Weighted load-imbalance (WLI) factor is estimated as "(Load_Max - Load_Ave)/Load_Ave",
//                   where Load_Max and Load_Ave are the maximum and average workload in all ranks
//                   --> WLI =   0.0% --> perfect balance
//                           = 100.0% --> estimated performance is only half of the maximum performance
//                2. Weighting at each level is assumed to be equal to "amr->NUpdateLv"
//                3. Call LB_EstimateWorkload_AllPatchGroup() to get the workload of all patch groups at a given level
//                   --> Note that LB_EstimateWorkload_AllPatchGroup() takes into account particles in the
//                       children patches
//                       --> WLI estimated here will be different from both Record__PatchCount and
//                           Record__ParticleCount. The latter only considers particles in the leaf patches
//                4. Invoked by main() to determine whether we should redistribute all patches
//                   (by calling LB_Init_LoadBalance()) to improve the load balance
//
// Return      :  amr->LB->WLI
//-------------------------------------------------------------------------------------------------------
double LB_EstimateLoadImbalance()
{

// 1. get the workload at each level for each rank
#  ifdef PARTICLE
   const double ParWeight = amr->LB->Par_Weight;
#  else
   const double ParWeight = 0.0;
#  endif

   double Load_ThisRank[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const int NPG = amr->NPatchComma[lv][1] / 8;

      double *Load_AllPG = new double [NPG];

      LB_EstimateWorkload_AllPatchGroup( lv, ParWeight, Load_AllPG );

      Load_ThisRank[lv] = 0.0;
      for (int t=0; t<NPG; t++)  Load_ThisRank[lv] += Load_AllPG[t];

//    multiply the weighting at different levels
      Load_ThisRank[lv] *= (double)amr->NUpdateLv[lv];

      delete [] Load_AllPG;
   }


// 2. collect the workload from all ranks
   double (*Load_AllRank)[NLEVEL] = ( MPI_Rank == 0 ) ? new double [MPI_NRank][NLEVEL] : NULL;

   MPI_Gather( Load_ThisRank, NLEVEL, MPI_DOUBLE, Load_AllRank, NLEVEL, MPI_DOUBLE, 0, MPI_COMM_WORLD );


   if ( MPI_Rank == 0 )
   {
//    3. estimate the weighted load-imbalance (WLI) factor
      double Load_Max[NLEVEL], Load_Ave[NLEVEL], Load_Imb[NLEVEL];   // at each level (Imb : Imbalance)
      double Load_Ave_AllLv, Load_Max_AllLv;                         // weighted sum over all levels

//    3-1. each level
      for (int lv=0; lv<NLEVEL; lv++)
      {
         Load_Max[lv] = -1.0;
         Load_Ave[lv] = 0.0;

         for (int r=0; r<MPI_NRank; r++)
         {
            Load_Max[lv]  = MAX( Load_Max[lv], Load_AllRank[r][lv] );
            Load_Ave[lv] += Load_AllRank[r][lv];
         }

         Load_Ave[lv] /= (double)MPI_NRank;
         Load_Imb[lv]  = ( Load_Max[lv] == 0.0 ) ? 0.0 : ( Load_Max[lv] - Load_Ave[lv] ) / Load_Ave[lv];
      }

//    3-2. weighted sum over all levels
      Load_Ave_AllLv = 0.0;
      Load_Max_AllLv = 0.0;

      for (int lv=0; lv<NLEVEL; lv++)
      {
         Load_Max_AllLv += Load_Max[lv];
         Load_Ave_AllLv += Load_Ave[lv];
      }

      amr->LB->WLI = ( Load_Max_AllLv - Load_Ave_AllLv ) / Load_Ave_AllLv;


//    4. write to the file "Record__LoadBalance"
      if ( OPT__RECORD_LOAD_BALANCE )
      {
         const char FileName[] = "Record__LoadBalance";
         static bool FirstTime = true;

         if ( FirstTime )
         {
            if ( Aux_CheckFileExist(FileName) )
               Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

            FirstTime = false;
         }

//       get the total number of patches and particles
         long NPatchAll=0, NParAll=0;
         for (int lv=0; lv<NLEVEL; lv++)  NPatchAll += NPatchTotal[lv];
#        ifdef PARTICLE
         NParAll = amr->Par->NPar_Active_AllRank;
#        endif

         FILE *File = fopen( FileName, "a" );

         fprintf( File, "Time %13.7e,  Step %7ld,  NPatch %12ld,  NPar %12ld\n\n",
                  Time[0], Step, NPatchAll, NParAll );

         fprintf( File, "%4s", "Rank" );
         for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, "%16s %-2d", "Level", lv );
         fprintf( File, "\n" );

         for (int r=0; r<MPI_NRank; r++)
         {
            fprintf( File, "%4d", r );
            for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, " %8.2e(%+7.2lf%%)", Load_AllRank[r][lv],
                                                      (Load_Ave[lv]==0.0)?0.0:100.0*(Load_AllRank[r][lv]-Load_Ave[lv])/Load_Ave[lv] );
            fprintf( File, "\n" );
         }

         fprintf( File, "-------------------------------------------------------------------------------------" );
         fprintf( File, "-------------------------------------------------------------------------------------\n" );

         fprintf( File, "%4s", "Sum:" );
         for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, " %8.2e%10s", Load_Ave[lv]*MPI_NRank, "" );
         fprintf( File, "\n" );

         fprintf( File, "%4s", "Ave:" );
         for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, " %8.2e%10s", Load_Ave[lv], "" );
         fprintf( File, "\n" );

         fprintf( File, "%4s", "Max:" );
         for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, " %8.2e%10s", Load_Max[lv], "" );
         fprintf( File, "\n" );

         fprintf( File, "%4s", "Imb:" );
         for (int lv=0; lv<NLEVEL; lv++)     fprintf( File, " %7.2lf%%%10s", 100.0*Load_Imb[lv], "" );
         fprintf( File, "\n" );

         fprintf( File, "Weighted load-imbalance factor = %6.2f%%\n", 100.0*amr->LB->WLI );

         fprintf( File, "-------------------------------------------------------------------------------------" );
         fprintf( File, "-------------------------------------------------------------------------------------\n" );
         fprintf( File, "\n\n" );

         fclose( File );
      } // if ( OPT__RECORD_LOAD_BALANCE )
   } // if ( MPI_Rank == 0 )


// broadcast WLI to all ranks
   MPI_Bcast( &amr->LB->WLI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// free memory
   if ( MPI_Rank == 0 )    delete [] Load_AllRank;


   return amr->LB->WLI;

} // FUNCTION : LB_EstimateLoadImbalance



#endif // #ifdef LOAD_BALANCE
