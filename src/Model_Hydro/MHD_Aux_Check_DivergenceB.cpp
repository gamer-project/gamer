#include "GAMER.h"

#ifdef MHD


// tolerated errors
#ifdef FLOAT8
#  define ERROR_TOLERANCE  1.0e-11
#else
#  define ERROR_TOLERANCE  1.0e-5f
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_Aux_Check_DivergenceB
// Description :  Verify the divergence of the magnetic field
//
// Note        :  1. Tolerated error is set by the symbolic constant ERROR_TOLERANCE
//                   --> Should be close to the round-off errors
//                2. Invoked by Aux_Check()
//                3. Results are recorded in the file "Record__DivB"
//
// Parameter   :  Verbose : Output warning messages for cells with large errors
//                comment : You can put the location where this function is invoked in this string
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MHD_Aux_Check_DivergenceB( const bool Verbose, const char *comment )
{

   real   DivB;
   double DivB_Ave=0.0, DivB_Max=0.0;
   long   NAve=0, NFailed=0;
   int    Pass=true;


// loop over all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
      {
         if ( MPI_Rank == TargetRank )
         {
            const int MagSg = amr->MagSg[lv];

//          we check both leaf and non-leaf patches
            for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
            {
               for (int k=0; k<PS1; k++)
               for (int j=0; j<PS1; j++)
               for (int i=0; i<PS1; i++)
               {
//                compute div(B)
                  DivB = MHD_GetCellCenteredDivBInPatch( lv, PID, i, j, k, MagSg );

//                record the maximum error
                  DivB_Max = MAX( DivB_Max, DivB );

//                record errors larger than the tolerance
                  if ( DivB > ERROR_TOLERANCE )
                  {
                     NFailed ++;

                     if ( Verbose )
                     {
                        if ( Pass )
                        {
                           Aux_Message( stderr, "\"%s\" : <%s> FAILED, Time = %13.7e, Step = %ld !!\n",
                                        comment, __FUNCTION__, Time[0], Step );
                           Aux_Message( stderr, "%4s %2s %7s %2s %2s %2s %21s %21s %21s %21s %21s %21s %21s\n",
                                        "Rank", "Lv", "PID", "i", "j", "k", "BxL", "BxR", "ByL", "ByR", "BzL", "BzR", "Error" );
                           Pass = false;
                        }

//                      record B field
                        const real (*B)[PS1P1*PS1*PS1] = amr->patch[MagSg][lv][PID]->magnetic;

                        const int idx_BxL = IDX321_BX( i, j, k, PS1, PS1 );
                        const int idx_ByL = IDX321_BY( i, j, k, PS1, PS1 );
                        const int idx_BzL = IDX321_BZ( i, j, k, PS1, PS1 );

                        real BxL, BxR, ByL, ByR, BzL, BzR;

                        BxL = B[MAGX][ idx_BxL            ];
                        BxR = B[MAGX][ idx_BxL + 1        ];
                        ByL = B[MAGY][ idx_ByL            ];
                        ByR = B[MAGY][ idx_ByL + PS1      ];
                        BzL = B[MAGZ][ idx_BzL            ];
                        BzR = B[MAGZ][ idx_BzL + SQR(PS1) ];

                        Aux_Message( stderr, "%4d %2d %7d %2d %2d %2d %21.14e %21.14e %21.14e %21.14e %21.14e %21.14e %21.14e\n",
                                     MPI_Rank, lv, PID, i, j, k, BxL, BxR, ByL, ByR, BzL, BzR, DivB );
                     } // if ( Verbose )
                  } // if ( DivB > ERROR_TOLERANCE )

//                record the average error
                  DivB_Ave += SQR( DivB );   // L2 norm
                  NAve ++;
               } // i,j,k
            } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         } // if ( MPI_Rank == TargetRank )

//       do one rank at a time only in the verbose mode
         if ( Verbose )
         {
            MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );
            MPI_Barrier( MPI_COMM_WORLD );
         }
      } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks
   double DivB_Ave_AllRank, DivB_Max_AllRank;
   long   NAve_AllRank, NFailed_AllRank;

   MPI_Reduce( &NAve,    &NAve_AllRank,      1, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &NFailed, &NFailed_AllRank,   1, MPI_LONG,   MPI_SUM, 0, MPI_COMM_WORLD );

   MPI_Reduce( &DivB_Ave, &DivB_Ave_AllRank, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &DivB_Max, &DivB_Max_AllRank, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )    DivB_Ave_AllRank = sqrt( DivB_Ave_AllRank / NAve_AllRank );


// record results in a file
   if ( MPI_Rank == 0 )
   {
      static bool FirstTime = true;
      char FileName[2*MAX_STRING];
      sprintf( FileName, "%s/Record__DivB", OUTPUT_DIR );

//    output header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         else
         {
            FILE *File = fopen( FileName, "w" );
            Aux_Message( File, "\n# Tolerated Error: %13.7e\n\n", ERROR_TOLERANCE );
            Aux_Message( File, "#-------------------------------------------------------------------------\n" );
            Aux_Message( File, "#%12s   %10s   %13s   %13s   %13s\n", "Time", "Step", "AveError", "MaxError", "FailedCells" );
            fclose( File );
         }

         FirstTime = false;
      }

      FILE *File = fopen( FileName, "a" );
      Aux_Message( File, "%13.7e   %10ld   %13.7e   %13.7e   %13ld\n", Time[0], Step, DivB_Ave_AllRank, DivB_Max_AllRank, NFailed_AllRank );
      fclose( File );
   } // if ( MPI_Rank == 0 )


   if ( Pass  &&  Verbose  &&  MPI_Rank == 0 )
      Aux_Message( stdout, "\"%s\" : <%s> PASSED, Time = %13.7e, Step = %ld\n",
                   comment, __FUNCTION__, Time[0], Step );

} // FUNCTION : MHD_Aux_Check_DivergenceB



#endif // #ifdef MHD
