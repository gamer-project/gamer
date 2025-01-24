#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_TextFile
// Description :  Output particle attributes in the text format
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Output_TextFile( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


// check
   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// header
   if ( MPI_Rank == 0 )
   {
      FILE *File = fopen( FileName, "w" );

      fprintf( File, "#Time %20.14e   Step %13ld   Active Particles %13ld\n\n",
               Time[0], Step, amr->Par->NPar_Active_AllRank );

      fprintf( File, "#" );

      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
      fprintf( File, " %*s", (v==0)?StrLen_Flt-1:StrLen_Flt, ParAttFltLabel[v] );

      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
      fprintf( File, " %*s", OPT__OUTPUT_TEXT_LENGTH_INT, ParAttIntLabel[v] );

      fprintf( File, "\n" );

      fclose( File );
   }


// data
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         FILE *File = fopen( FileName, "a" );

         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
//          skip inactive particles
            if ( amr->Par->Mass[p] < 0.0 )   continue;

            for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   fprintf( File, BlankPlusFormat_Flt, amr->Par->AttributeFlt[v][p] );
            for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   fprintf( File, " %*ld", OPT__OUTPUT_TEXT_LENGTH_INT, (long)amr->Par->AttributeInt[v][p] );

            fprintf( File, "\n" );
         }

         fclose( File );
      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ... done\n", __FUNCTION__, DumpID );

} // FUNCTION : Par_Output_TextFile



#endif // #ifdef PARTICLE
