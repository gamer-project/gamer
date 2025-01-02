#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_BinaryFile
// Description :  Output particle attributes in the binary format
//
// Note        :  1. Number of data elements: N_attribute*N_activate_particle
//                2. Data of all active particles with a selected attribute will be dumped consecutively,
//                   followed by the next attribute, so on and so forth
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Output_BinaryFile( const char *FileName )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );


   FILE *File;

// check
   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
   {
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );
      File = fopen( FileName, "wb" );
      fclose(File);
   }
   MPI_Barrier( MPI_COMM_WORLD );


// dump floating-point data
   real_par *attribute_flt_buff = (real_par*)malloc( sizeof(real_par)*amr->Par->NPar_AcPlusInac );

   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   {
      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "ab" );
            long counter = 0;

//          store particle data in a buffer
            for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
            {
//             skip inactive particles
               if ( amr->Par->Mass[p] < 0.0 )   continue;
               else                             attribute_flt_buff[ counter ++ ] = amr->Par->AttributeFlt[v][p];
            }

//          dump data from the buffer
            fwrite( attribute_flt_buff, sizeof(real_par), counter, File );
            fclose( File );
         } // if ( MPI_Rank == TargetMPIRank )

         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   } // for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)

   free( attribute_flt_buff );


// dump integer data
   long_par *attribute_int_buff = (long_par*)malloc( sizeof(long_par)*amr->Par->NPar_AcPlusInac );

   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)
   {
      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "ab" );
            long counter = 0;

//          store particle data in a buffer
            for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
            {
//             skip inactive particles
               if ( amr->Par->Mass[p] < 0.0 )   continue;
               else                             attribute_int_buff[ counter ++ ] = amr->Par->AttributeInt[v][p];
            }

//          dump data from the buffer
            fwrite( attribute_int_buff, sizeof(long_par), counter, File );
            fclose( File );
         } // if ( MPI_Rank == TargetMPIRank )

         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   } // for (int v=0; v<PAR_NATT_INT_TOTAL; v++)

   free( attribute_int_buff );

} // FUNCTION : Par_Output_BinaryFile



#endif // #ifdef PARTICLE
