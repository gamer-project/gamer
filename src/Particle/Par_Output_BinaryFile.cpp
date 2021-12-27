#include "GAMER.h"

#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_BinaryFile
// Description :  Output the particle position and velocity
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_Output_BinaryFile( const char *FileName )
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s (DumpID = %d) ...\n", __FUNCTION__, DumpID );
//   int N_att = PAR_NATT_TOTAL;
   FILE *File;

// check
   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
   {
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );
      File = fopen( FileName, "w" );
      fclose(File);
   }
   MPI_Barrier( MPI_COMM_WORLD );

   real *attribute_buff = (real*)malloc(sizeof(real)*amr->Par->NPar_AcPlusInac);
// data
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "a" );
            long counter = 0;
            for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
            {
//             skip inactive particles
               if ( amr->Par->Mass[p] < 0.0 )
                  continue;
               else
               {
                  attribute_buff[counter] = amr->Par->Attribute[v][p];
                  counter ++;
               }
            }
////          record number of activate particle and number of attribute
//            if (v==0)
//            {
//                MPI_Gather(&counter, 1, MPI_LONG, &counter_gather, 1, MPI_LONG, 0, MPI_COMM_WORLD);
//                counter_gather += counter;
//                if (MPI_Rank==0) 
//                {
//                    fwrite(&counter_gather, sizeof(long), 1, File);
//                    printf("%ld\n", counter_gather);
//                    fwrite(&N_att, sizeof(int), 1, File);
//                }
//            }
////
//          record particle data
            fwrite(attribute_buff, sizeof(real), counter, File);
            fclose( File );
         } // if ( MPI_Rank == TargetMPIRank )
         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   } // for (int v=0; v<PAR_NATT_TOTAL; v++)

   free(attribute_buff);
} // FUNCTION : Par_Output_BinaryFile
#endif // end of ifdef PARTICLE
