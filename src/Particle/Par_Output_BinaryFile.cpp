#include "GAMER.h"

#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Output_BinaryFile
// Description :  Output the particle position and velocity
// Format      :  1. Total data numbers: N_attribute*N_activate_particle
//                2. Data of all activate particles with a selected attribute will be dumped consecutively, followed by the next attribute, so on and so forth
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

   real *attribute_buff = (real*)malloc(sizeof(real)*amr->Par->NPar_AcPlusInac);
// data
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
      for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
      {
         if ( MPI_Rank == TargetMPIRank )
         {
            File = fopen( FileName, "ab" );
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
