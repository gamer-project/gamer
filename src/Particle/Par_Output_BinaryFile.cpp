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


// check
   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File;
   real *attribute_buff = (real*)malloc(sizeof(real)*amr->Par->NPar_AcPlusInac);
// data
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
          if ( MPI_Rank==0 )
              File = fopen( FileName, "w" );
          else
              File = fopen( FileName, "a" );
         for (int v=0; v<PAR_NATT_TOTAL; v++)
         {
            int counter = 0;
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
            fwrite(attribute_buff, sizeof(real), counter, File);

         }
         fclose( File );
         MPI_Barrier( MPI_COMM_WORLD );
      } // if ( MPI_Rank == TargetMPIRank )

   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   free(attribute_buff);
} // FUNCTION : Par_Output_BinaryFile
#endif // end of ifdef PARTICLE
