#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_PatchCorner
// Description :  Output the corner coordinates of all real and buffer patches, which can be visualized
//                by "splot" in gnuplot
//
// Note        :  An empty line is inserted between the data of real and buffer patches
//                --> For gnuplot plot, one can use (1) "every :::0::0" to plot the real   patches
//                                                  (2) "every :::1::1" to plot the buffer patches
//
// Parameter   :  lv       : Target refinement level 
//                comment  : String to attach to the end of the file name
//-------------------------------------------------------------------------------------------------------
void Output_PatchCorner( const int lv, const char *comment )
{

   char FileName[100];
   sprintf( FileName, "PatchCorner_%05d_%02d", MPI_Rank, lv );
   if ( comment != NULL )       
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


   const int NReal = amr->NPatchComma[lv][1];
   const int NBuff = amr->NPatchComma[lv][27] - amr->NPatchComma[lv][1];

   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Time %13.7e  Step %ld  Counter %ld  Rank %d  Level %d  NRealPatch %d  NBufferPatch %d\n", 
            Time[lv], Step, AdvanceCounter[lv], MPI_Rank, lv, NReal, NBuff );
   fprintf( File, "=========================================================================================\n" );
   fprintf( File, "%8s   %10s   %10s   %10s\n", "PID", "Corner[x]", "Corner[y]", "Corner[z]" );

// output real patches
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      fprintf( File, "%8d   %10d   %10d   %10d\n", PID, amr->patch[0][lv][PID]->corner[0], 
                                                        amr->patch[0][lv][PID]->corner[1], 
                                                        amr->patch[0][lv][PID]->corner[2] );
   fprintf( File, "\n" );

// output buffer patches
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][27]; PID++)
      fprintf( File, "%8d   %10d   %10d   %10d\n", PID, amr->patch[0][lv][PID]->corner[0], 
                                                        amr->patch[0][lv][PID]->corner[1], 
                                                        amr->patch[0][lv][PID]->corner[2] );

   fclose( File );

} // FUNCTION : Output_PatchCorner


