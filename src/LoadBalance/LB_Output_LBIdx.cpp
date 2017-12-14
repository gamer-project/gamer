#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Output_LBIdx
// Description :  Output the load-balance indices and their corresponding coordinates for all real patches
//
// Note        :  Output file can be directly plotted by gnuplot
//                For example, try "splot 'LBIdxMap_Lv00' u 2:3:4 every :::0::0 w lp"
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void LB_Output_LBIdx( const int lv )
{

   const int NP = amr->NPatchComma[lv][1];
   long LB_Idx, *List;
   int  Coord[3];
   char FileName[100];

// set the filename   
   sprintf( FileName, "LBIdxMap_Lv%02d", lv );


// check if the target file already exists
   if ( Aux_CheckFileExist(FileName)  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


// output the load-balance indices   
   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         FILE *File = fopen( FileName, "a" );

         if ( YourTurn == 0 )    fprintf( File, "%8s: %8s %8s %8s\n", "LB_Idx", "x", "y", "z" );

//       get LB_Idx
         List = new long [NP];
         for (int PID=0; PID<NP; PID++)
         {
            LB_Idx    = amr->patch[0][lv][PID]->LB_Idx;
            List[PID] = LB_Idx;
         }

//       sorting         
         Mis_Heapsort( NP, List, NULL );

//       get corner         
         for (int t=0; t<NP; t++)
         {
            LB_Index2Corner( lv, List[t], Coord, CHECK_ON );
            fprintf( File, "%8ld: %8d %8d %8d\n", List[t], Coord[0], Coord[1], Coord[2] );
         }

//       add an empty line between the outputs of different MPI ranks         
         fprintf( File, "\n" );

         delete [] List;
         fclose( File );

      } // if ( MPI_Rank == YourTurn )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)

} // FUNCTION : LB_Output_LBIdx



#endif // #ifdef LOAD_BALANCE
