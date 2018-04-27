#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_ExchangeDataPatchList
// Description :  Output SendP_IDList[lv] or RecvP_IDList[lv]
//
// Option      :  0 -> RecvP_IDList[lv]
//             :  1 -> SendP_IDList[lv]
//-------------------------------------------------------------------------------------------------------
void Output_ExchangeDataPatchList( const int option, const int lv, const char *comment )
{

// check
   if ( option != 0  &&  option != 1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "option", option );


   char FileName[100];
   if ( option )  sprintf( FileName, "SendDataPatchList_%d_%d", MPI_Rank, lv );
   else           sprintf( FileName, "RecvDataPatchList_%d_%d", MPI_Rank, lv );

   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Time = %13.7e  Step = %ld  Rank = %d  Level = %d\n\n", Time[0], Step, MPI_Rank, lv );


   for (int s=0; s<26; s++)
   {
      const int NP    = ( option ) ? amr->ParaVar->SendP_NList [lv][s] : amr->ParaVar->RecvP_NList [lv][s];
      const int *List = ( option ) ? amr->ParaVar->SendP_IDList[lv][s] : amr->ParaVar->RecvP_IDList[lv][s];

      fprintf( File, "Face = %d     Length = %d\n", s, NP );

      for (int P=0; P<NP; P++)   fprintf( File, "%5d ", List[P] );

      fprintf( File, "\n\n" );
   }

   fclose( File );

} // FUNCTION : Output_ExchangeDataPatchList



#endif // #ifndef SERIAL
