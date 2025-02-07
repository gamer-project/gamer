#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_ExchangeFluxPatchList
// Description :  Output SendF_IDList[lv] or RecvF_IDList[lv]
//
// Option      :  0 -> SendF_IDList[lv]
//             :  1 -> RecvF_IDList[lv]
//-------------------------------------------------------------------------------------------------------
void Output_ExchangeFluxPatchList( const int option, const int lv, const char *comment )
{

   if ( option < 0  ||  option > 1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "option", option );

   if ( lv < 0  ||  lv >= NLEVEL-1 )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );


   char FileName[2*MAX_STRING];
   switch ( option )
   {
      case 0:  sprintf( FileName, "%s/SendFluxPatchList_%d_%d", OUTPUT_DIR, MPI_Rank, lv );
               break;
      case 1:  sprintf( FileName, "%s/RecvFluxPatchList_%d_%d", OUTPUT_DIR, MPI_Rank, lv );
               break;
   }

   if ( comment != NULL )
   {
      strcat( FileName, "_" );
      strcat( FileName, comment );
   }

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );


   FILE *File = fopen( FileName, "w" );

   fprintf( File, "Time = %13.7e  Step = %ld  Rank = %d  Level = %d\n\n", Time[0], Step, MPI_Rank, lv );


   int NP    = 0;
   int *List = NULL;

   for (int s=0; s<6; s++)
   {
      switch ( option )
      {
         case 1:  NP    = amr->ParaVar->SendF_NList [lv][s];
                  List  = amr->ParaVar->SendF_IDList[lv][s];
                  break;

         case 2:  NP    = amr->ParaVar->RecvF_NList [lv][s];
                  List  = amr->ParaVar->RecvF_IDList[lv][s];
                  break;
      }

      fprintf( File, "Face = %d     Length = %d\n", s, NP );

      for (int P=0; P<NP; P++)   fprintf( File, "%5d ", List[P] );

      fprintf( File, "\n\n" );
   }

   fclose( File );

} // FUNCTION : Output_ExchangeFluxPatchList



#endif // #ifndef SERIAL
