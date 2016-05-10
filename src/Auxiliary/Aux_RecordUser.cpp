#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_RecordUser
// Description :  Record user-specified information
//
// Note        :  1. Please turn on the option "OPT__RECORD_USER" 
//                2. This function will be called both during the program initialization and after each full update
// 
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Aux_RecordUser( )
{

   const char FileName[] = "Record__User";
   static bool FirstTime = true;

   if ( FirstTime )
   {
//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "%14s%14s%3s%14s\n",  "Time", "Step", "", "dt" );
         fclose( File_User );
      }

      FirstTime = false;
   }

// user-specified info
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e%14ld%3s%14.7e\n", Time[0], Step, "", dTime_Base );
      fclose( File_User );
   }

} // FUNCTION : Aux_RecordUser


