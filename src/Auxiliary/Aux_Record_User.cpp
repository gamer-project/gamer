#include "GAMER.h"

// declare as static so that other functions cannot invoke it directly and must use the function pointer
static void Aux_Record_User_Template( );

// this function pointer must be set by a test problem initializer
void (*Aux_Record_User_Ptr)() = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_User_Template
// Description :  Template of recording user-specified information
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_User_Template()
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__User", OUTPUT_DIR );

   if ( FirstTime )
   {
//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%13s%14s%3s%14s\n",  "Time", "Step", "", "dt" );
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

} // FUNCTION : Aux_Record_User_Template


