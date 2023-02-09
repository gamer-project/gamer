#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Load_ExecuteTable
// Description :  Load the yt inline analysis execuation table from the file "Input__ExecuteYTTable"
//-------------------------------------------------------------------------------------------------------
void YT_Load_ExecuteTable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "YT_Load_ExecuteTable ...\n" );


   const char FileName[] = "Input__ExecuteYTTable";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   const int MaxLine = 10000;
   char *input_line = NULL;
   size_t len = 0;
   int Trash, line, n;


// allocate the dump table
   ExecuteYTTable = new double [MaxLine];


// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (line=0; line<MaxLine; line++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
         Aux_Error( ERROR_INFO, "incorrect reading at line %d of the file <%s> !!\n", line+2, FileName );

      sscanf( input_line, "%d%lf", &Trash, &ExecuteYTTable[line] );

//    stop the reading
      if ( input_line[0] == 42 )                   // '*' == 42
      {

//       ensure that at least one data dump time is loaded
         if ( line == 0 )
            Aux_Error( ERROR_INFO, "please provide at least one data dump time in the dump table !!\n" );

         ExecuteYTTable_NExecute   = line;                 // record the number of yt inline analysis execution
         ExecuteYTTable[line]      = __FLT_MAX__;          // set the next execution as an extremely large number

//    do not allow the END_T substituted by the last yt inline analysis execution time
//         if ( ExecuteYTTable[line-1] < END_T )
//         {
//            END_T          = ExecuteYTTable[line-1];    // reset the ending time as the time of the last dump
//
//            if ( MPI_Rank == 0 )
//               Aux_Message( stdout, "NOTE : the END_T is reset to the time of the last yt inline analysis execution = %13.7e\n",
//                            END_T );
//         }


//       verify the loaded yt inline analysis execution table
         for (int t=1; t<=line; t++)
         {
            if ( ExecuteYTTable[t] < ExecuteYTTable[t-1] )
               Aux_Error( ERROR_INFO, "values recorded in \"%s\" must be monotonically increasing !!\n",
                          FileName );
         }

         break;

      } // if ( input_line[0] == 42 )
   } // for (line=0; line<MaxLine; line++)


   if ( line == MaxLine )
      Aux_Error( ERROR_INFO, "please prepare a symbol * in the end of the file <%s> !!\n", FileName );


   fclose( File );

   if ( input_line != NULL )     free( input_line );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "YT_Load_ExecuteTable ... done\n" );

} // FUNCTION : YT_Load_ExecuteTable
