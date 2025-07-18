#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_DumpTable
// Description :  Load the dump table from the file "Input__DumpTable"
//-------------------------------------------------------------------------------------------------------
void Init_Load_DumpTable()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_DumpTable ...\n" );


   const char FileName[] = "Input__DumpTable";

   if ( !Aux_CheckFileExist(FileName) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   FILE *File = fopen( FileName, "r" );

   const int MaxLine = 10000;
   char *input_line = NULL;
   size_t len = 0;
   int Trash, line, n;


// allocate the dump table
   DumpTable = new double [MaxLine];


// skip the header
   getline( &input_line, &len, File );

// begin to read
   for (line=0; line<MaxLine; line++)
   {
      n = getline( &input_line, &len, File );

//    check
      if ( n <= 1 )
         Aux_Error( ERROR_INFO, "incorrect reading at line %d of the file <%s> !!\n", line+2, FileName );

      sscanf( input_line, "%d%lf", &Trash, &DumpTable[line] );

//    stop the reading
      if ( input_line[0] == 42 )                   // '*' == 42
      {
//       ensure that at least one data dump time is loaded
         if ( line == 0 )
            Aux_Error( ERROR_INFO, "please provide at least one data dump time in the dump table !!\n" );

         DumpTable_NDump = line;          // record the number of data dumps
         DumpTable[line] = __FLT_MAX__;   // set the next dump as an extremely large number

//       for END_T < 0.0, reset the ending time as the time of the last dump
         if ( END_T < 0.0 )
         {
            END_T = DumpTable[line-1];

            if ( MPI_Rank == 0 )
               Aux_Message( stdout, "NOTE : END_T is reset to the time of the last data dump in %s: %13.7e\n",
                            FileName, END_T );
         }

//       for END_T >= 0.0, reset the ending time to the largest data dump time not greater than the input END_T
         else
         {
            for (int t=line-1; t>=0; t--)
            {
               if ( END_T >= DumpTable[t] )
               {
                  END_T = DumpTable[t];

                  if ( MPI_Rank == 0 )
                     Aux_Message( stdout, "NOTE : END_T is reset to the largest time in %s not greater "
                                          "than the input END_T: %13.7e\n",
                                  FileName, END_T );
                  break;
               }
            }
         } // if ( END_T < 0.0 ) ... else ...


//       verify the loaded dump table
         for (int t=1; t<=line; t++)
         {
            if ( DumpTable[t] <= DumpTable[t-1] )
               Aux_Error( ERROR_INFO, "values recorded in \"%s\" must be monotonically increasing !!\n",
                          FileName );
         }

         break;

      } // if ( input_line[0] == 42 )
   } // for (line=0; line<MaxLine; line++)


   if ( line == MaxLine )
      Aux_Error( ERROR_INFO, "please prepare a * symbol at the end of the file <%s> !!\n", FileName );


   fclose( File );

   if ( input_line != NULL )     free( input_line );


// insert initial time Time[0]
   for (int t=0; t<DumpTable_NDump; t++)
   {
//    check if the initial time is already in the dump table
      if (   (  Time[0] != 0.0 && fabs( (Time[0]-DumpTable[t])/Time[0] ) < 1.0e-8   )  ||
             (  Time[0] == 0.0 && fabs(  Time[0]-DumpTable[t]          ) < 1.0e-12  )      )    break;

//    insert the initial time into the dump table
      if (  ( Time[0] < DumpTable[t] )  &&  ( t == 0 || Time[0] > DumpTable[t-1] )  )
      {
         for (int s=DumpTable_NDump; s>t; s--)  DumpTable[s] = DumpTable[s-1];
         DumpTable[t] = Time[0];
         DumpTable_NDump ++;

         if ( MPI_Rank == 0 )
            Aux_Message( stdout, "NOTE : initial time %13.7e is inserted into the dump table\n", Time[0] );

         break;
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_Load_DumpTable ... done\n" );

} // FUNCTION : Init_Load_DumpTable
