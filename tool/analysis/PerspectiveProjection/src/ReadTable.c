#include "../include/General.h"


#define COMMENT_SYMBOL '#'
#define MAX_STRING     100
#define DELIMITER      " "



//-------------------------------------------------------------------------------------------------------
// Function    :  ConntNumLine
// Description :  Count the number of rows in file
//                --> excluding the line beginning with comment symbol (#)
//                --> excluding the line composed with comment symbols and white spaces
// Note        :
// Parameter   :  table : File stream
// Return      :  The number of rows in table
//-------------------------------------------------------------------------------------------------------
int CountNumLine( FILE *table )
{
   uint64_t numLine = 0;
   char line[MAX_STRING];

   rewind(table);

// Get the number of charactors in each line
   while( fgets(line, MAX_STRING, table) )
   {
//   Skip the line beginning with comment symbol
     if ( line[0] == COMMENT_SYMBOL  ||  line[0] == '\n' )   continue;

//   Check the line beginning with white space
     if ( line[0] == ' ' )
     {
       bool skip = true;
       for (int i=0; line[i] != '\0'; i++)
       {
//       Keep the line when it has any meaningful characters
         if ( !(line[i] == COMMENT_SYMBOL  ||  line[i] == ' '  ||  line[i] == '\n') )
         {
            skip &= false;
            break;
         }

//       Break the loop immediately when encounter comment symbol
         if ( line[i] == COMMENT_SYMBOL )  break;
       } // for (int i=0; line[i] != '\0'; i++)

       if (skip)   continue;
     } // if ( line[0] == ' ' )

     numLine++;
   } // while( fgets(line, MAX_STRING, table) )

   return numLine;
} // FUNCTION : CountNumLine



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadTable
// Description :
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void ReadTable( int TargetColumn, int numRow, FILE *table, float *columnArray )
{
   rewind(table);

   char line[MAX_STRING];

   int r = 0;

   while( fgets(line, MAX_STRING, table) )
   {

//   Step 1: Skip the line beginning with comment symbol or newline
     if ( line[0] == COMMENT_SYMBOL || line[0] == '\n' ) continue;

//   Step 2: Check the line beginning with white space
     if ( line[0] == ' ' )
     {
       bool skip = true;

       for (int i=0; line[i] != '\0'; i++)
       {

//       Step 2-1: Kepp the line when it has any meaningful characters
         if ( !(line[i] == COMMENT_SYMBOL || line[i] == ' ' || line[i] == '\n') ){
            skip &= false;
            break;
         }

//       Step 2-2: Break the loop immediately when encounter comment symbol
         if ( line[i] == COMMENT_SYMBOL )  break;

       } // for (int i=0; line[i] != '\0'; i++)

       if (skip) continue;
     } // if ( line[0] == ' ' )

     // Step 3: parse a string into tokens
     char *token = strtok( line, DELIMITER );
     int c = 1;

     while( token != NULL )
     {
       if ( c++ == TargetColumn ) columnArray[r] = atof(token);
       token = strtok( NULL, DELIMITER );
     }

     if ( r == numRow-1 ) break;

     r++;
   } // while( fgets(line, MAX_STRING, table) )
} // FUNCTION : ReadTable
