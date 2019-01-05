#include "GAMER.h"

#define COMMENT_SYM     '#'      // comment symbol
#define DELIMITER       " \t"    // delimiter characters used by strtok()




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_LoadTable
// Description :  Load the target columns from the table
//
// Note        :  1. Overloaded with different types
//                2. Put the target columns in "TCol[]", which must be sorted into ascending numerical order
//                   in advance
//                3. Allocate memory for the pointer "Data" if AllocMem == true
//                   --> Must be freed manually
//                4. Delimiter characters for strtok() are defined by DELIMITER
//
// Parameter   :  Data        : Pointer to be allocated (if AllocMem == true) and to store the data
//                              --> call-by-reference
//                FileName    : Filename of the target table
//                NCol_Target : Total number of target columns
//                TCol        : Target columns (must be sorted into ascending numerical order in advance)
//                RowMajor    : true/false --> store data into "Data" in the row-/column-major order
//                              -->    Row-major: Data[Row][Column]
//                                  Column-major: Data[Column][Row]
//                AllocMem    : true/false --> allocate/do not allocate memory for the pointer "Data"
//
// Return      :  Total number of matched rows
//-------------------------------------------------------------------------------------------------------
template <typename T>
int Aux_LoadTable( T *&Data, const char *FileName, const int NCol_Target, const int TCol[], const bool RowMajor,
                   const bool AllocMem )
{

// TCol[] must be sorted into ascending numerical order in advance
   for (int t=1; t<NCol_Target; t++)
   {
      if ( TCol[t] <= TCol[t-1] )
         Aux_Error( ERROR_INFO, "TCol is not in ascending numerical order ([%d]=%d, [%d]=%d) !!\n",
                    t-1, TCol[t-1], t, TCol[t] );
   }


// count the number of rows
   const int NRow_Target = Aux_CountRow( FileName );


// allocate memory
   if ( AllocMem )   Data = new T [NCol_Target*NRow_Target];


// load data
   char  FirstItem[MAX_STRING];
   int   NCol_Check, NthCol, NCol_Match, NthRow=0, NRow_Match=0;
   char *Token=NULL;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NCol_Check = sscanf( Line, "%s", FirstItem );
      NthRow ++;

//    skip empty lines or lines starting with the comment symbol
//    --> must check NCheck < 0 as well since EOF is negative
      if ( NCol_Check <= 0  ||  FirstItem[0] == COMMENT_SYM )    continue;

//    loop over all tokens seperated by the delimiter characters
      NthCol     = 0;
      NCol_Match = 0;
      while (  NthCol <= TCol[NCol_Target-1]  &&  ( Token = (NthCol==0)?strtok(Line,DELIMITER):strtok(NULL,DELIMITER) ) != NULL  )
      {
         if ( NthCol == TCol[NCol_Match] )
         {
            if ( RowMajor )   Data[ NRow_Match*NCol_Target + NCol_Match ] = atof( Token );
            else              Data[ NCol_Match*NRow_Target + NRow_Match ] = atof( Token );

            NCol_Match ++;
         }

         NthCol ++;
      }

//    check if we find all target columns
      if ( NCol_Match != NCol_Target )
         Aux_Error( ERROR_INFO, "Number of matched columns (%d) != expect (%d) at row %d !!\n",
                    NCol_Match, NCol_Target, NthRow-1 );

      NRow_Match ++;
   } // while ( fgets(Line, MAX_STRING, File) != NULL )

// check if we find all target rows
   if ( NRow_Match != NRow_Target )
      Aux_Error( ERROR_INFO, "Number of matched rows (%d) != expect (%d) !!\n", NRow_Match, NRow_Target );

   fclose( File );
   delete [] Line;


// return the total number of matched rows
   return NRow_Match;

} // FUNCTION : Aux_LoadTable



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CountRow
// Description :  Count the total number of data rows in the target file
//
// Note        :  1. Empty lines and lines starting with the comment symbol will be skipped
//                   --> The comment symbol is defined by COMMENT_SYM
//
// Parameter   :  FileName : Filename of the target table
//
// Return      :  Total number of matched rows
//-------------------------------------------------------------------------------------------------------
int Aux_CountRow( const char *FileName )
{

   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "table \"%s\" does not exist !!\n", FileName );


   char FirstItem[MAX_STRING];
   int  NRow=0, NItem;

   char *Line = new char [MAX_STRING];
   FILE *File = fopen( FileName, "r" );

   while ( fgets(Line, MAX_STRING, File) != NULL )
   {
      NItem = sscanf( Line, "%s", FirstItem );

//    skip empty lines and lines starting with the comment symbol
      if ( NItem > 0  &&  FirstItem[0] != COMMENT_SYM )  NRow ++;
   }

   fclose( File );
   delete [] Line;


   return NRow;

} // FUNCTION : Aux_CountRow



// explicit template instantiation
template int Aux_LoadTable <float > ( float  *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <double> ( double *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <int   > ( int    *&, const char *, const int, const int [], const bool, const bool );
template int Aux_LoadTable <long  > ( long   *&, const char *, const int, const int [], const bool, const bool );

