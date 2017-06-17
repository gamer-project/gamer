#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined PARTICLE )

static bool CheckEmptyString( const char *InputString );
static int  CountRow( const char *FileName );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadProfile_Merger
// Description :  Load the radial profiles of radius, gas mass density, and pressure
//
// Note        :  1. Regard "#" as comment
//                2. Assumed format (starting from the 0th column)
//                      radius   : 7th column
//                      density  : 2th column
//                      pressure : 6th column
//                3. Assumed input units: cgs
//
// Parameter   :  FileName : Name of the target profile table
//                Profile  : Pointer to be allocated to store the profile data
//                           --> Must be freed manually
//                NBin     : Number of radial bins
//
// Return      :  Profile, NBin
//-------------------------------------------------------------------------------------------------------
void LoadProfile_Merger( const char *FileName, double **Profile, int &NBin )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      %s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char  *Line       = new char [MaxLine];
   char  *FirstChar  = NULL;
   FILE  *File       = NULL;
   double tmp;
   int    NLoad;


// check the input file
   if (  ( File = fopen(FileName, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );


// get the total number of mass bins
   NBin = CountRow( FileName );


// allocate data
   for (int v=0; v<3; v++)    Profile[v] = new double [NBin];


// loop over all rows in the input file
   NLoad = 0;
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
//    skip empty lines
      if (  !CheckEmptyString( Line )  )
      {
         FirstChar = Line;

//       find the first non-empty character
         while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//       skip lines starting with "#"
         if ( *FirstChar != '#' )
         {
            if ( NLoad >= NBin )    Aux_Error( ERROR_INFO, "NLoad (%d) >= NBin (%d) !!\n", NLoad, NBin );

            sscanf( Line, "%lf%lf%lf%lf%lf%lf%lf%lf",
                    &tmp, &tmp, Profile[1]+NLoad, &tmp, &tmp, &tmp, Profile[2]+NLoad, Profile[0]+NLoad );

            NLoad ++;
         }
      }
   } // while ( fgets( Line, MaxLine, File ) != NULL )

   fclose( File );
   delete [] Line;


// verify the number of loaded data
   if ( NLoad != NBin )    Aux_Error( ERROR_INFO, "total number of loaded data (%d) != expect (%d) !!\n", NLoad, NBin );


// convert to code units (assuming the input units are cgs)
   for (int b=0; b<NBin; b++)
   {
      Profile[0][b] /= UNIT_L;
      Profile[1][b] /= UNIT_D;
      Profile[2][b] /= UNIT_P;
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "      %s ... done\n", __FUNCTION__ );

} // FUNCTION : LoadProfile_Merger



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckEmptyString
// Description :  Check whether the input string is empty
//
// Note        :  Empty string is defined as a string containing only " ", "\n" and "\t"
//
// Return      :  true  : Input string is empty
//                false : Input string is NOT empty
//-------------------------------------------------------------------------------------------------------
bool CheckEmptyString( const char *InputString )
{
   static const char *EmptyChar = " \n\t";

   return strspn( InputString, EmptyChar ) == strlen( InputString );

} // FUNCTION : CheckEmptyString



//-------------------------------------------------------------------------------------------------------
// Function    :  CountRow
// Description :  Count the total number of rows in the input table
//
// Parameter   :  FileName : Name of the input file
//
// Return      :  NRow
//-------------------------------------------------------------------------------------------------------
int CountRow( const char *FileName )
{

   const int MaxLine = 1024;        // maximum number of characters per line

   int   NRow      = 0;             // number of data rows
   int   NEmpty    = 0;             // number of empty rows
   char *FirstChar = NULL;
   char *Line = new char [MaxLine];

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );

// get the number of data rows
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
      if (  !CheckEmptyString( Line )  )
      {
         FirstChar = Line;

//       find the first non-empty character
         while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//       skip lines starting with "#"
         if ( *FirstChar != '#' )   NRow   ++;
         else                       NEmpty ++;
      }

      else                          NEmpty ++;
   }

   if ( NRow < 1 )   Aux_Error( ERROR_INFO, "no data rows are found !!\n" );

   fclose( File );
   delete [] Line;


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "         Input file = %s\n", FileName );
      Aux_Message( stdout, "         Data  rows = %d\n", NRow     );
      Aux_Message( stdout, "         Empty rows = %d\n", NEmpty   );
   }

   return NRow;

} // FUNCTION : CountRow



#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )
