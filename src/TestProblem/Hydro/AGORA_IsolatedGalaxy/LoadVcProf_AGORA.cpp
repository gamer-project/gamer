#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined PARTICLE )

bool CheckEmptyString_AGORA( const char *InputString );
int  CountRow_AGORA( const char *Filename );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadVcProf_AGORA
// Description :  Load the circular velocity radial profile
//
// Note        :  1. Regard "#" as comment
//                2. Assumed format [radius in kpc] [velocity in km/s]
//                3. This function will allocate memory for the input pointer "Profile"
//                   --> One must deallocate it manually
//
// Parameter   :  Filename : Name of the target profile table
//                Profile  : Pointer to be allocated to store the profile data
//                NBin     : Number of radial bins
//
// Return      :  Profile, NBin
//-------------------------------------------------------------------------------------------------------
void LoadVcProf_AGORA( const char *Filename, double **Profile, int &NBin )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ...\n", __FUNCTION__ );


   const int MaxLine = 1024;                 // maximum number of characters per line
   char  *Line       = new char [MaxLine];
   char  *FirstChar  = NULL;
   FILE  *File       = NULL;
   int    NLoad;


// check the input file
   if (  ( File = fopen(Filename, "r") ) == NULL  )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", Filename );


// get the total number of radial bins
   NBin = CountRow_AGORA( Filename );


// allocate data
   for (int v=0; v<2; v++)    Profile[v] = new double [NBin];


// loop over all rows in the input file
   NLoad = 0;
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
//    skip empty lines
      if (  !CheckEmptyString_AGORA( Line )  )
      {
         FirstChar = Line;

//       find the first non-empty character
         while ( *FirstChar == ' '  ||  *FirstChar == '\t' )   FirstChar ++;

//       skip lines starting with "#"
         if ( *FirstChar != '#' )
         {
            if ( NLoad >= NBin )    Aux_Error( ERROR_INFO, "NLoad (%d) >= NBin (%d) !!\n", NLoad, NBin );

            sscanf( Line, "%lf%lf", Profile[0]+NLoad, Profile[1]+NLoad );

            NLoad ++;
         }
      }
   } // while ( fgets( Line, MaxLine, File ) != NULL )

   fclose( File );
   delete [] Line;


// verify the number of loaded data
   if ( NLoad != NBin )    Aux_Error( ERROR_INFO, "total number of loaded data (%d) != expect (%d) !!\n", NLoad, NBin );


// convert to code units
   for (int b=0; b<NBin; b++)
   {
      Profile[0][b] *= Const_kpc / UNIT_L;
      Profile[1][b] *= Const_km  / UNIT_V;
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   %s ... done\n", __FUNCTION__ );

} // FUNCTION : LoadVcProf_AGORA



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckEmptyString_AGORA
// Description :  Check whether the input string is empty
//
// Note        :  Empty string is defined as a string containing only " ", "\n" and "\t"
//
// Return      :  true  : Input string is empty
//                false : Input string is NOT empty
//-------------------------------------------------------------------------------------------------------
bool CheckEmptyString_AGORA( const char *InputString )
{
   static const char *EmptyChar = " \n\t";

   return strspn( InputString, EmptyChar ) == strlen( InputString );

} // FUNCTION : CheckEmptyString_AGORA



//-------------------------------------------------------------------------------------------------------
// Function    :  CountRow_AGORA
// Description :  Count the total number of rows in the input table
//
// Parameter   :  Filename : Name of the input file
//
// Return      :  NRow
//-------------------------------------------------------------------------------------------------------
int CountRow_AGORA( const char *Filename )
{

   const int MaxLine = 1024;        // maximum number of characters per line

   int   NRow      = 0;             // number of data rows
   int   NEmpty    = 0;             // number of empty rows
   char *FirstChar = NULL;
   char *Line = new char [MaxLine];

   FILE *File = fopen( Filename, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", Filename );

// get the number of data rows
   while ( fgets( Line, MaxLine, File ) != NULL )
   {
      if (  !CheckEmptyString_AGORA( Line )  )
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
      Aux_Message( stdout, "      Input file = %s\n", Filename );
      Aux_Message( stdout, "      Data  rows = %d\n", NRow     );
      Aux_Message( stdout, "      Empty rows = %d\n", NEmpty   );
   }

   return NRow;

} // FUNCTION : CountRow_AGORA



#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )
