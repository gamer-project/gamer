#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <unistd.h>

using namespace std;

void ReadOption( int argc, char **argv );
void Init();
void Load();
int CountRow( const char *FileName, const int NHeader );
bool CheckEmptyString( const char *InputString );
void End();
void CheckMonotonicity();


const int NHeader_Profile = 1;   // number of header lines in FileName_Profile

char *FileName_Profile = NULL;

double *r       = NULL;          // radius
double *MaxDens = NULL;          // peak density
double *AveDens = NULL;          // average density

int NRow_Profile;                // number of data rows in FileName_Profile
int NCheck = 20;                 // number of data rows to be checked




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command-line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while ( (c = getopt(argc, argv, "hp:n:")) != -1 )
   {
      switch ( c )
      {
         case 'p': FileName_Profile = optarg;
                   break;
         case 'n': NCheck           = atoi(optarg);
                   break;
         case 'h':
         case '?': cerr << endl << "usage: " << argv[0]
                        << " [-h (for help)] [-p profile file] [-n number of data rows to be checked [20]]"
                        << endl << endl;
                   exit( 1 );

      } // switch ( c ) ...
   } // while ...


// check the input command-line options
   if ( FileName_Profile == NULL )
   {
      fprintf( stderr, "ERROR : please provide the filename of the profile file [-p profile file] !!\n" );
      exit( 1 );
   }

   if ( NCheck < 0 )
   {
      fprintf( stderr, "ERROR : NCheck = %d < 0 !!\n", NCheck );
      exit( 1 );
   }


// record the command-line options
   fprintf( stdout, "\n" );
   fprintf( stdout, "-------------------------------------------------------------------------------------------------\n" );
   fprintf( stdout, "Command-line arguments :\n" );
   for (int v=0; v<argc; v++)    fprintf( stdout, " %s", argv[v] );
   fprintf( stdout, "\n" );
   fprintf( stdout, "-------------------------------------------------------------------------------------------------\n" );
   fprintf( stdout, "\n" );

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  Init
// Description :  Initialization
//-------------------------------------------------------------------------------------------------------
void Init()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   NRow_Profile = CountRow( FileName_Profile, NHeader_Profile );

   r       = new double [NRow_Profile];
   MaxDens = new double [NRow_Profile];
   AveDens = new double [NRow_Profile];


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : Init



//-------------------------------------------------------------------------------------------------------
// Function    :  End
// Description :  End the program
//-------------------------------------------------------------------------------------------------------
void End()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   delete [] r;
   delete [] MaxDens;
   delete [] AveDens;


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : End



//-------------------------------------------------------------------------------------------------------
// Function    :  Load
// Description :  Load the density profile
//-------------------------------------------------------------------------------------------------------
void Load()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   const int MaxLine = 1024;        // maximum number of characters per line
   int    NRow;
   double tmp;
   char *Line = new char[MaxLine];


// load the density profile
   FILE *File_Profile = fopen( FileName_Profile, "r" );

   if ( File_Profile == NULL )
   {
      fprintf( stderr, "ERROR : input profile file \"%s\" does not exist !!\n", FileName_Profile );
      exit( 1 );
   }

   NRow = 0;

   for (int t=0; t<NHeader_Profile; t++)  fgets( Line, MaxLine, File_Profile );  // skip headers

   while ( fgets( Line, MaxLine, File_Profile ) != NULL )
   {
      if (  !CheckEmptyString( Line )  )
      {
         sscanf( Line, "%lf%lf%lf%lf%lf%lf%lf%lf", r+NRow, &tmp, AveDens+NRow, &tmp, MaxDens+NRow, &tmp, &tmp, &tmp );
         NRow ++;
      }
   }

   if ( NRow != NRow_Profile )
   {
      fprintf( stderr, "ERROR : NRow (%d) != NRow_Profile (%d) !!\n", NRow, NRow_Profile );
      exit( 1 );
   }

   fclose( File_Profile );

   delete [] Line;


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : Load



//-------------------------------------------------------------------------------------------------------
// Function    :  CountRow
// Description :  Count the number of rows in the input file
//
// Parameter   :  FileName : Name of the input file
//                NHeader  : Number of header lines to be skipped
//
// Return      :  NRow
//-------------------------------------------------------------------------------------------------------
int CountRow( const char *FileName, const int NHeader )
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   const int MaxLine = 1024;        // maximum number of characters per line
   int   NRow   = 0;                // number of data rows
   int   NEmpty = 0;                // number of empty rows
   char *Line = new char [MaxLine];

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist !!\n", FileName );
      exit( 1 );
   }

// get the number of data rows
   for (int t=0; t<NHeader; t++)  fgets( Line, MaxLine, File );    // skip headers

   while ( fgets( Line, MaxLine, File ) != NULL )
   {
      if (  !CheckEmptyString( Line )  )  NRow   ++;
      else                                NEmpty ++;
   }

   if ( NRow < 1 )
   {
      fprintf( stderr, "ERROR : no data rows are found !!\n" );
      exit( 1 );
   }

   fclose( File );
   delete [] Line;


   fprintf( stdout, "   Input file  = %s\n", FileName );
   fprintf( stdout, "   Header rows = %d\n", NHeader  );
   fprintf( stdout, "   Data   rows = %d\n", NRow     );
   fprintf( stdout, "   Empty  rows = %d\n", NEmpty   );

   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

   return NRow;

} // FUNCTION : CountRow



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckEmptyString
// Description :  Check whether the input string is empty
//
// Note        :  Empty string is defined as a string containing only " ", "\n" and "\t"
//
// Return      :  true  : The input string is empty
//                false : The input string is NOT empty
//-------------------------------------------------------------------------------------------------------
bool CheckEmptyString( const char *InputString )
{
   static const char *EmptyChar = " \n\t";

   return strspn( InputString, EmptyChar ) == strlen( InputString );

} // FUNCTION : CheckEmptyString



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckMonotonicity
// Description :  Check whether the peak and average density in the input profile is monotonic
//
// Note        :  Only check the first NCheck data
//
//-------------------------------------------------------------------------------------------------------
void CheckMonotonicity()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   if ( NCheck > NRow_Profile )
   {
      fprintf( stderr, "ERROR : NCheck (%d) > NRow_Profile (%d) !!\n", NCheck, NRow_Profile );
      exit( 1 );
   }

   for (int t=1; t<NCheck; t++)
   {
      if ( MaxDens[t] > MaxDens[t-1] )
         fprintf( stderr, "   WARNING : Row %2d -> %2d,    peak density %13.7e -> %13.7e\n",
                  t-1, t, MaxDens[t-1], MaxDens[t] );

      if ( AveDens[t] > AveDens[t-1] )
         fprintf( stderr, "   WARNING : Row %2d -> %2d, average density %13.7e -> %13.7e\n",
                  t-1, t, AveDens[t-1], AveDens[t] );
   }


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : CheckMonotonicity



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   ReadOption( argc, argv );

   Init();

   Load();

   CheckMonotonicity();

   End();

   return EXIT_SUCCESS;

} // FUNCTION : main


