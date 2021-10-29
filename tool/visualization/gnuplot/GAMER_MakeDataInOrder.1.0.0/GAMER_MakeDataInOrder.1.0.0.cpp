#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <unistd.h>
using namespace std;

// single/double precision
#ifdef FLOAT8
   typedef double real;
#else
   typedef float  real;
#endif

// extreme values
#ifndef __INT_MAX__
#  define __INT_MAX__     2147483647
#endif


void Heapsort( const int N, real Array[], int IdxTable[] );
void Heapsort_SiftDown( const int L, const int R, real Array[], int IdxTable[] );
void ReadOption( int argc, char **argv );
void SortKey();
void LoadData();
void End();

const char *delim     = " \n\t";             // delimiter for "strtok"
const int   MaxColumn = 20;                  // maximum number of tokens
const int   MaxLine   = 1024;                // maximum number of characters per line

int    Header         = 1;                   // number of header lines
int    TColumn        = 0;                   // the targeted column to be sorted (start with 0)
real  *Key            = NULL;                // data to be sorted
int   *Key_IdxTable   = NULL;                // index table of Key
int   *Key_RankTable  = NULL;                // rank table of Key
char (*Data)[MaxLine] = NULL;                // string array to store the input data line by line
char  *FileName_In    = NULL;
char  *FileName_Out   = NULL;

int   NRow, NColumn;                         // number of data rows and columns (excluding header)





//-------------------------------------------------------------------------------------------------------
// Function    :  DumpData
// Description :  Dump data
//-------------------------------------------------------------------------------------------------------
void DumpData()
{

   fprintf( stdout, "%s ... ", __FUNCTION__ );   fflush( stdout );


   FILE *File = fopen( FileName_Out, "a" );

   for (int Row=0; Row<NRow; Row++)    fprintf( File, "%s", Data[Row] );

   fclose( File );


   fprintf( stdout, "done\n" );   fflush( stdout );

} // FUNCTION : DumpData



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load data from the input file and store them in the array "Data" according to the
//                rank table "Key_RankTable"
//-------------------------------------------------------------------------------------------------------
void LoadData()
{

   fprintf( stdout, "%s ... ", __FUNCTION__ );   fflush( stdout );


   FILE *File = fopen( FileName_In, "r" );

   if ( File == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }


   char *Line = new char [MaxLine];

// skip header
   for (int t=0; t<Header; t++)  fgets( Line, MaxLine, File );

// load data
   for (int Row=0; Row<NRow; Row++)
   {
      fgets( Line, MaxLine, File );
      strcpy( Data[ Key_RankTable[Row] ], Line );
   }

   fclose( File );
   delete [] Line;


   fprintf( stdout, "done\n" );   fflush( stdout );

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  SortKey
// Description :  Sort the targeted column and get the rank table
//-------------------------------------------------------------------------------------------------------
void SortKey()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   int   Column;
   char *str  = NULL, *Temp = NULL;
   char *Line = new char [MaxLine];


// output header
   FILE *File_In = fopen( FileName_In, "r" );
   if ( File_In == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }

   FILE *File_Check = fopen( FileName_Out, "r" );
   if ( File_Check != NULL )
   {
      fprintf( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName_Out );
      fclose( File_Check );
   }

   FILE *File_Out = fopen( FileName_Out, "w" );
   for (int t=0; t<Header; t++)
   {
      fgets( Line, MaxLine, File_In );
      fprintf( File_Out, "%s", Line );
   }
   fclose( File_Out );


// get the number of data columns
   fgets( Line, MaxLine, File_In );
   NColumn = 0;

   for (Column=0, str=Line; Column<=MaxColumn; Column++, str=NULL)
   {
      if ( strtok( str, delim ) == NULL  )  break;
      NColumn++;
   }

   if ( NColumn < 1  ||  NColumn > MaxColumn )
   {
      fprintf( stderr, "ERROR : incorrect data columns = %d (min = 1, max = %d) !!\n", NColumn, MaxColumn );
      exit( 1 );
   }

   if ( TColumn < 0  ||  TColumn >= NColumn )
   {
      fprintf( stderr, "ERROR : incorrect targeted column = %d (min = 0, max = %d) !!\n", TColumn, NColumn-1 );
      exit( 1 );
   }


// get the number of data rows
   NRow = 0;
   rewind( File_In );
   for (int t=0; t<Header; t++)  fgets( Line, MaxLine, File_In );    // skip headers

   while ( fgets( Line, MaxLine, File_In ) != NULL )  NRow ++;

   if ( NRow < 1 )
   {
      fprintf( stderr, "ERROR : no data rows are found !!\n" );
      exit( 2 );
   }

   printf( "   Data Columns    = %d\n", NColumn );
   printf( "   Data Rows       = %d\n", NRow  );
   printf( "   Targeted Column = %d\n", TColumn );


// allocate data
   Key           = new real [NRow];
   Key_IdxTable  = new int  [NRow];
   Key_RankTable = new int  [NRow];
   Data          = new char [NRow][MaxLine];


// load the variable to be sorted (the "key")
   rewind( File_In );
   for (int t=0; t<Header; t++)  fgets( Line, MaxLine, File_In );    // skip headers

   for (int Row=0; Row<NRow; Row++)
   {
      fgets( Line, MaxLine, File_In );

      for (Column=0, str=Line; Column<=TColumn; Column++, str=NULL)
      {
         Temp = strtok( str, delim );

         if ( Temp != NULL )
         {
            if ( Column == TColumn )   Key[Row] = atof( Temp );
         }
         else
         {
            fprintf( stderr, "ERROR : incorrect read at column %d, row %d!!\n", Column+1, Header+Row+1 );
            exit( 2 );
         }
      }
   }

   fclose( File_In );
   delete [] Line;


// sort the key
   Heapsort( NRow, Key, Key_IdxTable );


// get the rank table
   for (int Row=0; Row<NRow; Row++)    Key_RankTable[ Key_IdxTable[Row] ] = Row;


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : SortKey



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while ( (c = getopt(argc, argv, "hi:o:n:c:")) != -1 )
   {
      switch ( c )
      {
         case 'i':   FileName_In  = optarg;        break;
         case 'o':   FileName_Out = optarg;        break;
         case 'n':   Header       = atoi(optarg);  break;
         case 'c':   TColumn      = atoi(optarg);  break;
         case 'h':
         case '?':   cerr << endl << "usage: " << argv[0]
                          << " [-h (for help)] [-i input filename] [-o output filename]"
                          << endl << "                              "
                          << " [-n number of header lines [1]] [-c targeted column [0]]"
                          << endl << endl;
                     exit( 1 );

      } // switch ( c )
   } // while ..


// set the name of the output file if it is not specified
   if ( FileName_Out == NULL )
   {
      FileName_Out = new char [1024];
      sprintf( FileName_Out, "%s_%s", FileName_In, "Ordered" );
   }


// check parameters
   if ( FileName_In == NULL )
   {
      fprintf( stderr, "ERROR : please provide the input filename (-i input filename) !!\n" );
      exit( 1 );
   }

   if ( FileName_Out == NULL )
   {
      fprintf( stderr, "ERROR : please provide the output filename (-o output filename) !!\n" );
      exit( 1 );
   }

   if ( fopen( FileName_In, "r" ) == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist (-i input filename) !!\n", FileName_In );
      exit( 1 );
   }

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  End
// Description :  End program
//-------------------------------------------------------------------------------------------------------
void End()
{

   delete [] Key;
   delete [] Key_IdxTable;
   delete [] Key_RankTable;
   delete [] Data;

} // FUNCTION : End



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort
// Description :  Use the Heapsort algorithm to sort the input array into ascending numerical order
//                --> An index table will also be constructed if "IdxTable != NULL"
//
// Note        :  Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  N        :  Size of Array
//                Array    :  Array to be sorted
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Heapsort( const int N, real Array[], int IdxTable[] )
{

// initialize the IdxTable
   if ( IdxTable != NULL )
      for (int t=0; t<N; t++)    IdxTable[t] = t;

// heap creation
   for (int L=N/2-1; L>=0; L--)  Heapsort_SiftDown( L, N-1, Array, IdxTable );

// retirement-and-promotion
   real Tmp_real;
   int  Tmp_int;
   for (int R=N-1; R>0; R--)
   {
      Tmp_real = Array[R];
      Array[R] = Array[0];
      Array[0] = Tmp_real;

      if ( IdxTable != NULL )
      {
         Tmp_int     = IdxTable[R];
         IdxTable[R] = IdxTable[0];
         IdxTable[0] = Tmp_int;
      }

      Heapsort_SiftDown( 0, R-1, Array, IdxTable );
   }

} // FUNCTION : Mis_Heapsort



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort_SiftDown
// Description :  Sift-down process for the Heapsort algorithm
//
// Note        :  Ref : Numerical Recipes Chapter 8.3 - 8.4
//
// Parameter   :  L        :  Left  range of the sift-down
//                R        :  Right range of the sift-down
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
void Heapsort_SiftDown( const int L, const int R, real Array[], int IdxTable[] )
{

   int  Idx_up    = L;
   int  Idx_down  = 2*Idx_up + 1;
   real Target    = Array[Idx_up];
   int  TargetIdx = ( IdxTable == NULL ) ? -1 : IdxTable[Idx_up];

   while ( Idx_down <= R )
   {
//    find the better employee
      if ( Idx_down < R  &&  Array[Idx_down+1] > Array[Idx_down] )   Idx_down ++;

//    terminate the sift-down process if the target (supervisor) is better than both its employees
      if ( Target >= Array[Idx_down] )    break;

//    otherwise, promote the better employee
      Array[Idx_up] = Array[Idx_down];
      if ( IdxTable != NULL )    IdxTable[Idx_up] = IdxTable[Idx_down];

//    prepare the next sift-down operation
      Idx_up   = Idx_down;
      Idx_down = 2*Idx_up + 1;
   }

// put target at its best position
   Array[Idx_up] = Target;
   if ( IdxTable != NULL )    IdxTable[Idx_up] = TargetIdx;

} // FUNCTION : Heapsort_SiftDown



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   fprintf( stdout, "Invoking %s ...\n", "GAMER_MakeDataInOrder" );

   ReadOption( argc, argv );

   SortKey();

   LoadData();

   DumpData();

   End();

   fprintf( stdout, "Invoking %s ... done\n", "GAMER_MakeDataInOrder" );

   return 0;

} // FUNCTION : main
