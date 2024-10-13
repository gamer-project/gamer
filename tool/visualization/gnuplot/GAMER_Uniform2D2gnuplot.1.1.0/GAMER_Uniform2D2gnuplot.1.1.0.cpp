#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
using namespace std;

// single/double precision
#ifdef FLOAT8
   typedef double real;
#else
   typedef float  real;
#endif

// extreme values
#ifndef __INT_MAX__
   #define __INT_MAX__     2147483647
#endif


void ReadOption( int argc, char **argv );
void GetRange();
void LoadData();
void End();

const char *delim   = " \n\t";      // delimiter for "strtok"
const int  MaxToken = 20;           // maximum number of tokens
const int  MaxLine  = 1024;         // maximum number of characters per line

char *FileName_In  = NULL;
char *FileName_Out = NULL;

int   Max[3] = { -__INT_MAX__, -__INT_MAX__, -__INT_MAX__ };
int   Min[3] = { +__INT_MAX__, +__INT_MAX__, +__INT_MAX__ };
int   Scale, Dim1, Dim2, NGrid1, NGrid2, NComp, Target;
int   Header = 1;
real  ***Data;




//-------------------------------------------------------------------------------------------------------
// Function    :  DumpData
// Description :  Dump data
//-------------------------------------------------------------------------------------------------------
void DumpData()
{

   fprintf( stdout, "%s ... ", __FUNCTION__ );   fflush( stdout );


   FILE *File_Check = fopen( FileName_Out, "r" );
   if ( File_Check != NULL )
   {
      fprintf( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n", FileName_Out );
      fclose( File_Check );
   }


   FILE *File = fopen( FileName_Out, "w" );
   int  mn[3];
   char Dir1, Dir2;
   real *Out;

   switch ( Target )
   {
      case 0 : Dir1 = 'j';   Dir2 = 'k';   break;
      case 1 : Dir1 = 'i';   Dir2 = 'k';   break;
      case 2 : Dir1 = 'i';   Dir2 = 'j';   break;
   }

   fprintf( File, "%4c %4c", Dir1, Dir2 );
   for (int v=0; v<NComp; v++)   fprintf( File, "     COLUMN %2d", v );
   fprintf( File, "\n" );

   for (int n=0; n<NGrid1; n++)
   {
      mn[0] = n*Scale;

      for (int m=0; m<NGrid2; m++)
      {
         mn[1] = m*Scale;
         Out   = Data[m][n];

         fprintf( File, "%4d %4d", mn[0], mn[1] );
         for (int v=0; v<NComp; v++)   fprintf( File, " %13.6e", Out[v] );
         fprintf( File, "\n" );
      }

      fprintf( File, "\n" );
   }

   fclose( File );


   fprintf( stdout, "done\n" );   fflush( stdout );

} // FUNCTION : DumpData



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load data from the input file
//-------------------------------------------------------------------------------------------------------
void LoadData()
{

   fprintf( stdout, "%s ... ", __FUNCTION__ );   fflush( stdout );


   const int MaxSize = 1024;
   int   NLine = 0;
   int   m, n, dim, ijk[3];
   char *Line = new char [MaxSize];
   real *In   = new real [NComp];
   char *Temp  = NULL, *str = NULL;


   FILE *File = fopen( FileName_In, "r" );

   if ( File == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }


// skip header
   for (int t=0; t<Header; t++)  fgets( Line, MaxSize, File );


// load data
   while ( fgets( Line, MaxSize, File ) != NULL )
   {
      NLine ++;

//    a. load the cell indices
      for (dim=0, str=Line; dim<3; dim++, str=NULL)
      {
         Temp = strtok( str, delim );

         if ( Temp != NULL )
            ijk[dim] = atoi( Temp );
         else
         {
            fprintf( stderr, "ERROR : incorrect read at line %d, token %d!!\n", NLine, dim+1 );
            exit(2);
         }
      }


//    b. load the physical variables
      for (int v=0; v<NComp; v++)
      {
         Temp = strtok( NULL, delim );

         if ( Temp != NULL )
            In[v] = atof( Temp );
         else
         {
            fprintf( stderr, "ERROR : incorrect read at line %d, token %d!!\n", NLine, v+4 );
            exit(2);
         }
      }

      m = ( ijk[Dim2] - Min[Dim2] ) / Scale;
      n = ( ijk[Dim1] - Min[Dim1] ) / Scale;

      for (int v=0; v<NComp; v++)   Data[m][n][v] = In[v];

   } // while ( fgets( Line, MaxSize, File ) != NULL )


   fclose( File );
   delete [] Line;
   delete [] In;


   fprintf( stdout, "done\n" );   fflush( stdout );

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  GetRange
// Description :  Determine the data range by scanning the input file
//-------------------------------------------------------------------------------------------------------
void GetRange()
{

   fprintf( stdout, "%s ... \n", __FUNCTION__ );   fflush( stdout );


   const int MaxSize = 1024;
   int   Token, NToken, ijk[3], NRow;
   int   Min2[3] = { +__INT_MAX__, +__INT_MAX__, +__INT_MAX__ };
   char *Line = new char [MaxSize];
   char *str  = NULL;


   FILE *File = fopen( FileName_In, "r" );

   if ( File == NULL )
   {
      fprintf( stderr, "ERROR : input file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }


// get the number of data columns
   for (int t=0; t<Header; t++)  fgets( Line, MaxSize, File );    // skip headers

   fgets( Line, MaxLine, File );
   NToken = 0;

   for (Token=0, str=Line; Token<=MaxToken; Token++, str=NULL)
   {
      if ( strtok( str, delim ) == NULL  )  break;
      NToken++;
   }

   if ( NToken < 4  ||  NToken > MaxToken )
   {
      fprintf( stderr, "ERROR : incorrect data columns = %d (min = 4, max = %d) !!\n", NToken, MaxToken );
      exit(1);
   }

   NComp = NToken - 3;


// get the data range and interval
   NRow = 0;
   rewind( File );
   for (int t=0; t<Header; t++)  fgets( Line, MaxSize, File );    // skip headers

   while ( fgets( Line, MaxSize, File ) != NULL )
   {
      NRow ++;
      sscanf( Line, "%d%d%d", &ijk[0], &ijk[1], &ijk[2] );

      for (int d=0; d<3; d++)
      {
         if      ( ijk[d] <= Min [d] )  Min [d] = ijk[d];
         else if ( ijk[d] <  Min2[d] )  Min2[d] = ijk[d];
         if      ( ijk[d] >  Max [d] )  Max [d] = ijk[d];
      }

   } // while ( fgets( Line, MaxSize, File ) != NULL )

   for (int d=0; d<3; d++)
   {
      if ( Min[d] == Max[d] )
      {
         Target = d;
         Dim1   = ( d == 0 ) ? 1 : 0;
         Dim2   = ( d == 2 ) ? 1 : 2;
         break;
      }
   }

   Scale  = Min2[Dim1] - Min[Dim1];
   NGrid1 = ( Max[Dim1] - Min[Dim1] )/Scale + 1;
   NGrid2 = ( Max[Dim2] - Min[Dim2] )/Scale + 1;


// allocate memory
   Data       = new real** [NGrid2];
   Data[0]    = new real*  [NGrid2*NGrid1];
   Data[0][0] = new real   [NGrid2*NGrid1*NComp];

   for (int m=1; m<NGrid2; m++)   Data[m]    = Data[0  ]    + m*NGrid1;

   for (int n=1; n<NGrid1; n++)   Data[0][n] = Data[0  ][0] + n*NComp;

   for (int m=1; m<NGrid2; m++)   Data[m][0] = Data[m-1][0] + NGrid1*NComp;

   for (int m=1; m<NGrid2; m++)
   for (int n=1; n<NGrid1; n++)   Data[m][n] = Data[m  ][0] + n*NComp;


   printf( "   Data Columns = %d\n", NComp );
   printf( "   Data Rows    = %d\n", NRow  );
   printf( "   Target       = %c slice\n", 88+Target );
   printf( "   Range 1      = [%d ... %d]\n", Min[Dim1], Max[Dim1] );
   printf( "   Range 2      = [%d ... %d]\n", Min[Dim2], Max[Dim2] );
   printf( "   NGrid 1      = %d\n", NGrid1 );
   printf( "   NGrid 2      = %d\n", NGrid2 );
   printf( "   Grid Scale   = %d\n", Scale );


   fclose( File );
   delete [] Line;


   fprintf( stdout, "%s ... done\n", __FUNCTION__ );   fflush( stdout );

} // FUNCTION : GetRange



//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   int c;

   while ( (c = getopt(argc, argv, "hi:o:n:")) != -1 )
   {
      switch ( c )
      {
         case 'i':   FileName_In  = optarg;        break;
         case 'o':   FileName_Out = optarg;        break;
         case 'n':   Header       = atoi(optarg);  break;
         case 'h':
         case '?':   cerr << endl << "usage: " << argv[0]
                          << " [-h (for help)] [-i input filename] [-o output filename]"
                          << endl << "                                "
                          << " [-n header lines [1]]"
                          << endl << endl;
                     exit( 1 );

      } // switch ( c )
   } // while ..


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

   delete [] Data[0][0];
   delete [] Data[0];
   delete [] Data;

} // FUNCTION : End



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   cout << "Invoking GAMER_Uniform2D2gnuplot ..." << endl;

   ReadOption( argc, argv );

   GetRange();

   LoadData();

   DumpData();

   End();

   cout << "GAMER_Uniform2D2gnuplot terminated successfully" << endl;

   return 0;

} // FUNCTION : main



