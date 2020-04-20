#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_GetCPUInfo
// Description :  Record the CPU information
//
// Parameter   :  FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Aux_GetCPUInfo( const char *FileName )
{

   FILE *Note = fopen( FileName, "a" );
   char *line = NULL;
   size_t len = 0;
   char String[2][100];


// 1. get the CPU info
   const char *CPUInfo_Path = "/proc/cpuinfo";

   if ( !Aux_CheckFileExist(CPUInfo_Path) )
   {
      Aux_Message( stderr, "WARNING : CPU information file \"%s\" does not exist !!\n", CPUInfo_Path );
      return;
   }

   FILE *CPUInfo = fopen( CPUInfo_Path, "r" );

   while ( getline(&line, &len, CPUInfo) != -1 )
   {
      sscanf( line, "%s%s", String[0], String[1] );

      if (  strcmp( String[0], "model" ) == 0  &&  strcmp( String[1], "name" ) == 0  )
      {
         strncpy( line, "CPU Type  ", 10 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cpu" ) == 0  &&  strcmp( String[1], "MHz" ) == 0  )
      {
         strncpy( line, "CPU MHz", 7 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cache" ) == 0  &&  strcmp( String[1], "size" ) == 0  )
      {
         strncpy( line, "Cache Size", 10 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cpu" ) == 0  &&  strcmp( String[1], "cores" ) == 0  )
      {
         strncpy( line, "CPU Cores", 9 );
         fprintf( Note, "%s", line );
         break;
      }
   }

   if ( line != NULL )
   {
      free( line );
      line = NULL;
   }

   fclose( CPUInfo );


// 2. get the memory info
   const char *MemInfo_Path = "/proc/meminfo";

   if ( !Aux_CheckFileExist(MemInfo_Path) )
   {
      Aux_Message( stderr, "WARNING : memory information file \"%s\" does not exist !!\n", MemInfo_Path );
      return;
   }

   FILE *MemInfo = fopen( MemInfo_Path, "r" );

   while ( getline(&line, &len, MemInfo) != -1 )
   {
      sscanf( line, "%s%s", String[0], String[1] );

      if (  strncmp( String[0], "MemTotal", 8 ) == 0  )
      {
         fprintf( Note, "Total Memory    : %4.1f GB\n", atof( String[1] )/(1<<20) );
         break;
      }
   }

   if ( line != NULL )
   {
      free( line );
      line = NULL;
   }

   fclose( MemInfo );
   fclose( Note );

} // FUNCTION : Aux_GetCPUInfo
