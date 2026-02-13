#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_GetCPUInfo
// Description :  Record the CPU information
//
// Parameter   :  FileName : Name of the output file
//-------------------------------------------------------------------------------------------------------
void Aux_GetCPUInfo( const char *FileName )
{

// CPU info reporting is not currently supported on macOS
#  ifdef __APPLE__
   Aux_Message( stderr, "WARNING : function \"%s\" is not supported on macOS !!\n", __FUNCTION__ );
   return;
#  endif


   FILE *Note = fopen( FileName, "a" );
   char *line = NULL;
   size_t len = 0;
   char String[2][MAX_STRING];
   char Trash[MAX_STRING];
   int SocketNow = -1, SocketPrevious = -1;
   int CorePerSocket = 0, NSocket = 0;
   bool GotFirstCPUInfo = false;


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

      if (  strcmp( String[0], "physical" ) == 0  &&  strcmp( String[1], "id" ) == 0 )
      {
         sscanf( line, "%s%s%s%d", String[0], String[1], Trash, &SocketNow );
         if ( SocketNow != SocketPrevious )
         {
            SocketPrevious = SocketNow;
            NSocket++;
         }
      }

      if ( GotFirstCPUInfo )   continue;

      if (  strcmp( String[0], "model" ) == 0  &&  strcmp( String[1], "name" ) == 0  )
      {
         memcpy( line, "CPU Type  ", 10 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cpu" ) == 0  &&  strcmp( String[1], "MHz" ) == 0  )
      {
         memcpy( line, "CPU MHz", 7 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cache" ) == 0  &&  strcmp( String[1], "size" ) == 0  )
      {
         memcpy( line, "Cache Size", 10 );
         fprintf( Note, "%s", line );
      }

      if (  strcmp( String[0], "cpu" ) == 0  &&  strcmp( String[1], "cores" ) == 0  )
      {
         memcpy( line, "CPU Cores", 9 );
         fprintf( Note, "%s", line );
         sscanf( line, "%s%s%s%d", String[0], String[1], Trash, &CorePerSocket );
         GotFirstCPUInfo = true;
      }
   }

   if ( line != NULL )
   {
      free( line );
      line = NULL;
   }

   fprintf( Note, "%-16s: %d\n", "Socket(s)", NSocket );
// assuming the CPUs in the node are the same
   fprintf( Note, "%-16s: %d\n", "Core(s) per Node", CorePerSocket*NSocket );

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
