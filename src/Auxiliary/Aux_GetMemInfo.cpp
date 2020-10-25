#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_GetMemInfo
// Description :  Record the memory consumption
//
// Note        :  1. This function will record the following information from the file "/proc/[pid]/status"
//                   (1) VmSize : current virtual memory size
//                   (2) VmRSS  : current resident set size
//                2. Only the maximum values among all MPI ranks will be recorded
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_GetMemInfo()
{

   const char FileName_Record[] = "Record__MemInfo";
   const int  StrSize           = 128;
   const int  PID               = getpid();

   static bool FirstTime=true;
   char   FileName_Status[StrSize], Useless[2][StrSize], *line=NULL;
   char   VmSize[StrSize], VmRSS[StrSize];
   bool   GetVmSize=false, GetVmRSS=false;
   double Vm_double[2], Vm_max[2], Vm_sum[2];
   size_t len=0;


// 1. read memory information from the file "FileName_Status"
   sprintf( FileName_Status, "/proc/%d/status", PID );

   if ( !Aux_CheckFileExist(FileName_Status) )
   {
      Aux_Message( stderr, "WARNING : PID status file \"%s\" does not exist (Rank %d) !!\n",
                   FileName_Status, MPI_Rank );
      return;
   }

   FILE *StatusFile = fopen( FileName_Status, "r" );

   while (  !GetVmSize  ||  !GetVmRSS  )
   {
      if ( getline( &line, &len, StatusFile ) == -1 )
      {
         Aux_Message( stderr, "WARNING : some memory information is not found at Rank %d ", MPI_Rank );
         Aux_Message( stderr, "(VmSize: %s, VmRSS: %s)\n", (GetVmSize) ? "OK" : "NO", (GetVmRSS ) ? "OK" : "NO" );
         break;
      }

      if      ( strncmp( line, "VmSize:", 7 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmSize, Useless[1] );
         GetVmSize = true;
      }

      else if ( strncmp( line, "VmRSS:", 6 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmRSS, Useless[1] );
         GetVmRSS = true;
      }
   } // while (  !GetVmSize  ||  !GetVmRSS  )

   fclose( StatusFile );

   if ( line != NULL )  free( line );


// 2. gather information from all ranks
   Vm_double[0] = atof( VmSize );
   Vm_double[1] = atof( VmRSS  );

   MPI_Reduce( Vm_double, Vm_max, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Vm_double, Vm_sum, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// 3. record memory information
   if ( MPI_Rank == 0 )
   {
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName_Record) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName_Record );

         FirstTime = false;

         FILE *File_Record = fopen( FileName_Record, "a" );
         fprintf( File_Record, "#%13s%14s%s%20s%20s%20s%20s\n", "Time", "Step", " ", "Virtual_Max (MB)",
                  "Virtual_Sum (MB)", "Resident_Max (MB)", "Resident_Sum (MB)" );
         fclose( File_Record );
      }

      FILE *File_Record = fopen( FileName_Record, "a" );
      fprintf( File_Record, "%14.7e%14ld%20.2f%20.2f%20.2f%20.2f\n",
               Time[0], Step, Vm_max[0]/1024.0, Vm_sum[0]/1024.0, Vm_max[1]/1024.0, Vm_sum[1]/1024.0 );
      fclose( File_Record );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_GetMemInfo


