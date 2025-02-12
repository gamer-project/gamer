#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_GetMemInfo
// Description :  Record the memory consumption
//
// Note        :  1. This function will record the following information from the file "/proc/[pid]/status"
//                   (1) VmSize/Peak : current/peak virtual  memory size
//                   (2) VmRSS/HWM   : current/peak physical memory size
//                2. Only the maximum values among all MPI ranks will be recorded
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_GetMemInfo()
{

// memory reporting is not currently supported on macOS
#  ifdef __APPLE__
   return;
#  endif


   const int PID   = getpid();
   const int NInfo = 4; // number of memory information to be recorded (VmSize/Peak, VmRSS/HWM)

   static bool FirstTime=true;
   char   FileName_Record[2*MAX_STRING], FileName_Status[MAX_STRING], Useless[NInfo][MAX_STRING], *line=NULL;
   char   VmSize[MAX_STRING], VmPeak[MAX_STRING], VmRSS[MAX_STRING], VmHWM[MAX_STRING];
   bool   GetVmSize=false, GetVmPeak=false, GetVmRSS=false, GetVmHWM=false;
   double Vm_double[NInfo], Vm_max[NInfo], Vm_sum[NInfo];
   size_t len=0;

   sprintf( FileName_Record, "%s/Record__MemInfo", OUTPUT_DIR );


// 1. read memory information from the file "FileName_Status"
   sprintf( FileName_Status, "/proc/%d/status", PID );

   if ( !Aux_CheckFileExist(FileName_Status) )
   {
      Aux_Message( stderr, "WARNING : PID status file \"%s\" does not exist (Rank %d) !!\n",
                   FileName_Status, MPI_Rank );
      return;
   }

   FILE *StatusFile = fopen( FileName_Status, "r" );

   while (  !GetVmSize  ||  !GetVmPeak  ||  !GetVmRSS  ||  !GetVmHWM )
   {
      if ( getline( &line, &len, StatusFile ) == -1 )
      {
         Aux_Message( stderr, "WARNING : some memory information is not found at Rank %d ", MPI_Rank );
         Aux_Message( stderr, "(VmSize: %s, VmPeak %s, VmRSS: %s, VmHWM %s)\n",
                      (GetVmSize) ? "OK" : "NO", (GetVmPeak) ? "OK" : "NO",
                      (GetVmRSS ) ? "OK" : "NO", (GetVmHWM ) ? "OK" : "NO" );
         break;
      }

      if      ( strncmp( line, "VmSize:", 7 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmSize, Useless[1] );
         GetVmSize = true;
      }

      else if ( strncmp( line, "VmPeak:", 7 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmPeak, Useless[1] );
         GetVmPeak = true;
      }

      else if ( strncmp( line, "VmRSS:", 6 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmRSS, Useless[1] );
         GetVmRSS = true;
      }

      else if ( strncmp( line, "VmHWM:", 6 ) == 0 )
      {
         sscanf( line, "%s%s%s", Useless[0], VmHWM, Useless[1] );
         GetVmHWM = true;
      }
   } // while (  !GetVmSize  ||  !GetVmRSS  )

   fclose( StatusFile );

   if ( line != NULL )  free( line );


// 2. gather information from all ranks
   Vm_double[0] = atof( VmSize );
   Vm_double[1] = atof( VmPeak );
   Vm_double[2] = atof( VmRSS  );
   Vm_double[3] = atof( VmHWM  );

   MPI_Reduce( Vm_double, Vm_max, NInfo, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
   MPI_Reduce( Vm_double, Vm_sum, NInfo, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


// 3. record memory information
   if ( MPI_Rank == 0 )
   {
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName_Record) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName_Record );

         FirstTime = false;

         FILE *File_Record = fopen( FileName_Record, "a" );
         fprintf( File_Record, "# Vir_Max  : maximum virtual  memory size of a single process at the present\n" );
         fprintf( File_Record, "# Vir_Sum  : total   virtual  memory size of all processes    at the present\n" );
         fprintf( File_Record, "# Vir_Peak : maximum virtual  memory size of a single process during the entire simulation\n" );
         fprintf( File_Record, "# Phy_Max  : maximum physical memory size of a single process at the present\n" );
         fprintf( File_Record, "# Phy_Sum  : total   physical memory size of all processes    at the present\n" );
         fprintf( File_Record, "# Phy_Peak : maximum physical memory size of a single process during the entire simulation\n" );
         fprintf( File_Record, "#------------------------------------------------------------------------------------------\n\n" );
         fprintf( File_Record, "#%13s%14s%s%20s%20s%20s%20s%20s%20s\n",
                  "Time", "Step", " ",
                  "Vir_Max (MB)", "Vir_Sum (MB)", "Vir_Peak (MB)",
                  "Phy_Max (MB)", "Phy_Sum (MB)", "Phy_Peak (MB)" );
         fclose( File_Record );
      }

      FILE *File_Record = fopen( FileName_Record, "a" );
      fprintf( File_Record, "%14.7e%14ld%20.2f%20.2f%20.2f%20.2f%20.2f%20.2f\n",
               Time[0], Step,
               Vm_max[0]/1024.0, Vm_sum[0]/1024.0, Vm_max[1]/1024.0,
               Vm_max[2]/1024.0, Vm_sum[2]/1024.0, Vm_max[3]/1024.0 );
      fclose( File_Record );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_GetMemInfo


