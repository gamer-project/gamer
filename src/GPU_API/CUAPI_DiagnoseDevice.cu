#include "CUAPI.h"

void Aux_GetCPUInfo( const char *FileName );

#ifdef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_DiagnoseDevice
// Description :  Take a diagnosis of each GPU
//-------------------------------------------------------------------------------------------------------
void CUAPI_DiagnoseDevice()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// get the hostname and PID of each process
   const int PID = getpid();
   char Host[1024];
   gethostname( Host, 1024 );


// get the number of devices
   int DeviceCount;
   CUDA_CHECK_ERROR(  cudaGetDeviceCount( &DeviceCount )  );

   if ( DeviceCount == 0 )
      Aux_Error( ERROR_INFO, "no devices supporting CUDA at MPI_Rank %2d (host = %8s) !!\n", MPI_Rank, Host );


// get the device ID
   int GetDeviceID = 999;
   CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );


// load the device properties
   cudaDeviceProp DeviceProp;
   CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );


// get the number of cores per multiprocessor
   int NCorePerMP;
   if      ( DeviceProp.major == 2  &&  DeviceProp.minor == 0 )  NCorePerMP =  32;
   else if ( DeviceProp.major == 2  &&  DeviceProp.minor == 1 )  NCorePerMP =  48;
   else if ( DeviceProp.major == 3 )                             NCorePerMP = 192;
   else if ( DeviceProp.major == 5 )                             NCorePerMP = 128;
   else if ( DeviceProp.major == 6 )                             NCorePerMP =  64;
   else if ( DeviceProp.major == 7 )                             NCorePerMP =  64;
   else if ( DeviceProp.major == 8  &&  DeviceProp.minor == 0 )  NCorePerMP =  64;
   else if ( DeviceProp.major == 8  &&  DeviceProp.minor == 6 )  NCorePerMP = 128;
   else if ( DeviceProp.major == 8  &&  DeviceProp.minor == 9 )  NCorePerMP = 128;
   else if ( DeviceProp.major == 9 )                             NCorePerMP = 128;
   else
      fprintf( stderr, "WARNING : unable to determine the number of cores per multiprocessor for version %d.%d ...\n",
               DeviceProp.major, DeviceProp.minor );


// record the device properties
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Note", OUTPUT_DIR );

   if ( MPI_Rank == 0 )
   {
       FILE *Note = fopen( FileName, "a" );
       fprintf( Note, "Device Diagnosis\n" );
       fprintf( Note, "***********************************************************************************\n" );
       fclose( Note );
   }

   for (int YourTurn=0; YourTurn<MPI_NRank; YourTurn++)
   {
      if ( MPI_Rank == YourTurn )
      {
         int DriverVersion = 0, RuntimeVersion = 0;
         CUDA_CHECK_ERROR(  cudaDriverGetVersion( &DriverVersion )  );
         CUDA_CHECK_ERROR(  cudaRuntimeGetVersion( &RuntimeVersion )  );

         FILE *Note = fopen( FileName, "a" );
         if ( MPI_Rank != 0 )   fprintf( Note, "\n\n" );
         fprintf( Note, "MPI_Rank = %3d, hostname = %10s, PID = %5d\n\n", MPI_Rank, Host, PID );
         fprintf( Note, "CPU Info :\n" );
         fflush( Note );

         Aux_GetCPUInfo( FileName );

         fprintf( Note, "\n" );
         fprintf( Note, "GPU Info :\n" );
         fprintf( Note, "Number of GPUs                          : %d\n"     , DeviceCount );
         fprintf( Note, "GPU ID                                  : %d\n"     , GetDeviceID );
         fprintf( Note, "GPU Name                                : %s\n"     , DeviceProp.name );
         fprintf( Note, "CUDA Driver Version                     : %d.%d\n"  , DriverVersion/1000, DriverVersion%100 );
         fprintf( Note, "CUDA Runtime Version                    : %d.%d\n"  , RuntimeVersion/1000, RuntimeVersion%100 );
         fprintf( Note, "CUDA Major Revision Number              : %d\n"     , DeviceProp.major );
         fprintf( Note, "CUDA Minor Revision Number              : %d\n"     , DeviceProp.minor );
         fprintf( Note, "Clock Rate                              : %f GHz\n" , DeviceProp.clockRate/1.0e6 );
         fprintf( Note, "Global Memory Size                      : %ld MB\n" , (long)DeviceProp.totalGlobalMem/1024/1024 );
         fprintf( Note, "Constant Memory Size                    : %ld KB\n" , (long)DeviceProp.totalConstMem/1024 );
         fprintf( Note, "Shared Memory Size per Block            : %ld KB\n" , (long)DeviceProp.sharedMemPerBlock/1024 );
         fprintf( Note, "Number of Registers per Block           : %d\n"     , DeviceProp.regsPerBlock );
         fprintf( Note, "Warp Size                               : %d\n"     , DeviceProp.warpSize );
         fprintf( Note, "Number of Multiprocessors:              : %d\n"     , DeviceProp.multiProcessorCount );
         fprintf( Note, "Number of FP32 Cores per Multiprocessor : %d\n"     , NCorePerMP );
         fprintf( Note, "Total Number of Cores:                  : %d\n"     , DeviceProp.multiProcessorCount*NCorePerMP );
         fprintf( Note, "Max Number of Threads per Block         : %d\n"     , DeviceProp.maxThreadsPerBlock );
         fprintf( Note, "Max Size of the Block X-Dimension       : %d\n"     , DeviceProp.maxThreadsDim[0] );
         fprintf( Note, "Max Size of the Grid X-Dimension        : %d\n"     , DeviceProp.maxGridSize[0] );
         fprintf( Note, "Concurrent Copy and Execution           : %s\n"     , DeviceProp.asyncEngineCount>0  ? "Yes" : "No" );
         fprintf( Note, "Concurrent Up/Downstream Copies         : %s\n"     , DeviceProp.asyncEngineCount==2 ? "Yes" : "No" );
#        if ( CUDART_VERSION >= 3000 )
         fprintf( Note, "Concurrent Kernel Execution             : %s\n"     , DeviceProp.concurrentKernels ? "Yes" : "No" );
#        endif
#        if ( CUDART_VERSION >= 3010 )
         fprintf( Note, "GPU has ECC Support Enabled             : %s\n"     , DeviceProp.ECCEnabled ? "Yes" : "No" );
#        endif

         fclose( Note );
       }

       MPI_Barrier( MPI_COMM_WORLD );

   } // for (int YourTurn=0; YourTurn<NGPU; YourTurn++)

   if ( MPI_Rank == 0 )
   {
      FILE *Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fclose( Note );

      Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
   }

} // FUNCTION : CUAPI_DiagnoseDevice



#endif // #ifdef GPU
