#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif
#ifdef LAOHU
extern "C" { int GetFreeGpuDevID( int, int ); }
#endif

#ifdef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetDevice
// Description :  Set the active device
//
// Parameter   :  Mode :    -3 --> set by the gpudevmgr library on the NAOC Laohu cluster
//                          -2 --> set automatically by CUDA (must work with the "compute-exclusive mode")
//                          -1 --> set by MPI ranks : SetDeviceID = MPI_Rank % DeviceCount
//                       >=  0 --> set to "Mode"
//-------------------------------------------------------------------------------------------------------
void CUAPI_SetDevice( const int Mode )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "CUAPI_SetDevice ... \n" );


// check
#  ifdef LAOHU
   if ( Mode < -3 )     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Mode", Mode );
   if ( Mode != -3  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : \"OPT__GPUID_SELECT != -3\" on the Laohu cluster !?\n" );
#  else
   if ( Mode < -2 )     Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Mode", Mode );
#  endif


// get the hostname of each MPI process
   char Host[1024];
   gethostname( Host, 1024 );


// verify that there are GPU supporting CUDA
   int DeviceCount;
   CUDA_CHECK_ERROR(  cudaGetDeviceCount( &DeviceCount )  );

   if ( DeviceCount == 0 )
      Aux_Error( ERROR_INFO, "no devices support CUDA at MPI_Rank %2d (host = %8s) !!\n", MPI_Rank, Host );


// set the device ID
   void **d_TempPtr = NULL;
   int SetDeviceID, GetDeviceID = 999;
   cudaDeviceProp DeviceProp;

   switch ( Mode )
   {
#     ifdef LAOHU
      case -3:
         SetDeviceID = GetFreeGpuDevID( DeviceCount, MPI_Rank );

         if ( SetDeviceID < DeviceCount )
            CUDA_CHECK_ERROR(  cudaSetDevice( SetDeviceID )  );

         else
            Aux_Error( ERROR_INFO, "SetDeviceID (%d) >= DeviceCount (%d) at MPI_Rank %2d (host = %8s) !!\n",
                       SetDeviceID, DeviceCount, MPI_Rank, Host );
         break;
#     endif


      case -2:
         CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_TempPtr, sizeof(int) )  );  // to set the GPU ID
         CUDA_CHECK_ERROR(  cudaFree( d_TempPtr )  );

//       make sure that the "exclusive" compute mode is adopted
         CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );
         CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );

         if ( DeviceProp.computeMode != cudaComputeModeExclusive )
         {
            Aux_Message( stderr, "WARNING : \"exclusive\" compute mode is NOT enabled for \"%s\" at Rank %2d",
                         "OPT__GPUID_SELECT == -2", MPI_Rank );
            Aux_Message( stderr, " (host=%8s) !!\n", Host );
         }
         break;


      case -1:
         SetDeviceID = MPI_Rank % DeviceCount;
         CUDA_CHECK_ERROR(  cudaSetDevice( SetDeviceID )  );

         if ( MPI_NRank > 1  &&  MPI_Rank == 0 )
         {
            Aux_Message( stderr, "WARNING : please make sure that different MPI ranks will use different GPUs " );
            Aux_Message( stderr, "for \"%s\" !!\n", "OPT__GPUID_SELECT == -1" );
         }
         break;


      default:
         SetDeviceID = Mode;

         if ( SetDeviceID < DeviceCount )
            CUDA_CHECK_ERROR(  cudaSetDevice( SetDeviceID )  );

         else
            Aux_Error( ERROR_INFO, "SetDeviceID (%d) >= DeviceCount (%d) at MPI_Rank %2d (host = %8s) !!\n",
                       SetDeviceID, DeviceCount, MPI_Rank, Host );

         if ( MPI_NRank > 1  &&  MPI_Rank == 0 )
         {
            Aux_Message( stderr, "WARNING : please make sure that different MPI ranks will use different GPUs " );
            Aux_Message( stderr, "for \"%s\" !!\n", "OPT__GPUID_SELECT == -1" );
         }
         break;
   } // switch ( Mode )


// check
// (0) load the device properties and the versions of CUDA and driver
   int DriverVersion = 0, RuntimeVersion = 0;
   CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );
   CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );
   CUDA_CHECK_ERROR(  cudaDriverGetVersion( &DriverVersion )  );
   CUDA_CHECK_ERROR(  cudaRuntimeGetVersion( &RuntimeVersion )  );


// (1) verify the device version
   if ( DeviceProp.major < 1 )
      Aux_Error( ERROR_INFO, "\ndevice major version < 1 at MPI_Rank %2d (host = %8s) !!\n", MPI_Rank, Host );

   if ( Mode >= -1 )
   {
//    (2) verify that the device ID is properly set
      if ( GetDeviceID != SetDeviceID )
         Aux_Error( ERROR_INFO, "GetDeviceID (%d) != SetDeviceID (%d) at MPI_Rank %2d (host = %8s) !!\n",
                    GetDeviceID, SetDeviceID, MPI_Rank, Host );

//    (3) verify that the adopted ID is accessible
      CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_TempPtr, sizeof(int) )  );
      CUDA_CHECK_ERROR(  cudaFree( d_TempPtr )  );
   }


// (4) verify the capability of double precision
#  ifdef FLOAT8
   if ( DeviceProp.major < 2  &&  DeviceProp.minor < 3 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" at MPI_Rank %2d (host = %8s) does not support FLOAT8 !!\n",
                 DeviceProp.name, MPI_Rank, Host );
#  endif


// (5) verify the GPU architecture
#  if   ( GPU_ARCH == FERMI )
   if ( DeviceProp.major != 2 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Fermi architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == KEPLER )
   if ( DeviceProp.major != 3 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Kepler architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == MAXWELL )
   if ( DeviceProp.major != 5 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Maxwell architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == PASCAL )
   if ( DeviceProp.major != 6 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Pascal architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  else
#  error : UNKNOWN GPU_ARCH !!
#  endif // GPU_ARCH


// (6) some options are not supported
// (6-1) fluid solver
#  if ( MODEL == HYDRO )
#  if ( FLU_SCHEME == WAF  &&  defined FLOAT8 )
   if ( RuntimeVersion < 3020 )
      Aux_Error( ERROR_INFO, "double-precision WAF scheme is not supported in CUDA version < 3.2 !!\n" );
#  endif

#  if ( FLU_SCHEME == WAF  &&  RSOLVER == EXACT )
#     error : ERROR : Currently WAF scheme does not support the exact Riemann solver !!;
#  endif

#  if (  defined FLOAT8  &&  CHECK_INTERMEDIATE == EXACT  && \
         ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
      if ( RuntimeVersion < 3020 )
         Aux_Error( ERROR_INFO, "CHECK_INTERMEDIATE == EXACT + FLOAT8 is not supported in CUDA < 3.2 !!" );
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  endif // #if ( MODEL == HYDRO )

// (6-2) SOR Poisson solver
#  if ( POT_SCHEME == SOR )
#     ifdef SOR_USE_SHUFFLE
      if ( DeviceProp.warpSize != 32 )
         Aux_Error( ERROR_INFO, "warp size (%d) != 32 !!\n", DeviceProp.warpSize );

      if ( DeviceProp.maxThreadsPerBlock > 1024 )
         Aux_Error( ERROR_INFO, "maximum number of threads per block (%d) > 1024 !!\n", DeviceProp.maxThreadsPerBlock );
#     endif

#     ifdef SOR_USE_PADDING
      if ( DeviceProp.warpSize != 32 )
         Aux_Error( ERROR_INFO, "warp size (%d) != 32 !!\n", DeviceProp.warpSize );

      if ( POT_GHOST_SIZE != 5 )
         Aux_Error( ERROR_INFO, "POT_GHOST_SIZE (%d) != 5 !!\n", POT_GHOST_SIZE );
#     endif
#  endif // if ( POT_SCHEME == SOR )


// (7) warp size
   if ( DeviceProp.warpSize != WARP_SIZE )
      Aux_Error( ERROR_INFO, "inconsistent warp size (warpSize %d, WARP_SIZE %d) !!\n",
                 DeviceProp.warpSize, WARP_SIZE );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "CUAPI_SetDevice ... done\n" );

} // FUNCTION : CUAPI_SetDevice



#endif // #ifdef GPU
