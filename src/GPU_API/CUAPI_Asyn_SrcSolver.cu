#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



__global__
void CUSRC_SrcSolver_IterateAllCells(
   const real g_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
         real g_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
   const real g_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
   const double g_Corner_Array[][3],
   const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
   const double TimeNew, const double TimeOld,
   const real MinDens, const real MinPres, const real MinEint, const EoS_t EoS );

// device pointers
extern real (*d_Flu_Array_S_In )[FLU_NIN_S ][ CUBE(SRC_NXT)           ];
extern real (*d_Flu_Array_S_Out)[FLU_NOUT_S][ CUBE(PS1)               ];
#ifdef MHD
extern real (*d_Mag_Array_S_In)[NCOMP_MAG  ][ SRC_NXT_P1*SQR(SRC_NXT) ];
#else
static real (*d_Mag_Array_S_In)[NCOMP_MAG  ][ SRC_NXT_P1*SQR(SRC_NXT) ] = NULL;
#endif
extern double (*d_Corner_Array_S)[3];

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_SrcSolver
// Description :  Invoke the source-term solvers
//
//                ***********************************************************
//                **                Asynchronous Function                  **
//                **                                                       **
//                **  will return before the execution in GPU is complete  **
//                ***********************************************************
//
// Note        :  1. Use streams for the asychronous memory copy between device and host
//                2. Prefix "d" : for pointers pointing to the "Device" memory space
//                   Prefix "h" : for pointers pointing to the "Host"   memory space
//
// Parameter   :  h_Flu_Array_In    : Host array storing the input fluid variables
//                h_Flu_Array_Out   : Host array to store the output fluid variables
//                h_Mag_Array_In    : Host array storing the input B field (for MHD only)
//                h_Corner_Array    : Host array storing the physical corner coordinates of each patch
//                SrcTerms          : Structure storing all source-term variables
//                NPatchGroup       : Number of patch groups to be evaluated
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                GPU_NStream       : Number of CUDA streams for the asynchronous memory copy
//
// Return      :  h_Flu_Array_Out[]
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_SrcSolver( const real h_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
                                 real h_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
                           const real h_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
                           const double h_Corner_Array[][3],
                           const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
                           const double TimeNew, const double TimeOld,
                           const real MinDens, const real MinPres, const real MinEint,
                           const int GPU_NStream )
{

// check
#  ifdef GAMER_DEBUG
   if ( h_Flu_Array_In  == NULL )   Aux_Error( ERROR_INFO, "h_Flu_Array_In = NULL !!\n" );
   if ( h_Flu_Array_Out == NULL )   Aux_Error( ERROR_INFO, "h_Flu_Array_Out = NULL !!\n" );
#  ifdef MHD
   if ( h_Mag_Array_In  == NULL )   Aux_Error( ERROR_INFO, "h_Mag_Array_In = NULL !!\n" );
#  endif
   if ( h_Corner_Array  == NULL )   Aux_Error( ERROR_INFO, "h_Corner_Array = NULL !!\n" );
#  endif


// EoS is not defined when MODEL != HYDRO
#  if ( MODEL != HYDRO )
   EoS_t EoS;
#  endif


// set the block size
   dim3 BlockDim_SrcSolver( SRC_BLOCK_SIZE, 1, 1 );


// set the number of patches and the corresponding data size to be transferred into GPU in each stream
   const int NPatch = NPatchGroup*8;

   int *NPatch_per_Stream  = new int [GPU_NStream];
   int *UsedPatch          = new int [GPU_NStream];
   int *Flu_MemSize_In     = new int [GPU_NStream];
   int *Flu_MemSize_Out    = new int [GPU_NStream];
#  ifdef MHD
   int *Mag_MemSize_In     = new int [GPU_NStream];
#  endif
   int *Corner_MemSize     = new int [GPU_NStream];

// number of patches in each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatch;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatch / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatch - UsedPatch[GPU_NStream-1];
   }

// corresponding data size to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Flu_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN_S *CUBE(SRC_NXT);
      Flu_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NOUT_S*CUBE(PS1);
#     ifdef MHD
      Mag_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_MAG*SRC_NXT_P1*SQR(SRC_NXT);
#     endif
      Corner_MemSize [s] = sizeof(double)*NPatch_per_Stream[s]*3;
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_S_In + UsedPatch[s], h_Flu_Array_In + UsedPatch[s],
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );

#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Mag_Array_S_In + UsedPatch[s], h_Mag_Array_In + UsedPatch[s],
                         Mag_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_S + UsedPatch[s], h_Corner_Array + UsedPatch[s],
                         Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUSRC_SrcSolver_IterateAllCells <<< NPatch_per_Stream[s], BlockDim_SrcSolver, 0, Stream[s] >>>
                                      ( d_Flu_Array_S_In  + UsedPatch[s],
                                        d_Flu_Array_S_Out + UsedPatch[s],
                                        d_Mag_Array_S_In  + UsedPatch[s],
                                        d_Corner_Array_S  + UsedPatch[s],
                                        SrcTerms, NPatchGroup, dt, dh, TimeNew, TimeOld,
                                        MinDens, MinPres, MinEint, EoS );

      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array_Out + UsedPatch[s], d_Flu_Array_S_Out + UsedPatch[s],
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] Flu_MemSize_In;
   delete [] Flu_MemSize_Out;
#  ifdef MHD
   delete [] Mag_MemSize_In;
#  endif
   delete [] Corner_MemSize;

} // FUNCTION : CUAPI_Asyn_SrcSolver



#endif // #ifdef GPU
