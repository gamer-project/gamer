#include "CUAPI.h"

#ifdef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetMemSize
// Description :  Set parameters controlling the size of GPU global memory allocation
//                --> GPU_NSTREAM, *_GPU_NPGROUP
//
// Parameter   :  GPU_NStream     : Number of streams for the asynchronous memory copy in GPU
//                Flu_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the fluid solver
//                Pot_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the Poisson solver
//                Che_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the Grackle solver
//                Src_GPU_NPGroup : Number of patch groups sent into GPU simultaneously for the source-term solver
//-------------------------------------------------------------------------------------------------------
void CUAPI_SetMemSize( int &GPU_NStream, int &Flu_GPU_NPGroup, int &Pot_GPU_NPGroup, int &Che_GPU_NPGroup,
                       int &Src_GPU_NPGroup )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// get the device ID
   int GetDeviceID = 999;
   CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );


// load the device properties
   cudaDeviceProp DeviceProp;
   CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );


// (1) GPU_NSTREAM
   if ( GPU_NStream <= 0 )
   {
      if ( DeviceProp.deviceOverlap )
      {
#        if   ( MODEL == HYDRO )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == MAXWELL )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == PASCAL )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == VOLTA )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == TURING )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == AMPERE )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == ADA_LOVELACE )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == HOPPER )
            GPU_NStream = 4;
#           else
#           error : UNKNOWN GPU_ARCH !!
#           endif

#        elif ( MODEL == ELBDM )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == MAXWELL )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == PASCAL )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == VOLTA )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == TURING )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == AMPERE )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == ADA_LOVELACE )
            GPU_NStream = 4;
#           elif ( GPU_ARCH == HOPPER )
            GPU_NStream = 4;
#           else
#           error : ERROR : UNKNOWN GPU_ARCH !!
#           endif
#        else
#           error : ERROR : UNKNOWN MODEL !!
#        endif // MODEL
      } // if ( DeviceProp.deviceOverlap )

      else
         GPU_NStream = 1;

      PRINT_RESET_PARA( GPU_NStream, FORMAT_INT, "--> might be further fine-tuned (GPU_NSTREAM in Input__Parameter)" );
   } // if ( GPU_NStream <= 0 )


// (2) XXX_GPU_NPGROUP
// (2-1) FLU_GPU_NPGROUP
   if ( Flu_GPU_NPGroup <= 0 )
   {
#     if   ( MODEL == HYDRO )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == MAXWELL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == PASCAL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == VOLTA )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == TURING )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == AMPERE )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == ADA_LOVELACE )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == HOPPER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif

#     elif ( MODEL == ELBDM )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == MAXWELL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == PASCAL )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == VOLTA )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == TURING )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == AMPERE )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == ADA_LOVELACE )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == HOPPER )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif
#     else
#        error : ERROR : UNKNOWN MODEL !!
#     endif // MODEL

      PRINT_RESET_PARA( Flu_GPU_NPGroup, FORMAT_INT, "--> might be further fine-tuned (FLU_GPU_NPGROUP in Input__Parameter)" );
   } // if ( Flu_GPU_NPGroup <= 0 )

// (2-2) POT_GPU_NPGROUP
#  ifdef GRAVITY
   if ( Pot_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == AMPERE )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == ADA_LOVELACE )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == HOPPER )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      PRINT_RESET_PARA( Pot_GPU_NPGroup, FORMAT_INT, "--> might be further fine-tuned (POT_GPU_NPGROUP in Input__Parameter)" );
   } // if ( Pot_GPU_NPGroup <= 0 )
#  endif

// (2-3) CHE_GPU_NPGROUP
#  ifdef SUPPORT_GRACKLE
   if ( Che_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == AMPERE )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == ADA_LOVELACE )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == HOPPER )
      Che_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      PRINT_RESET_PARA( Che_GPU_NPGroup, FORMAT_INT, "--> might be further fine-tuned (CHE_GPU_NPGROUP in Input__Parameter)" );
   } // if ( Che_GPU_NPGroup <= 0 )
#  endif

// (2-4) SRC_GPU_NPGROUP
   if ( Src_GPU_NPGroup <= 0 )
   {
#     if   ( GPU_ARCH == FERMI )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == MAXWELL )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == PASCAL )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == VOLTA )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == TURING )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == AMPERE )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == ADA_LOVELACE )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == HOPPER )
      Src_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      PRINT_RESET_PARA( Src_GPU_NPGroup, FORMAT_INT, "--> might be further fine-tuned (SRC_GPU_NPGROUP in Input__Parameter)" );
   } // if ( Src_GPU_NPGroup <= 0 )


// check
   if ( GPU_NStream < 1 )  Aux_Error( ERROR_INFO, "GPU_NSTREAM (%d) < 1 !!\n", GPU_NStream );

   if ( Flu_GPU_NPGroup % GPU_NStream != 0 )
      Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 Flu_GPU_NPGroup, GPU_NStream );

#  ifdef GRAVITY
   if ( Pot_GPU_NPGroup % GPU_NStream != 0 )
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 Pot_GPU_NPGroup, GPU_NStream );
#  endif

#  ifdef SUPPORT_GRACKLE
   /*
   if ( Che_GPU_NPGroup % GPU_NStream != 0 )
      Aux_Error( ERROR_INFO, "CHE_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 Che_GPU_NPGroup, GPU_NStream );
                 */
#  endif

   if ( Src_GPU_NPGroup % GPU_NStream != 0 )
      Aux_Error( ERROR_INFO, "SRC_GPU_NPGROUP (%d) %% GPU_NSTREAM (%d) != 0 !!\n",
                 Src_GPU_NPGroup, GPU_NStream );

#  ifdef OPENMP
   if ( Flu_GPU_NPGroup < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "FLU_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", Flu_GPU_NPGroup, OMP_NTHREAD );

#  ifdef GRAVITY
   if ( Pot_GPU_NPGroup < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "POT_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", Pot_GPU_NPGroup, OMP_NTHREAD );
#  endif

#  ifdef SUPPORT_GRACKLE
   if ( Che_GPU_NPGroup < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "CHE_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", Che_GPU_NPGroup, OMP_NTHREAD );
#  endif

   if ( Src_GPU_NPGroup < OMP_NTHREAD )
      Aux_Error( ERROR_INFO, "SRC_GPU_NPGROUP (%d) < OMP_NTHREAD (%d) !!\n", Src_GPU_NPGroup, OMP_NTHREAD );
#  endif // #ifdef OPENMP


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_SetMemSize



#endif // #ifdef GPU
