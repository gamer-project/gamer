#include "GAMER.h"

// function pointers
//void (*EoS_InitAuxArray_Ptr)( double AuxArray[] ) = NULL;
EoS_t CPUEoS_Ptr                                  = NULL;
void (*SetCPUEoS_Ptr)( EoS_t &CPUEoS_Ptr )        = NULL;
#ifdef GPU
EoS_t GPUEoS_Ptr                                  = NULL;
void (*SetGPUEoS_Ptr)( EoS_t &GPUEoS_Ptr )        = NULL;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Init
// Description :  Initialize the equation of state (EoS) routines and auxiliary arrays
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set the auxiliary CPU arrays by invoking EoS_InitAuxArray_Ptr()
//                   --> Corresponding GPU arrays are set by CUAPI_SetConstMemory()
//                3. Set the CPU/GPU EoS routines by invoking SetCPU/GPUEoS_Ptr()
//                4. For a user-specified EoS, all function pointers here must be set in advance by
//                   a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void EoS_Init()
{

// initialize the auxiliary CPU array
   /*
   if ( EoS_InitAuxArray_Ptr != NULL )
      EoS_InitAuxArray_Ptr( EoS_AuxArray );
   else
      Aux_Error( ERROR_INFO, "EoS_InitAuxArray_Ptr == NULL for external acceleration !!\n" );
   */

// set the CPU routine
   if ( SetCPUEoS_Ptr != NULL )
   {
      SetCPUEoS_Ptr( CPUEoS_Ptr );

      if ( CPUEoS_Ptr == NULL )
         Aux_Error( ERROR_INFO, "CPUEoS_Ptr == NULL for external acceleration !!\n" );
   }
   else
      Aux_Error( ERROR_INFO, "SetCPUEoS_Ptr == NULL for external acceleration !!\n" );

// set the GPU routine
#  ifdef GPU
   if ( SetGPUEoS_Ptr != NULL )
   {
      SetGPUEoS_Ptr( GPUEoS_Ptr );

      if ( GPUEoS_Ptr == NULL )
         Aux_Error( ERROR_INFO, "GPUEoS_Ptr == NULL for external acceleration !!\n" );
   }
   else
      Aux_Error( ERROR_INFO, "SetGPUEoS_Ptr == NULL !!\n" );
#  endif

} // FUNCTION : EoS_nit
