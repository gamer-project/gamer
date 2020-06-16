#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtAccPot
// Description :  Initialize the external acceleration and potential routines and auxiliary arrays
//
// Note        :  1. Invoked by Init_GAMER() and EvolveLevel()
//                2. Enabled by the runtime options "OPT__GRAVITY_TYPE == 2/3" and "OPT__EXTERNAL_POT"
//                3. Set the auxiliary CPU arrays of the external acceleration and potential by invoking
//
//                      Init_ExtAccAuxArray_Ptr() and Init_ExtPotAuxArray_Ptr()
//
//                   --> Corresponding GPU arrays are set by CUAPI_SetConstMemory()
//                4. Set the CPU/GPU external acceleration and potential routines by invoking
//
//                      SetCPU/GPUExtAcc_Ptr() and SetCPU/GPUExtPot_Ptr()
//
//                5. Function pointers used here must be set in advance by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_ExtAccPot()
{

// initialize the external acceleration
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
//    initialize the auxiliary CPU array
      if ( Init_ExtAccAuxArray_Ptr != NULL )
         Init_ExtAccAuxArray_Ptr( ExtAcc_AuxArray );
      else
         Aux_Error( ERROR_INFO, "Init_ExtAccAuxArray_Ptr == NULL for external acceleration !!\n" );

//    set the CPU routine
      if ( SetCPUExtAcc_Ptr != NULL )
      {
         SetCPUExtAcc_Ptr( CPUExtAcc_Ptr );

         if ( CPUExtAcc_Ptr == NULL )
            Aux_Error( ERROR_INFO, "CPUExtAcc_Ptr == NULL for external acceleration !!\n" );
      }
      else
         Aux_Error( ERROR_INFO, "SetCPUExtAcc_Ptr == NULL for external acceleration !!\n" );

//    set the GPU routine
#     ifdef GPU
      if ( SetGPUExtAcc_Ptr != NULL )
      {
         SetGPUExtAcc_Ptr( GPUExtAcc_Ptr );

         if ( GPUExtAcc_Ptr == NULL )
            Aux_Error( ERROR_INFO, "GPUExtAcc_Ptr == NULL for external acceleration !!\n" );
      }
      else
         Aux_Error( ERROR_INFO, "SetGPUExtAcc_Ptr == NULL !!\n" );
#     endif
   } // if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )


// initialize the external potential
   if ( OPT__EXTERNAL_POT )
   {
//    initialize the auxiliary CPU array
      if ( Init_ExtPotAuxArray_Ptr != NULL )
         Init_ExtPotAuxArray_Ptr( ExtPot_AuxArray );
      else
         Aux_Error( ERROR_INFO, "Init_ExtPotAuxArray_Ptr == NULL for external potential !!\n" );

//    set the CPU routine
      if ( SetCPUExtPot_Ptr != NULL )
      {
         SetCPUExtPot_Ptr( CPUExtPot_Ptr );

         if ( CPUExtPot_Ptr == NULL )
            Aux_Error( ERROR_INFO, "CPUExtPot_Ptr == NULL for external potential !!\n" );
      }
      else
         Aux_Error( ERROR_INFO, "SetCPUExtPot_Ptr == NULL for external potential !!\n" );

//    set the GPU routine
#     ifdef GPU
      if ( SetGPUExtPot_Ptr != NULL )
      {
         SetGPUExtPot_Ptr( GPUExtPot_Ptr );

         if ( GPUExtPot_Ptr == NULL )
            Aux_Error( ERROR_INFO, "GPUExtPot_Ptr == NULL for external potential !!\n" );
      }
      else
         Aux_Error( ERROR_INFO, "SetGPUExtPot_Ptr == NULL !!\n" );
#     endif
   } // if ( OPT__EXTERNAL_POT )

} // FUNCTION : Init_ExtAccPot



#endif // #ifdef GRAVITY
