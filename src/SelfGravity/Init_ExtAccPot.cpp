#include "GAMER.h"

#ifdef GRAVITY



#include "CUPOT.h"
// auxiliary arrays to be set by Init_ExtAcc/PotAuxArray_Ptr
double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
double ExtPot_AuxArray[EXT_POT_NAUX_MAX];


// function pointers --> must be set by a test problem initializer
void (*Init_ExtAccAuxArray_Ptr)( double AuxArray[] ) = NULL;
ExtAcc_t CPUExtAcc_Ptr                               = NULL;
void (*SetCPUExtAcc_Ptr)( ExtAcc_t &CPUExtAcc_Ptr )  = NULL;
#ifdef GPU
ExtAcc_t GPUExtAcc_Ptr                               = NULL;
void (*SetGPUExtAcc_Ptr)( ExtAcc_t &GPUExtAcc_Ptr )  = NULL;
#endif

void (*Init_ExtPotAuxArray_Ptr)( double AuxArray[] ) = NULL;
ExtPot_t CPUExtPot_Ptr                               = NULL;
void (*SetCPUExtPot_Ptr)( ExtPot_t &CPUExtPot_Ptr )  = NULL;
#ifdef GPU
ExtPot_t GPUExtPot_Ptr                               = NULL;
void (*SetGPUExtPot_Ptr)( ExtPot_t &GPUExtPot_Ptr )  = NULL;
#endif




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
//                      SetCPU/GPUExtAcc_PointMass() and SetCPU/GPUExtPot_PointMass()
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
      if ( Init_ExtAccAuxArray_Ptr == NULL )   Aux_Error( ERROR_INFO, "Init_ExtAccAuxArray_Ptr == NULL !!\n" );
      Init_ExtAccAuxArray_Ptr( ExtAcc_AuxArray );

//    set the CPU routine
      if ( SetCPUExtAcc_Ptr == NULL )   Aux_Error( ERROR_INFO, "SetCPUExtAcc_Ptr == NULL !!\n" );
      SetCPUExtAcc_Ptr( CPUExtAcc_Ptr );

//    set the GPU routine
#     ifdef GPU
      if ( SetGPUExtAcc_Ptr == NULL )   Aux_Error( ERROR_INFO, "SetGPUExtAcc_Ptr == NULL !!\n" );
      SetGPUExtAcc_Ptr( GPUExtAcc_Ptr );
#     endif
   }


// initialize the external potential
   if ( OPT__EXTERNAL_POT )
   {
//    initialize the auxiliary CPU array
      if ( Init_ExtPotAuxArray_Ptr == NULL )   Aux_Error( ERROR_INFO, "Init_ExtPotAuxArray_Ptr == NULL !!\n" );
      Init_ExtPotAuxArray_Ptr( ExtPot_AuxArray );

//    set the CPU routine
      if ( SetCPUExtPot_Ptr == NULL )   Aux_Error( ERROR_INFO, "SetCPUExtPot_Ptr == NULL !!\n" );
      SetCPUExtPot_Ptr( CPUExtPot_Ptr );

//    set the GPU routine
#     ifdef GPU
      if ( SetGPUExtPot_Ptr == NULL )   Aux_Error( ERROR_INFO, "SetGPUExtPot_Ptr == NULL !!\n" );
      SetGPUExtPot_Ptr( GPUExtPot_Ptr );
#     endif
   }

} // FUNCTION : Init_ExtAccPot



#endif // #ifdef GRAVITY
