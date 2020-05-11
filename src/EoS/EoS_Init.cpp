#include "GAMER.h"

#if ( MODEL == HYDRO )


// auxiliary array to be set by EoS_InitAuxArray_Ptr
double EoS_AuxArray[EOS_NAUX_MAX];

// function pointers
void (*EoS_InitAuxArray_Ptr)( double [] )           = NULL;
EoS_DE2P_t CPUEoS_DensEint2Pres_Ptr                 = NULL;
void (*SetCPUEoS_DensEint2Pres_Ptr)( EoS_DE2P_t & ) = NULL;
#ifdef GPU
EoS_DE2P_t GPUEoS_DensEint2Pres_Ptr                 = NULL;
void (*SetGPUEoS_DensEint2Pres_Ptr)( EoS_DE2P_t & ) = NULL;
#endif

// prototypes of built-in EoS
#if   ( EOS == EOS_IDEALGAS )
void EoS_InitAuxArray_IdealGas( double [] );
void SetCPUEoS_DensEint2Pres_IdealGas( EoS_DE2P_t & );
# ifdef GPU
void SetGPUEoS_DensEint2Pres_IdealGas( EoS_DE2P_t & );
# endif

#elif ( EOS == EOS_NUCLEAR )
# error : ERROR : EOS_NUCLEAR is NOT supported yet !!

#endif // # EOS




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Init
// Description :  Initialize the equation of state (EoS) routines and auxiliary arrays
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set the auxiliary CPU arrays by invoking EoS_InitAuxArray_*_Ptr()
//                   --> Corresponding GPU arrays are set by CUAPI_SetConstMemory()
//                3. Set the CPU/GPU EoS routines by invoking SetCPU/GPUEoS_*_Ptr()
//                4. For a tabular or user-specified EoS, all function pointers must be set in advance
//                   by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void EoS_Init()
{

// set function pointers for the built-in EoS
#  if   ( EOS == EOS_IDEALGAS )
   EoS_InitAuxArray_Ptr        = EoS_InitAuxArray_IdealGas;
   SetCPUEoS_DensEint2Pres_Ptr = SetCPUEoS_DensEint2Pres_IdealGas;
#  ifdef GPU
   SetGPUEoS_DensEint2Pres_Ptr = SetGPUEoS_DensEint2Pres_IdealGas;
#  endif

#  elif ( EOS == EOS_NUCLEAR )
#  error : ERROR : EOS_NUCLEAR is NOT supported yet !!

#  endif // # EOS


// set the auxiliary CPU array
   if ( EoS_InitAuxArray_Ptr != NULL )
      EoS_InitAuxArray_Ptr( EoS_AuxArray );
   else
      Aux_Error( ERROR_INFO, "EoS_InitAuxArray_Ptr == NULL for EoS %d !!\n", EOS );


// set the CPU routines
   if ( SetCPUEoS_DensEint2Pres_Ptr != NULL )
   {
      SetCPUEoS_DensEint2Pres_Ptr( CPUEoS_DensEint2Pres_Ptr );

      if ( CPUEoS_DensEint2Pres_Ptr == NULL )
         Aux_Error( ERROR_INFO, "CPUEoS_DensEint2Pres_Ptr == NULL for EoS %d !!\n", EOS );
   }
   else
      Aux_Error( ERROR_INFO, "SetCPUEoS_DensEint2Pres_Ptr == NULL for EoS %d !!\n", EOS );


// set the GPU routines
#  ifdef GPU
   if ( SetGPUEoS_DensEint2Pres_Ptr != NULL )
   {
      SetGPUEoS_DensEint2Pres_Ptr( GPUEoS_DensEint2Pres_Ptr );

      if ( GPUEoS_DensEint2Pres_Ptr == NULL )
         Aux_Error( ERROR_INFO, "GPUEoS_DensEint2Pres_Ptr == NULL for EoS %d !!\n", EOS );
   }
   else
      Aux_Error( ERROR_INFO, "SetGPUEoS_DensEint2Pres_Ptr == NULL for EoS %d !!\n", EOS );
#  endif

} // FUNCTION : EoS_Init



#endif // #if ( MODEL == HYDRO )
