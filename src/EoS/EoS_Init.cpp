#include "GAMER.h"

#if ( MODEL == HYDRO )



// prototypes of built-in EoS
#if   ( EOS == EOS_IDEALGAS )
void EoS_InitAuxArray_IdealGas( double [] );
void EoS_InitCPUFunc_IdealGas( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
# ifdef GPU
void EoS_InitGPUFunc_IdealGas( EoS_DE2P_t &, EoS_DP2E_t &, EoS_DP2C_t & );
# endif

#elif ( EOS == EOS_NUCLEAR )
# error : ERROR : EOS_NUCLEAR is NOT supported yet !!

#endif // # EOS




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Init
// Description :  Initialize the equation of state (EoS) routines and auxiliary arrays
//
// Note        :  1. Invoked by Init_GAMER()
//                2. Set the auxiliary CPU arrays by invoking EoS_InitAuxArray_Ptr()
//                   --> Corresponding GPU arrays are set by CUAPI_SetConstMemory()
//                3. Set the CPU/GPU EoS routines by invoking EoS_InitCPU/GPUFunc_Ptr()
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
   EoS_InitAuxArray_Ptr = EoS_InitAuxArray_IdealGas;
   EoS_InitCPUFunc_Ptr  = EoS_InitCPUFunc_IdealGas;
#  ifdef GPU
   EoS_InitGPUFunc_Ptr  = EoS_InitGPUFunc_IdealGas;
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
   if ( EoS_InitCPUFunc_Ptr != NULL )
   {
      EoS_InitCPUFunc_Ptr( EoS_DensEint2Pres_CPUPtr, EoS_DensPres2Eint_CPUPtr, EoS_DensPres2CSqr_CPUPtr );

      if ( EoS_DensEint2Pres_CPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensEint2Pres_CPUPtr == NULL for EoS %d !!\n", EOS );

      if ( EoS_DensPres2Eint_CPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensPres2Eint_CPUPtr == NULL for EoS %d !!\n", EOS );

      if ( EoS_DensPres2CSqr_CPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensPres2CSqr_CPUPtr == NULL for EoS %d !!\n", EOS );
   }
   else
      Aux_Error( ERROR_INFO, "EoS_InitCPUFunc_Ptr == NULL for EoS %d !!\n", EOS );


// set the GPU routines
#  ifdef GPU
   if ( EoS_InitGPUFunc_Ptr != NULL )
   {
      EoS_InitGPUFunc_Ptr( EoS_DensEint2Pres_GPUPtr, EoS_DensPres2Eint_GPUPtr, EoS_DensPres2CSqr_GPUPtr );

      if ( EoS_DensEint2Pres_GPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensEint2Pres_GPUPtr == NULL for EoS %d !!\n", EOS );

      if ( EoS_DensPres2Eint_GPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensPres2Eint_GPUPtr == NULL for EoS %d !!\n", EOS );

      if ( EoS_DensPres2CSqr_GPUPtr == NULL )
         Aux_Error( ERROR_INFO, "EoS_DensPres2CSqr_GPUPtr == NULL for EoS %d !!\n", EOS );
   }
   else
      Aux_Error( ERROR_INFO, "EoS_InitGPUFunc_Ptr == NULL for EoS %d !!\n", EOS );
#  endif // #ifdef GPU

} // FUNCTION : EoS_Init



#endif // #if ( MODEL == HYDRO )
