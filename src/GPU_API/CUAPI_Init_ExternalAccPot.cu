#include "CUAPI.h"

#if ( defined GPU  &&  defined GRAVITY )


#include "CUPOT.h"

extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];


// declare all GPU kernels requiring ExtPot_AuxArray[] and/or ExtAcc_AuxArray[]
#if ( MODEL == HYDRO )
int CUPOT_HydroGravitySolver_SetConstMem( double ExtAcc_AuxArray_h[] );
int CUPOT_dtSolver_HydroGravity_SetConstMem( double ExtAcc_AuxArray_h[] );
#if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
int CUFLU_FluidSolver_SetConstMem_ExtAcc( double ExtAcc_AuxArray_h[] );
#endif
#endif // #if ( MODEL == HYDRO )

#if ( MODEL == ELBDM )
int CUPOT_ELBDMGravitySolver_SetConstMem( double ExtPot_AuxArray_h[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Init_ExternalAccPot
// Description :  Set the auxiliary GPU constant-memory arrays for the external acceleration and potential

// Note        :  1. Invoked by Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  ExtPot_AuxArray[] and ExtAcc_AuxArray[] in various CUDA kernels
//-------------------------------------------------------------------------------------------------------
void CUAPI_Init_ExternalAccPot()
{

#  if ( MODEL == HYDRO )
   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
   {
      if (  CUPOT_HydroGravitySolver_SetConstMem( ExtAcc_AuxArray ) != 0  )
         Aux_Error( ERROR_INFO, "CUPOT_HydroGravitySolver_SetConstMem failed ...\n" );

      if (  CUPOT_dtSolver_HydroGravity_SetConstMem( ExtAcc_AuxArray ) != 0  )
         Aux_Error( ERROR_INFO, "CUPOT_dtSolver_HydroGravity_SetConstMem failed ...\n" );

#     if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
      if (  CUFLU_FluidSolver_SetConstMem_ExtAcc( ExtAcc_AuxArray ) != 0  )
         Aux_Error( ERROR_INFO, "CUFLU_FluidSolver_SetConstMem_ExtAcc failed ...\n" );
#     endif
   }
#  endif // if ( MODEL == HYDRO )

#  if ( MODEL == ELBDM )
   if (  OPT__EXTERNAL_POT  &&  CUPOT_ELBDMGravitySolver_SetConstMem( ExtPot_AuxArray ) != 0  )
      Aux_Error( ERROR_INFO, "CUPOT_ELBDMGravitySolver_SetConstMem failed ...\n" );
#  endif

} // FUNCTION : CUAPI_Init_ExternalAccPot



#endif // #if ( defined GPU  &&  defined GRAVITY )
