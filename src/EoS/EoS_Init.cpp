#include "GAMER.h"

#if ( MODEL == HYDRO )



// prototypes of built-in EoS
#if   ( EOS == EOS_GAMMA )
void EoS_Init_Gamma();
#elif ( EOS == EOS_ISOTHERMAL )
void EoS_Init_Isothermal();
#elif ( EOS == EOS_NUCLEAR )
# error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#endif // # EOS

// this function pointer must be set by a test problem initializer for non-built-in EoS
void (*EoS_Init_Ptr)() = NULL;


#ifdef GPU
extern real *d_EoS_Table[EOS_NTABLE_MAX];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  EoS_Init
// Description :  Initialize the equation of state (EoS)
//
// Note        :  1. Invoked by Init_GAMER()
//                2. For a non-built-in EoS, "EoS_Init_Ptr" must be set by a test problem initializer
//                   in advance
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void EoS_Init()
{

// check if EoS has been initialized already
// --> necessary since some test problem initializers may also call EoS_Init()
   static bool EoS_Initialized = false;

   if ( EoS_Initialized )  return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// initialize the EoS table pointers as NULL
   for (int t=0; t<EOS_NTABLE_MAX; t++)
   {
      h_EoS_Table[t] = NULL;
#     ifdef GPU
      d_EoS_Table[t] = NULL;
#     endif
   }


// set the initialization function pointer for the built-in EoS
#  if   ( EOS == EOS_GAMMA )
   EoS_Init_Ptr = EoS_Init_Gamma;
#  elif ( EOS == EOS_ISOTHERMAL )
   EoS_Init_Ptr = EoS_Init_Isothermal;
#  elif ( EOS == EOS_NUCLEAR )
#  error : ERROR : EOS_NUCLEAR is NOT supported yet !!
#  endif // # EOS


// initialize EoS
   if ( EoS_Init_Ptr != NULL )
      EoS_Init_Ptr();
   else
      Aux_Error( ERROR_INFO, "EoS_Init_Ptr == NULL for EoS %d !!\n", EOS );


// store relevant variables in the object "EoS" for the CPU/GPU solvers
#  ifdef GPU
   EoS.DensEint2Pres_FuncPtr = EoS_DensEint2Pres_GPUPtr;
   EoS.DensPres2Eint_FuncPtr = EoS_DensPres2Eint_GPUPtr;
   EoS.DensPres2CSqr_FuncPtr = EoS_DensPres2CSqr_GPUPtr;
   EoS.DensEint2Temp_FuncPtr = EoS_DensEint2Temp_GPUPtr;
   EoS.DensTemp2Pres_FuncPtr = EoS_DensTemp2Pres_GPUPtr;
   EoS.DensEint2Entr_FuncPtr = EoS_DensEint2Entr_GPUPtr;
   EoS.General_FuncPtr       = EoS_General_GPUPtr;

   CUAPI_SetConstMemory_EoS();

#  else

   EoS.DensEint2Pres_FuncPtr = EoS_DensEint2Pres_CPUPtr;
   EoS.DensPres2Eint_FuncPtr = EoS_DensPres2Eint_CPUPtr;
   EoS.DensPres2CSqr_FuncPtr = EoS_DensPres2CSqr_CPUPtr;
   EoS.DensEint2Temp_FuncPtr = EoS_DensEint2Temp_CPUPtr;
   EoS.DensTemp2Pres_FuncPtr = EoS_DensTemp2Pres_CPUPtr;
   EoS.DensEint2Entr_FuncPtr = EoS_DensEint2Entr_CPUPtr;
   EoS.General_FuncPtr       = EoS_General_CPUPtr;

   EoS.AuxArrayDevPtr_Flt    = EoS_AuxArray_Flt;
   EoS.AuxArrayDevPtr_Int    = EoS_AuxArray_Int;
   EoS.Table                 = h_EoS_Table;
#  endif // #ifdef GPU ... else ...


   EoS_Initialized = true;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : EoS_Init



#endif // #if ( MODEL == HYDRO )
