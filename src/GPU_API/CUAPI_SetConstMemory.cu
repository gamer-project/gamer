// define DEFINE_GLOBAL to declare all constant variables here
// --> must define it BEFORE including CUAPI.h since the latter will include "Macro.h" to set SET_GLOBAL()
#define DEFINE_GLOBAL
#include "CUAPI.h"
#include "CUDA_ConstMemory.h"
#undef DEFINE_GLOBAL

#ifdef GPU

extern real *d_EoS_Table[EOS_NTABLE_MAX];

#ifdef GRAVITY
void CUAPI_SetConstMemory_ExtAccPot();
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetConstMemory
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Init_GAMER()
//                3. Invoke CUAPI_SetConstMemory_ExtAccPot()
//
// Parameter   :  None
//
// Return      :  c_NormIdx[], c_Mp[], c_Mm[], c_ExtAcc_AuxArray[], c_ExtPot_AuxArray[]
//---------------------------------------------------------------------------------------------------
void CUAPI_SetConstMemory()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// cudaMemcpyToSymbol() has the default parameters "offset=0, kind=cudaMemcpyHostToDevice" and so
// we do not repeat them below

#  if ( MODEL == HYDRO )
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_AuxArray_Flt, EoS_AuxArray_Flt, EOS_NAUX_MAX  *sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_AuxArray_Int, EoS_AuxArray_Int, EOS_NAUX_MAX  *sizeof(int   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_Table,        d_EoS_Table,      EOS_NTABLE_MAX*sizeof(real* ) )  );
#  endif

#  if ( NCOMP_PASSIVE > 0 )
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_NormIdx, PassiveNorm_VarIdx, NCOMP_PASSIVE*sizeof(int) )  );
#  endif

#  ifdef GRAVITY
   CUAPI_SetConstMemory_ExtAccPot();

   const real h_Mp[3] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
   const real h_Mm[3] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };

   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Mp, h_Mp, 3*sizeof(real) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Mm, h_Mm, 3*sizeof(real) )  );
#  endif // #ifdef GRAVITY


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_SetConstMemory



#endif // #ifdef GPU
