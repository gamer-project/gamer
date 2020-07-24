// define DEFINE_GLOBAL to declare all constant variables here
// --> must define it BEFORE including CUAPI.h since the latter will include "Macro.h" to set SET_GLOBAL()
#define DEFINE_GLOBAL
#include "CUAPI.h"
#include "CUDA_ConstMemory.h"
#undef DEFINE_GLOBAL

#ifdef GPU




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetConstMemory
// Description :  Set the constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by Init_GAMER()
//
// Parameter   :  None
//
// Return      :  None
//---------------------------------------------------------------------------------------------------
void CUAPI_SetConstMemory()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// cudaMemcpyToSymbol() has the default parameters "offset=0, kind=cudaMemcpyHostToDevice" and so
// we do not repeat them below

#if ( MODEL == HYDRO )

#  if ( NCOMP_PASSIVE > 0 )
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_NormIdx,         PassiveNorm_VarIdx, NCOMP_PASSIVE*sizeof(int) )  );
#  endif

#  ifdef GRAVITY
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_ExtAcc_AuxArray, ExtAcc_AuxArray,    EXT_ACC_NAUX_MAX*sizeof(double) )  );
#  endif

#elif ( MODEL == ELBDM )

#  ifdef GRAVITY
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_ExtPot_AuxArray, ExtPot_AuxArray,    EXT_POT_NAUX_MAX*sizeof(double) )  );
#  endif

#endif // MODEL


#  ifdef GRAVITY
   const real h_Mp[3] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
   const real h_Mm[3] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };

   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Mp,              h_Mp,               3*sizeof(real) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_Mm,              h_Mm,               3*sizeof(real) )  );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_SetConstMemory



#endif // #ifdef GPU
