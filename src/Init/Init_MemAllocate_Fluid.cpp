#include "Copyright.h"
#ifndef GPU

#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Fluid
// Description :  Allocate memory for the fluid solver
//
// Note        :  Work when using CPUs only
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup )
{

   for (int t=0; t<2; t++)
   {
      h_Flu_Array_F_In       [t] = new real [Flu_NPatchGroup][FLU_NIN ][  FLU_NXT   *FLU_NXT   *FLU_NXT   ];
      h_Flu_Array_F_Out      [t] = new real [Flu_NPatchGroup][FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE];

      if ( amr->WithFlux )
      h_Flux_Array           [t] = new real [Flu_NPatchGroup][9][NFLUX][4*PATCH_SIZE*PATCH_SIZE];

#     ifdef UNSPLIT_GRAVITY
      if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
      h_Corner_Array_F       [t] = new double [Flu_NPatchGroup][3];
#     endif

      if ( OPT__ADAPTIVE_DT )
      h_MinDtInfo_Fluid_Array[t] = new real [Flu_NPatchGroup];
   }

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
