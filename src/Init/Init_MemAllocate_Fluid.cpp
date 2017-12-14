#ifndef GPU

#include "GAMER.h"



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Fluid
// Description :  Allocate memory for the fluid solver
//
// Note        :  Work when using CPUs only
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup, const int Pot_NPatchGroup )
{

   const int Flu_NPatch = 8*Flu_NPatchGroup;
#  ifdef GRAVITY
   const int Pot_NPatch = 8*Pot_NPatchGroup;
   const int dt_NPatch  = MAX( Flu_NPatch, Pot_NPatch );
#  else
   const int dt_NPatch  = Flu_NPatch;
#  endif

   for (int t=0; t<2; t++)
   {
      h_Flu_Array_F_In [t] = new real [Flu_NPatchGroup][FLU_NIN ][   FLU_NXT   *FLU_NXT   *FLU_NXT    ];
      h_Flu_Array_F_Out[t] = new real [Flu_NPatchGroup][FLU_NOUT][ 8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ];

      if ( amr->WithFlux )
      h_Flux_Array     [t] = new real [Flu_NPatchGroup][9][NFLUX_TOTAL][ 4*PATCH_SIZE*PATCH_SIZE ];

#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_F[t] = new real [Flu_NPatchGroup][USG_NXT_F][USG_NXT_F][USG_NXT_F];

      if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
      h_Corner_Array_F [t] = new double [Flu_NPatchGroup][3];
#     endif

      h_dt_Array_T     [t] = new real [dt_NPatch];
      h_Flu_Array_T    [t] = new real [Flu_NPatch][NCOMP_FLUID][ CUBE(PS1) ];

#     ifdef DUAL_ENERGY
      h_DE_Array_F_Out [t] = new char [Flu_NPatchGroup][ 8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ];
#     endif
   }

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
