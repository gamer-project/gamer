#ifndef GPU

#include "GAMER.h"
#include "CUFLU.h"


#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*h_PriVar)      [NCOMP_TOTAL][ CUBE(FLU_NXT)     ];
extern real (*h_Slope_PPM)[3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ];
extern real (*h_FC_Var)   [6][NCOMP_TOTAL][ CUBE(N_FC_VAR)    ];
extern real (*h_FC_Flux)  [3][NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ];
#endif




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
      h_Flu_Array_F_In [t] = new real [Flu_NPatchGroup][FLU_NIN ][ CUBE(FLU_NXT) ];
      h_Flu_Array_F_Out[t] = new real [Flu_NPatchGroup][FLU_NOUT][ CUBE(PS2) ];

      if ( amr->WithFlux )
      h_Flux_Array     [t] = new real [Flu_NPatchGroup][9][NFLUX_TOTAL][ SQR(PS2) ];

#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_F[t] = new real [Flu_NPatchGroup][ CUBE(USG_NXT_F) ];

      if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
      h_Corner_Array_F [t] = new double [Flu_NPatchGroup][3];
#     endif

      h_dt_Array_T     [t] = new real [dt_NPatch];
      h_Flu_Array_T    [t] = new real [Flu_NPatch][NCOMP_FLUID][ CUBE(PS1) ];

#     ifdef DUAL_ENERGY
      h_DE_Array_F_Out [t] = new char [Flu_NPatchGroup][ CUBE(PS2) ];
#     endif
   } // for (int t=0; t<2; t++)


#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   h_FC_Var    = new real [Flu_NPatchGroup][6][NCOMP_TOTAL][ CUBE(N_FC_VAR)    ];
   h_FC_Flux   = new real [Flu_NPatchGroup][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ];
   h_PriVar    = new real [Flu_NPatchGroup]   [NCOMP_TOTAL][ CUBE(FLU_NXT)     ];
#  if ( LR_SCHEME == PPM )
   h_Slope_PPM = new real [Flu_NPatchGroup][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ];
#  endif
#  endif // FLU_SCHEME

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
