#ifndef GPU

#include "GAMER.h"
#include "CUFLU.h"


#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*h_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*h_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*h_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*h_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*h_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*h_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#endif
#endif // FLU_SCHEME




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Fluid
// Description :  Allocate memory for the fluid solver
//
// Note        :  Work when using CPUs only
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Fluid( const int Flu_NPatchGroup, const int Pot_NPatchGroup, const int Src_NPatchGroup )
{

   const int Flu_NPatch = 8*Flu_NPatchGroup;
#  ifdef GRAVITY
   const int Pot_NPatch = 8*Pot_NPatchGroup;
   const int dt_NPatch  = MAX( Flu_NPatch, Pot_NPatch );
#  else
   const int dt_NPatch  = Flu_NPatch;
#  endif
   const int Src_NPatch = 8*Src_NPatchGroup;

   for (int t=0; t<2; t++)
   {
      h_Flu_Array_F_In     [t] = new real [Flu_NPatchGroup][FLU_NIN ][ CUBE(FLU_NXT) ];
      h_Flu_Array_F_Out    [t] = new real [Flu_NPatchGroup][FLU_NOUT][ CUBE(PS2) ];

      if ( amr->WithFlux )
      h_Flux_Array         [t] = new real [Flu_NPatchGroup][9][NFLUX_TOTAL][ SQR(PS2) ];

#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_F    [t] = new real [Flu_NPatchGroup][ CUBE(USG_NXT_F) ];

      if ( OPT__EXT_ACC     )
      h_Corner_Array_F     [t] = new double [Flu_NPatchGroup][3];
#     endif

      h_dt_Array_T         [t] = new real [dt_NPatch];
      h_Flu_Array_T        [t] = new real [Flu_NPatch][FLU_NIN_T][ CUBE(PS1) ];

#     ifdef DUAL_ENERGY
      h_DE_Array_F_Out     [t] = new char [Flu_NPatchGroup][ CUBE(PS2) ];
#     endif

#     ifdef MHD
      h_Mag_Array_F_In     [t] = new real [Flu_NPatchGroup][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
      h_Mag_Array_F_Out    [t] = new real [Flu_NPatchGroup][NCOMP_MAG][ PS2P1*SQR(PS2) ];

      if ( amr->WithElectric )
      h_Ele_Array          [t] = new real [Flu_NPatchGroup][9][NCOMP_ELE][ PS2P1*PS2 ];

      h_Mag_Array_T        [t] = new real [Flu_NPatch][NCOMP_MAG][ PS1P1*SQR(PS1) ];
#     endif

      if ( SrcTerms.Any ) {
      h_Flu_Array_S_In     [t] = new real [Src_NPatch][FLU_NIN_S ][ CUBE(SRC_NXT)           ];
      h_Flu_Array_S_Out    [t] = new real [Src_NPatch][FLU_NOUT_S][ CUBE(PS1)               ];
#     ifdef MHD
      h_Mag_Array_S_In     [t] = new real [Src_NPatch][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ];
#     endif
      h_Corner_Array_S     [t] = new double [Src_NPatch][3];
      }

#     if ( MODEL == ELBDM )
      h_IsCompletelyRefined[t] = new bool [Flu_NPatchGroup];
#     endif

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      h_HasWaveCounterpart [t] = new bool [Flu_NPatchGroup][ CUBE(HYB_NXT) ];
#     endif
   } // for (int t=0; t<2; t++)


#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   h_FC_Var         = new real [Flu_NPatchGroup][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
   h_FC_Flux        = new real [Flu_NPatchGroup][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
   h_PriVar         = new real [Flu_NPatchGroup]   [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
#  if ( LR_SCHEME == PPM )
   h_Slope_PPM      = new real [Flu_NPatchGroup][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
#  endif
#  ifdef MHD
   h_FC_Mag_Half    = new real [Flu_NPatchGroup][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
   h_EC_Ele         = new real [Flu_NPatchGroup][NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#  endif
#  endif // FLU_SCHEME
#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   h_GramFE_TimeEvo = new gramfe_matmul_float [PS2][ 2*FLU_NXT ];
#  endif

} // FUNCTION : Init_MemAllocate_Fluid



#endif // #ifndef GPU
