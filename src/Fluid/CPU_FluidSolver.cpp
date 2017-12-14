#ifndef GPU



#include "GAMER.h"

#ifdef UNSPLIT_GRAVITY
#include "CUPOT.h"
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
#endif

#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
void CPU_FluidSolver_RTVD( real Flu_Array_In [][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                           real Flu_Array_Out[][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                           real Flux_Array[][9][NCOMP_TOTAL][ PS2*PS2 ],
                           const double Corner_Array[][3],
                           const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                           const int NPatchGroup, const real dt, const real dh, const real Gamma,
                           const bool StoreFlux, const bool XYZ, const real MinDens, const real MinPres );
#elif ( FLU_SCHEME == WAF )
void CPU_FluidSolver_WAF( real Flu_Array_In [][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out[][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          real Flux_Array[][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup, const real dt, const real dh, const real Gamma,
                          const bool StoreFlux, const bool XYZ, const WAF_Limiter_t WAF_Limiter,
                          const real MinDens, const real MinPres );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
void CPU_FluidSolver_MHM( const real Flu_Array_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out     [][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          char DE_Array_Out      [][ PS2*PS2*PS2 ],
                          real Flux_Array     [][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup, const real dt, const real dh, const real Gamma,
                          const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                          const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
                          const double ExtAcc_AuxArray[], const real MinDens, const real MinPres,
                          const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                          const bool JeansMinPres, const real JeansMinPres_Coeff );
#elif ( FLU_SCHEME == CTU )
void CPU_FluidSolver_CTU( const real Flu_Array_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out     [][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          char DE_Array_Out      [][ PS2*PS2*PS2 ],
                          real Flux_Array     [][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup, const real dt, const real dh, const real Gamma,
                          const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                          const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
                          const double ExtAcc_AuxArray[], const real MinDens, const real MinPres,
                          const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                          const bool JeansMinPres, const real JeansMinPres_Coeff );
#endif // FLU_SCHEME

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
void CPU_ELBDMSolver( real Flu_Array_In [][FLU_NIN    ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real Flu_Array_Out[][FLU_NOUT   ][ PS2*PS2*PS2 ],
                      real Flux_Array[][9][NFLUX_TOTAL][ PS2*PS2 ],
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                      const real Taylor3_Coeff, const bool XYZ, const real MinDens );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver
// Description :  1. MODEL == HYDRO : use CPU to solve the Euler equations by different schemes
//                2. MODEL == ELBDM : use CPU to solve the kinematic operator in the Schrodinger's equation
//
// Note        :  Currently five hydro schemes are supported in HYDRO :
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. Weighted-Average-Flux scheme                   (WAF   ) -->   split
//                   3. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   4. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   5. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
//
// Parameter   :  h_Flu_Array_In       : Host array storing the input fluid variables
//                h_Flu_Array_Out      : Host array to store the output fluid variables
//                h_DE_Array_Out       : Host array to store the dual-energy status
//                h_Flux_Array         : Host array to store the output fluxes (useful only if StoreFlux == true)
//                h_Corner_Array       : Host array storing the physical corner coordinates of each patch group
//                h_Pot_Array_USG      : Host array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup          : Number of patch groups to be evaluated
//                dt                   : Time interval to advance solution
//                dh                   : Grid size
//                Gamma                : Ratio of specific heats
//                StoreFlux            : true --> store the coarse-fine fluxes
//                XYZ                  : true   : x->y->z ( forward sweep)
//                                       false1 : z->y->x (backward sweep)
//                                       --> only useful for the RTVD and WAF schemes
//                LR_Limiter           : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                       (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                      vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff         : Coefficient of the generalized MinMod limiter
//                EP_Coeff             : Coefficient of the extrema-preserving limiter
//                WAF_Limiter          : Flux limiter for the WAF scheme
//                                       (0/1/2/3) = (SuperBee/vanLeer/vanAlbada/MinBee)
//                ELBDM_Eta            : Particle mass / Planck constant
//                ELBDM_Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion for ELBDM
//                ELBDM_Taylor3_Auto   : true --> Determine ELBDM_Taylor3_Coeff automatically by invoking the
//                                                function "ELBDM_SetTaylor3Coeff"
//                Time                 : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType          : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                MinDens/Pres         : Minimum allowed density and pressure
//                DualEnergySwitch     : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive          : true --> normalize passive scalars so that the sum of their mass density
//                                                is equal to the gas mass density
//                NNorm                : Number of passive scalars to be normalized
//                                       --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx              : Target variable indices to be normalized
//                                       --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres         : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff   : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//
// Useless parameters in HYDRO : ELBDM_Eta, ELBDM_Taylor3_Coeff, ELBDM_Taylor3_Auto
// Useless parameters in ELBDM : Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, WAF_Limiter, MinPres
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver( real h_Flu_Array_In [][FLU_NIN    ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                      real h_Flu_Array_Out[][FLU_NOUT   ][ PS2*PS2*PS2 ],
                      char h_DE_Array_Out[][ PS2*PS2*PS2 ],
                      real h_Flux_Array[][9][NFLUX_TOTAL][ PS2*PS2 ],
                      const double h_Corner_Array[][3],
                      const real h_Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                      const int NPatchGroup, const real dt, const real dh, const real Gamma, const bool StoreFlux,
                      const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                      const WAF_Limiter_t WAF_Limiter, const real ELBDM_Eta, real ELBDM_Taylor3_Coeff,
                      const bool ELBDM_Taylor3_Auto, const double Time, const OptGravityType_t GravityType,
                      const real MinDens, const real MinPres, const real DualEnergySwitch,
                      const bool NormPassive, const int NNorm, const int NormIdx[],
                      const bool JeansMinPres, const real JeansMinPres_Coeff )
{

// check
#  ifdef GAMER_DEBUG
   if ( StoreFlux  &&  h_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "Flux_Array is not allocated !!\n" );

#  ifdef UNSPLIT_GRAVITY
   if (  ( GravityType == GRAVITY_SELF || GravityType == GRAVITY_BOTH )  &&  h_Pot_Array_USG == NULL  )
   Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );

   if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  &&  h_Corner_Array == NULL  )
   Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
#  endif
#  endif

#  ifndef UNSPLIT_GRAVITY
   double ExtPot_AuxArray[1];
   double ExtAcc_AuxArray[1];
#  endif


#  if   ( MODEL == HYDRO )

#     if   ( FLU_SCHEME == RTVD )

      CPU_FluidSolver_RTVD( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, h_Corner_Array, h_Pot_Array_USG,
                            NPatchGroup, dt, dh, Gamma, StoreFlux, XYZ, MinDens, MinPres );

#     elif ( FLU_SCHEME == WAF )

      CPU_FluidSolver_WAF ( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, h_Corner_Array, h_Pot_Array_USG,
                            NPatchGroup, dt, dh, Gamma, StoreFlux, XYZ, WAF_Limiter, MinDens, MinPres );

#     elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

      CPU_FluidSolver_MHM ( h_Flu_Array_In, h_Flu_Array_Out, h_DE_Array_Out, h_Flux_Array, h_Corner_Array, h_Pot_Array_USG,
                            NPatchGroup, dt, dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff, Time,
                            GravityType, ExtAcc_AuxArray, MinDens, MinPres, DualEnergySwitch, NormPassive, NNorm, NormIdx,
                            JeansMinPres, JeansMinPres_Coeff );

#     elif ( FLU_SCHEME == CTU )

      CPU_FluidSolver_CTU ( h_Flu_Array_In, h_Flu_Array_Out, h_DE_Array_Out, h_Flux_Array, h_Corner_Array, h_Pot_Array_USG,
                            NPatchGroup, dt, dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff, Time,
                            GravityType, ExtAcc_AuxArray, MinDens, MinPres, DualEnergySwitch, NormPassive, NNorm, NormIdx,
                            JeansMinPres, JeansMinPres_Coeff );

#     else

#     error : unsupported CPU hydro scheme

#     endif


#  elif ( MODEL == MHD )
#     warning : WAIT MHD !!!


#  elif ( MODEL == ELBDM )
//    evaluate the optimized Taylor expansion coefficient
      if ( ELBDM_Taylor3_Auto )  ELBDM_Taylor3_Coeff = ELBDM_SetTaylor3Coeff( dt, dh, ELBDM_Eta );

      CPU_ELBDMSolver( h_Flu_Array_In, h_Flu_Array_Out, h_Flux_Array, NPatchGroup, dt, dh, ELBDM_Eta, StoreFlux,
                       ELBDM_Taylor3_Coeff, XYZ, MinDens );

#  else
#     error : ERROR : unsupported MODEL !!
#  endif // MODEL

} // FUNCTION : CPU_FluidSolver



#endif // #ifndef GPU
