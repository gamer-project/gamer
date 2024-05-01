#include "CUAPI.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef GPU



// fluid solver prototypes in different models
#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
__global__ void CUFLU_FluidSolver_RTVD(
   real g_Fluid_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   real g_Fluid_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   real g_Flux     [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner[][3],
   const real g_Pot_USG[][ CUBE(USG_NXT_F) ],
   const real dt, const real _dh, const bool StoreFlux,
   const bool XYZ, const real MinDens, const real MinPres, const real MinEint,
   const EoS_t EoS );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
__global__
void CUFLU_FluidSolver_MHM(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool FracPassive, const int NFrac,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS, const MicroPhy_t MicroPhy );
#elif ( FLU_SCHEME == CTU )
__global__
void CUFLU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool FracPassive, const int NFrac,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS );
#endif // FLU_SCHEME
__global__ void CUFLU_dtSolver_HydroCFL( real g_dt_Array[], const real g_Flu_Array[][FLU_NIN_T][ CUBE(PS1) ],
                                         const real g_Mag_Array[][NCOMP_MAG][ PS1P1*SQR(PS1) ],
                                         const real dh, const real Safety, const real MinPres, const EoS_t EoS,
                                         const MicroPhy_t MicroPhy );
#ifdef GRAVITY
__global__
void CUPOT_dtSolver_HydroGravity( real g_dt_Array[], const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                                  const double g_Corner_Array[][3],
                                  const real dh, const real Safety, const bool P5_Gradient,
                                  const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                                  const double ExtAcc_Time );
#endif

#elif ( MODEL == ELBDM )
# if   ( WAVE_SCHEME == WAVE_FD )
__global__ void CUFLU_ELBDMSolver_FD( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                      const real Taylor3_Coeff, const bool XYZ, const real MinDens );
# elif ( WAVE_SCHEME == WAVE_GRAMFE )
#  if   ( GRAMFE_SCHEME == GRAMFE_FFT )
__global__ void CUFLU_ELBDMSolver_GramFE_FFT( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                              real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                              real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                              const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                              const bool XYZ, const real MinDens,
                                              typename FFT::workspace_type  Workspace,
                                              typename IFFT::workspace_type WorkspaceInv );
#  elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
__global__ void CUFLU_ELBDMSolver_GramFE_MATMUL( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                                 real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                                 real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                                 gramfe_matmul_float g_TimeEvo[][ FLU_NXT*2 ],
                                                 const real dt, const real dh, const real Eta, const bool StoreFlux,
                                                 const bool XYZ, const real MinDens );
#  else // GRAMFE_SCHEME
#     error : ERROR : unsupported GRAMFE_SCHEME !!
#  endif // GRAMFE_SCHEME
# else // WAVE_SCHEME
#  error : ERROR : unsupported WAVE_SCHEME !!
# endif // WAVE_SCHEME

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
__global__ void CUFLU_ELBDMSolver_HamiltonJacobi( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
                                                  #ifdef GAMER_DEBUG
                                                  real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                                  #else
                                                  real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                                                  #endif
                                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                                  const bool g_IsCompletelyRefined[],
                                                  const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                                  const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                                  const bool XYZ, const real MinDens );
#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )

#else // MODEL
#error : ERROR : unsupported MODEL !!
#endif // MODEL


#ifdef GRAVITY

// Poisson solver prototypes
#if   ( POT_SCHEME == SOR )
__global__ void CUPOT_PoissonSolver_SOR( const real g_Rho_Array    [][ CUBE(RHO_NXT) ],
                                         const real g_Pot_Array_In [][ CUBE(POT_NXT) ],
                                               real g_Pot_Array_Out[][ CUBE(GRA_NXT) ],
                                         const int Min_Iter, const int Max_Iter, const real Omega_6,
                                         const real Const, const IntScheme_t IntScheme );
#elif ( POT_SCHEME == MG )
__global__ void CUPOT_PoissonSolver_MG( const real g_Rho_Array    [][ CUBE(RHO_NXT) ],
                                        const real g_Pot_Array_In [][ CUBE(POT_NXT) ],
                                              real g_Pot_Array_Out[][ CUBE(GRA_NXT) ],
                                        const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                                        const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                                        const IntScheme_t IntScheme );
#endif // POT_SCHEME


// Gravity solver prototypes in different models
#if   ( MODEL == HYDRO )
__global__
void CUPOT_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN-1][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real   g_Emag_Array   [][ CUBE(PS1) ],
   const real dt, const real dh, const bool P5_Gradient,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double TimeNew, const double TimeOld, const real MinEint );

#elif ( MODEL == ELBDM )
__global__
void CUPOT_ELBDMGravitySolver(       real g_Flu_Array[][GRA_NIN][ CUBE(PS1) ],
                               const real g_Pot_Array[][ CUBE(GRA_NXT) ],
                               const real EtaDt, const real dh, const real Lambda );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

#endif // GRAVITY


// source-term solver prototype
__global__
void CUSRC_SrcSolver_IterateAllCells(
   const real g_Flu_Array_In [][FLU_NIN_S ][ CUBE(SRC_NXT)           ],
         real g_Flu_Array_Out[][FLU_NOUT_S][ CUBE(PS1)               ],
   const real g_Mag_Array_In [][NCOMP_MAG ][ SRC_NXT_P1*SQR(SRC_NXT) ],
   const double g_Corner_Array[][3],
   const SrcTerms_t SrcTerms, const int NPatchGroup, const real dt, const real dh,
   const double TimeNew, const double TimeOld,
   const real MinDens, const real MinPres, const real MinEint, const EoS_t EoS );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetCache
// Description :  Set cache preference
//
// Parameter   :
//-------------------------------------------------------------------------------------------------------
void CUAPI_SetCache()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. fluid solver
#  if   ( MODEL == HYDRO )
#  if   ( FLU_SCHEME == RTVD )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_RTVD,             cudaFuncCachePreferShared )  );
#  elif ( FLU_SCHEME == MHM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,              cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == MHM_RP )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,              cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == CTU )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_CTU,              cudaFuncCachePreferL1     )  );
#  endif
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_dtSolver_HydroCFL,            cudaFuncCachePreferShared )  );
#  ifdef GRAVITY
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_dtSolver_HydroGravity,        cudaFuncCachePreferShared )  );
#  endif

#  elif ( MODEL == ELBDM )
#  if   ( WAVE_SCHEME == WAVE_FD )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver_FD,               cudaFuncCachePreferShared )  );
#  elif ( WAVE_SCHEME == WAVE_GRAMFE )
#   if   ( GRAMFE_SCHEME == GRAMFE_FFT )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver_GramFE_FFT,       cudaFuncCachePreferShared )  );
#   elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver_GramFE_MATMUL,    cudaFuncCachePreferShared )  );
#   else // GRAMFE_SCHEME
#   error : ERROR : unsupported GRAMFE_SCHEME !!
#   endif // GRAMFE_SCHEME
#  else // WAVE_SCHEME
#  error : ERROR : unsupported WAVE_SCHEME !!
#  endif // WAVE_SCHEME
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver_HamiltonJacobi,   cudaFuncCachePreferShared )  );
#  endif

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY

// 2. Poisson solver
#  if   ( POT_SCHEME == SOR )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_PoissonSolver_SOR,            cudaFuncCachePreferShared )  );
#  elif ( POT_SCHEME == MG )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_PoissonSolver_MG,             cudaFuncCachePreferShared )  );
#  endif // POT_SCHEME


// 3. gravity solver
#  if   ( MODEL == HYDRO )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_HydroGravitySolver,           cudaFuncCachePreferShared )  );

#  elif ( MODEL == ELBDM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUPOT_ELBDMGravitySolver,           cudaFuncCachePreferL1     )  );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  endif // GRAVITY


// 4. source-term solver
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUSRC_SrcSolver_IterateAllCells,   cudaFuncCachePreferL1      )  );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_SetCache



#endif // #ifdef GPU
