#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef GPU



#ifdef GRAVITY
extern double ExtPot_AuxArray[EXT_POT_NAUX_MAX];
extern double ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];
#endif


// fluid solver prototypes in different models
#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
__global__ void CUFLU_FluidSolver_RTVD( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                        real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                        real g_Flux[][9][5][ PS2*PS2 ], 
                                        const double g_Corner[][3],
                                        const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                        const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                        const bool XYZ );
#elif ( FLU_SCHEME == WAF )
__global__ void CUFLU_FluidSolver_WAF( real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                       real g_Flux[][9][5][ PS2*PS2 ], 
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const bool XYZ, const WAF_Limiter_t WAF_Limiter );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
__global__ void CUFLU_FluidSolver_MHM( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[] [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff, const double Time, const OptGravityType_t GravityType );
#ifdef UNSPLIT_GRAVITY
int CUFLU_FluidSolver_SetConstMem( double ExtAcc_AuxArray_h[] );
#endif
#elif ( FLU_SCHEME == CTU )
__global__ void CUFLU_FluidSolver_CTU( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[] [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff, const double Time, const OptGravityType_t GravityType );
#ifdef UNSPLIT_GRAVITY
int CUFLU_FluidSolver_SetConstMem( double ExtAcc_AuxArray_h[] );
#endif
#endif
//__global__ void CUFLU_GetMaxCFL( real g_Fluid[][5][ PS2*PS2*PS2 ], real g_MaxCFL[], const real Gamma );

#elif ( MODEL == MHD )
#warning :: WAIT MHD !!!

#elif ( MODEL == ELBDM )
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   real g_Flux     [][9][NFLUX][ PS2*PS2 ], 
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux, 
                                   const real Taylor3_Coeff, const bool XYZ );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL


#ifdef GRAVITY

// Poisson solver prototypes
#if   ( POT_SCHEME == SOR )
#ifdef USE_PSOLVER_10TO14
__global__ void CUPOT_PoissonSolver_SOR_10to14cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ], 
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ], 
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6,
                                                    const real Const, const IntScheme_t IntScheme );
#else
__global__ void CUPOT_PoissonSolver_SOR_16to18cube( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ], 
                                                    const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ], 
                                                          real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                                    const int Min_Iter, const int Max_Iter, const real Omega_6, 
                                                    const real Const, const IntScheme_t IntScheme );
#endif // #ifdef USE_PSOLVER_10TO14 ... else ...
#elif ( POT_SCHEME == MG )
__global__ void CUPOT_PoissonSolver_MG( const real g_Rho_Array    [][ RHO_NXT*RHO_NXT*RHO_NXT ], 
                                        const real g_Pot_Array_In [][ POT_NXT*POT_NXT*POT_NXT ], 
                                              real g_Pot_Array_Out[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                        const real dh_Min, const int Max_Iter, const int NPre_Smooth,
                                        const int NPost_Smooth, const real Tolerated_Error, const real Poi_Coeff,
                                        const IntScheme_t IntScheme );
#endif // POT_SCHEME


// Gravity solver prototypes in different models
#if   ( MODEL == HYDRO )
__global__ void CUPOT_HydroGravitySolver(       real g_Flu_Array_New[][NCOMP][ PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ],
                                          const real g_Pot_Array_New[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                          const double g_Corner_Array[][3],
                                          const real g_Pot_Array_USG[][ USG_NXT_G*USG_NXT_G*USG_NXT_G ],
                                          const real g_Flu_Array_USG[][GRA_NIN-1][ PS1*PS1*PS1 ],
                                          const real Gra_Const, const bool P5_Gradient, const OptGravityType_t GravityType, 
                                          const double TimeNew, const double TimeOld, const real dt, const real dh );
int CUPOT_HydroGravitySolver_SetConstMem( double ExtAcc_AuxArray_h[] );

#elif ( MODEL == MHD )
#warning :: WAIT MHD !!!

#elif ( MODEL == ELBDM )
__global__ void CUPOT_ELBDMGravitySolver(       real g_Flu_Array[][GRA_NIN][ PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ],
                                          const real g_Pot_Array[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                          const double g_Corner_Array[][3],
                                          const real EtaDt, const real dh, const real Lambda, const bool ExtPot,
                                          const double Time );
int CUPOT_ELBDMGravitySolver_SetConstMem( double ExtPot_AuxArray_h[] );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

int CUPOT_PoissonSolver_SetConstMem();

#endif // GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Set_Default_GPU_Parameter
// Description :  Set several GPU parameters to the default values if they are not set in the input file
//
// Parameter   :  GPU_NStream       : Number of streams for the asynchronous memory copy in GPU 
//                Flu_GPU_NPGroup   : Number of patch groups sent into GPU simultaneously for the fluid solver
//                Pot_GPU_NPGroup   : Number of patch groups sent into GPU simultaneously for the Poisson solver
//-------------------------------------------------------------------------------------------------------
void CUAPI_Set_Default_GPU_Parameter( int &GPU_NStream, int &Flu_GPU_NPGroup, int &Pot_GPU_NPGroup )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... \n", __FUNCTION__ );


// get the device ID
   int GetDeviceID = 999;
   CUDA_CHECK_ERROR(  cudaGetDevice( &GetDeviceID )  );


// load the device properties
   cudaDeviceProp DeviceProp;
   CUDA_CHECK_ERROR(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );


// set the default GPU parameters
// (1) GPU_NSTREAM  
   if ( GPU_NStream < 0 )
   {
      if ( DeviceProp.deviceOverlap )     
      {
#        if   ( MODEL == HYDRO )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 8;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 48;             // optimized for K40
#           else
#           error : UNKNOWN GPU_ARCH !!
#           endif

#        elif ( MODEL == MHD )
#        warning :: WAIT MHD !!!

#        elif ( MODEL == ELBDM )
#           if   ( GPU_ARCH == FERMI )
            GPU_NStream = 8;
#           elif ( GPU_ARCH == KEPLER )
            GPU_NStream = 48;             // optimized for K40
#           else
#           error : ERROR : UNKNOWN GPU_ARCH !!
#           endif
#        else
#           error : ERROR : UNKNOWN MODEL !!
#        endif // MODEL
      }
      else                                
         GPU_NStream = 1;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "GPU_NSTREAM", GPU_NSTREAM );
   } // if ( GPU_NStream < 0 )

// (2) FLU_GPU_NPGROUP
   if ( Flu_GPU_NPGroup < 0 )  
   {
#     if   ( MODEL == HYDRO )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 32*DeviceProp.multiProcessorCount;  // optimized for K40 (for GPU_NStream=48, NPGroup=480)
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif

#     elif ( MODEL == MHD )
#        warning :: WAIT MHD !!!

#     elif ( MODEL == ELBDM )
#        if   ( GPU_ARCH == FERMI )
         Flu_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#        elif ( GPU_ARCH == KEPLER )
         Flu_GPU_NPGroup = 32*DeviceProp.multiProcessorCount;  // optimized for K40 (for GPU_NStream=48, NPGroup=480)
#        else
#        error : UNKNOWN GPU_ARCH !!
#        endif
#     else
#        error : ERROR : UNKNOWN MODEL !!
#     endif // MODEL

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "FLU_GPU_NPGROUP", Flu_GPU_NPGroup );
   } // if ( Flu_GPU_NPGroup < 0 )

// (3) POT_GPU_NPGROUP
#  ifdef GRAVITY
   if ( Pot_GPU_NPGroup < 0 )  
   {
#     if   ( GPU_ARCH == FERMI )
      Pot_GPU_NPGroup = 1*GPU_NStream*DeviceProp.multiProcessorCount;
#     elif ( GPU_ARCH == KEPLER )
      Pot_GPU_NPGroup = 32*DeviceProp.multiProcessorCount;  // optimized for K40 (for GPU_NStream=48, NPGroup=480)
#     else
#     error : UNKNOWN GPU_ARCH !!
#     endif

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "POT_GPU_NPGROUP", Pot_GPU_NPGroup );
   } // if ( Pot_GPU_NPGroup < 0 )
#  endif

// (4) cache preference
// (4-1) fluid solver
#  if   ( MODEL == HYDRO )
#  if   ( FLU_SCHEME == RTVD )   
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_RTVD, cudaFuncCachePreferShared )  );
#  elif ( FLU_SCHEME == WAF )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_WAF,  cudaFuncCachePreferShared )  );
#  elif ( FLU_SCHEME == MHM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,  cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == MHM_RP )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_MHM,  cudaFuncCachePreferL1     )  );
#  elif ( FLU_SCHEME == CTU )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_FluidSolver_CTU,  cudaFuncCachePreferL1     )  );
#  endif

// CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_GetMaxCFL,        cudaFuncCachePreferShared )  );

#  elif ( MODEL == MHD )
#  warning :: WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   CUDA_CHECK_ERROR(  cudaFuncSetCacheConfig( CUFLU_ELBDMSolver,      cudaFuncCachePreferShared )  );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


#  ifdef GRAVITY      

// (4-2) Poisson solver
#  if   ( POT_SCHEME == SOR )
#  ifdef USE_PSOLVER_10TO14
   CUDA_CHECK_ERROR( cudaFuncSetCacheConfig( CUPOT_PoissonSolver_SOR_10to14cube, cudaFuncCachePreferShared ) );
#  else
   CUDA_CHECK_ERROR( cudaFuncSetCacheConfig( CUPOT_PoissonSolver_SOR_16to18cube, cudaFuncCachePreferShared ) );
#  endif
#  elif ( POT_SCHEME == MG )
   CUDA_CHECK_ERROR( cudaFuncSetCacheConfig( CUPOT_PoissonSolver_MG,             cudaFuncCachePreferShared ) );
#  endif // POT_SCHEME


// (4-3) gravity solver
#  if   ( MODEL == HYDRO )
   CUDA_CHECK_ERROR( cudaFuncSetCacheConfig( CUPOT_HydroGravitySolver,           cudaFuncCachePreferShared ) );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   CUDA_CHECK_ERROR( cudaFuncSetCacheConfig( CUPOT_ELBDMGravitySolver,           cudaFuncCachePreferL1     ) );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

#  endif // GRAVITY


// (5) set the constant variables
#  ifdef GRAVITY      
   if ( CUPOT_PoissonSolver_SetConstMem() != 0 )
   Aux_Error( ERROR_INFO, "CUPOT_PoissonSolver_SetConstMem failed ...\n" );

#  if ( MODEL == HYDRO )
   if (  ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&
         CUPOT_HydroGravitySolver_SetConstMem( ExtAcc_AuxArray ) != 0  )
   Aux_Error( ERROR_INFO, "CUPOT_HydroGravitySolver_SetConstMem failed ...\n" );

#  if (  defined UNSPLIT_GRAVITY  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )
   if (  ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL || OPT__GRAVITY_TYPE == GRAVITY_BOTH )  &&
         CUFLU_FluidSolver_SetConstMem( ExtAcc_AuxArray ) != 0  )
   Aux_Error( ERROR_INFO, "CUFLU_FluidSolver_SetConstMem failed ...\n" );
#  endif
#  endif // if ( MODEL == HYDRO )

#  if ( MODEL == ELBDM )
   if (  OPT__EXTERNAL_POT  &&  CUPOT_ELBDMGravitySolver_SetConstMem( ExtPot_AuxArray ) != 0  )
   Aux_Error( ERROR_INFO, "CUPOT_ELBDMGravitySolver_SetConstMem failed ...\n" );
#  endif
#  endif
   

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_Set_Default_GPU_Parameter



#endif // #ifdef GPU
