#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU

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

#elif ( MODEL == ELBDM )

#if   ( WAVE_SCHEME == WAVE_FD )
__global__ void CUFLU_ELBDMSolver_FD( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                      const real Taylor3_Coeff, const bool XYZ, const real MinDens );
real ELBDM_SetTaylor3Coeff( const real dt, const real dh, const real Eta );
#elif ( WAVE_SCHEME == WAVE_GRAMFE )
#if   ( GRAMFE_SCHEME == GRAMFE_FFT )
__launch_bounds__(FFT::max_threads_per_block)
__global__
void CUFLU_ELBDMSolver_GramFE_FFT( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                   real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const bool XYZ, const real MinDens,
                                   typename FFT::workspace_type workspace,
                                   typename IFFT::workspace_type workspace_inverse );
#elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
void ELBDM_GramFE_ComputeTimeEvolutionMatrix( gramfe_matmul_float (*output)[ 2*FLU_NXT ], const real dt, const real dh, const real Eta );
__global__
void CUFLU_ELBDMSolver_GramFE_MATMUL( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      gramfe_matmul_float g_Evolve[][ FLU_NXT*2 ],
                                      const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                      const bool XYZ, const real MinDens );
#else
#  error : ERROR : unsupported GRAMFE_SCHEME !!
#endif // GRAMFE_SCHEME
#else
#  error : ERROR : unsupported WAVE_SCHEME !!
#endif // WAVE_SCHEME

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
__global__ void CUFLU_ELBDMSolver_HamiltonJacobi( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
#                                                 ifdef GAMER_DEBUG
                                                  real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
#                                                 else
                                                  real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
#                                                 endif
                                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                                  const bool h_IsCompletelyRefined[],
                                                  const bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                                  const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                                  const bool XYZ, const real MinDens );
#endif

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL

#if ( !defined GRAVITY  &&  MODEL == HYDRO )
static ExtAcc_t GPUExtAcc_Ptr = NULL;
#endif


// device pointers
extern real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ];
extern double (*d_Corner_Array_F)[3];
#if ( MODEL == HYDRO )
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ CUBE(PS2) ];
#else
static char (*d_DE_Array_F_Out)[ CUBE(PS2) ] = NULL;
#endif
#ifdef MHD
extern real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)         ];
extern real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ];
#else
static real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
static real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)          ] = NULL;
static real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ]            = NULL;
#endif
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*d_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*d_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*d_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#else
static real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
static real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ] = NULL;
#endif // MHD
#endif // FLU_SCHEME
#endif // #if ( MODEL == HYDRO )

#if ( MODEL == ELBDM )
extern bool (*d_IsCompletelyRefined);
#endif

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
extern bool (*d_HasWaveCounterpart)[ CUBE(HYB_NXT) ];
#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )

#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
extern gramfe_matmul_float (*d_Flu_TimeEvo)[2 * FLU_NXT];
#endif // #if ( GRAMFE_SCHEME == GRAMFE_MATMUL )

#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ];
#elif ( MODEL == HYDRO )
static real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ] = NULL;
#endif

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_FluidSolver
// Description :  1. MODEL == HYDRO : use GPU to solve the Euler equations by different schemes
//                                    --> invoke the kernel "CUFLU_FluidSolver_XXX"
//                2. MODEL == ELBDM : use GPU to solve the kinematic operator in the Schrodinger's equations
//                                    --> invoke the kernel "CUFLU_ELBDMSolver_XXX"
//
//                ***********************************************************
//                **                Asynchronous Function                  **
//                **                                                       **
//                **  will return before the execution in GPU is complete  **
//                ***********************************************************
//
// Note        :  1. Use streams for the asychronous memory copy between device and host
//                2. Prefix "d" : for pointers pointing to the "Device" memory space
//                   Prefix "h" : for pointers pointing to the "Host"   memory space
//                3. Use the input pamameter "XYZ" to control the order of update for dimensional-splitting
//                   method (currently only RTVD)
//                4. Currently five hydro schemes are supported :
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   3. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   4. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
// Parameter   :  h_Flu_Array_In        : Host array to store the input fluid variables
//                h_Flu_Array_Out       : Host array to store the output fluid variables
//                h_Mag_Array_In        : Host array storing the input B field (for MHD only)
//                h_Mag_Array_Out       : Host array to store the output B field (for MHD only)
//                h_DE_Array_Out        : Host array to store the dual-energy status
//                h_Flux_Array          : Host array to store the output fluxes
//                h_Ele_Array           : Host array to store the output electric field (for MHD only)
//                h_Corner_Array        : Host array storing the physical corner coordinates of each patch group
//                h_Pot_Array_USG       : Host array storing the input potential for UNSPLIT_GRAVITY
//                h_IsCompletelyRefined : Host array storing which patch groups are completely refined ( ELBDM only )
//                h_HasWaveCounterpart  : Host array storing which cells have wave counterpart ( ELBDM_HYBRID only )
//                NPatchGroup           : Number of patch groups evaluated simultaneously by GPU
//                dt                    : Time interval to advance solution
//                dh                    : Cell size
//                StoreFlux             : true --> store the coarse-fine fluxes
//                StoreElectric         : true --> store the coarse-fine electric field
//                XYZ                   : true  : x->y->z ( forward sweep)
//                                        false : z->y->x (backward sweep)
//                                        ~ useless in directionally unsplit schemes
//                LR_Limiter            : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                        (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                       vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff          : Coefficient of the generalized MinMod limiter
//                MinMod_MaxIter        : Maximum number of iterations to reduce MinMod_Coeff
//                ELBDM_Eta             : Particle mass / Planck constant
//                ELBDM_Taylor3_Coeff   : Coefficient in front of the third term in the Taylor expansion for ELBDM
//                ELBDM_Taylor3_Auto    : true --> Determine ELBDM_Taylor3_Coeff automatically by invoking the
//                                                 function "ELBDM_SetTaylor3Coeff"
//                Time                  : Current physical time                      (for UNSPLIT_GRAVITY only)
//                UsePot                : Add self-gravity and/or external potential (for UNSPLIT_GRAVITY only)
//                ExtAcc                : Add external acceleration                  (for UNSPLIT_GRAVITY only)
//                MicroPhy              : Microphysics object
//                MinDens/Pres/Eint     : Density, pressure, and internal energy floors
//                DualEnergySwitch      : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive           : true --> normalize passive scalars so that the sum of their mass density
//                                                 is equal to the gas mass density
//                NNorm                 : Number of passive scalars to be normalized
//                                        --> Should be set to the global variable "PassiveNorm_NVar"
//                FracPassive           : true --> convert passive scalars to mass fraction during data reconstruction
//                NFrac                 : Number of passive scalars for the option "FracPassive"
//                                        --> Should be set to the global variable "PassiveIntFrac_NVar"
//                JeansMinPres          : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff    : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                GPU_NStream           : Number of CUDA streams for the asynchronous memory copy
//                UseWaveFlag           : Determine whether to use wave or phase scheme
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In[][FLU_NIN ][ CUBE(FLU_NXT) ],
                             real h_Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                             real h_Mag_Array_In[][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
                             real h_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
                             char h_DE_Array_Out[][ CUBE(PS2) ],
                             real h_Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                             real h_Ele_Array[][9][NCOMP_ELE][ PS2P1*PS2 ],
                             const double h_Corner_Array[][3],
                             real h_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
                             const bool h_IsCompletelyRefined[],
                             const bool h_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                             const int NPatchGroup, const real dt, const real dh,
                             const bool StoreFlux, const bool StoreElectric,
                             const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter,
                             const real ELBDM_Eta, real ELBDM_Taylor3_Coeff, const bool ELBDM_Taylor3_Auto,
                             const double Time, const bool UsePot, const OptExtAcc_t ExtAcc, const MicroPhy_t MicroPhy,
                             const real MinDens, const real MinPres, const real MinEint,
                             const real DualEnergySwitch,
                             const bool NormPassive, const int NNorm,
                             const bool FracPassive, const int NFrac,
                             const bool JeansMinPres, const real JeansMinPres_Coeff,
                             const int GPU_NStream, const bool UseWaveFlag )
{

// check
#  ifdef GAMER_DEBUG
#  if   ( MODEL == HYDRO )

#  ifdef UNSPLIT_GRAVITY
   if ( UsePot )
   {
      if ( h_Pot_Array_USG   == NULL )    Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );
      if ( d_Pot_Array_USG_F == NULL )    Aux_Error( ERROR_INFO, "d_Pot_Array_USG_F == NULL !!\n" );
   }

   if ( ExtAcc )
   {
      if ( h_Corner_Array   == NULL )     Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
      if ( d_Corner_Array_F == NULL )     Aux_Error( ERROR_INFO, "d_Corner_Array_F == NULL !!\n" );
   }
#  endif

#  elif ( MODEL == ELBDM )
   if ( h_IsCompletelyRefined == NULL )   Aux_Error( ERROR_INFO, "h_IsCompletelyRefined == NULL !!\n" );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( h_HasWaveCounterpart == NULL  &&  !UseWaveFlag )
                                          Aux_Error( ERROR_INFO, "h_HasWaveCounterpart == NULL !!\n" );
#  endif

#  else
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif

   if ( StoreFlux )
   {
      if ( d_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "d_Flux_Array == NULL !!\n" );
      if ( h_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "h_Flux_Array == NULL !!\n" );
   }

#  ifdef MHD
   if ( h_Mag_Array_In    == NULL ) Aux_Error( ERROR_INFO, "h_Mag_Array_In == NULL !!\n" );
   if ( d_Mag_Array_F_In  == NULL ) Aux_Error( ERROR_INFO, "d_Mag_Array_F_In == NULL !!\n" );

   if ( h_Mag_Array_Out   == NULL ) Aux_Error( ERROR_INFO, "h_Mag_Array_Out == NULL !!\n" );
   if ( d_Mag_Array_F_Out == NULL ) Aux_Error( ERROR_INFO, "d_Mag_Array_F_Out == NULL !!\n" );

   if ( d_FC_Mag_Half     == NULL ) Aux_Error( ERROR_INFO, "d_FC_Mag_Half == NULL !!\n" );
   if ( d_EC_Ele          == NULL ) Aux_Error( ERROR_INFO, "d_EC_Ele == NULL !!\n" );

   if ( StoreElectric )
   {
      if ( d_Ele_Array == NULL )   Aux_Error( ERROR_INFO, "d_Ele_Array == NULL !!\n" );
      if ( h_Ele_Array == NULL )   Aux_Error( ERROR_INFO, "h_Ele_Array == NULL !!\n" );
   }
#  endif
#  endif // #ifdef GAMER_DEBUG


// thread block size
#  if (  !( MODEL == ELBDM  &&  WAVE_SCHEME == WAVE_GRAMFE  &&  GRAMFE_SCHEME == GRAMFE_FFT )  )
   const dim3 BlockDim_FluidSolver    ( FLU_BLOCK_SIZE_X, FLU_BLOCK_SIZE_Y,    1 ); // for the fluid solvers
#  endif
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   const dim3 BlockDim_FluidSolver_HJ ( FLU_BLOCK_SIZE_X, FLU_HJ_BLOCK_SIZE_Y, 1 ); // for the HJ solver
#  endif


// model-dependent operations
#  if   ( MODEL == HYDRO )

#  elif ( MODEL == ELBDM )

#  if ( WAVE_SCHEME == WAVE_GRAMFE  &&  GRAMFE_SCHEME == GRAMFE_FFT )
   uint cufftdx_shared_memory_size = NULL_INT;
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( UseWaveFlag ) {
#  endif

#  if   ( WAVE_SCHEME == WAVE_FD )

// evaluate the optimized Taylor expansion coefficient
   if ( ELBDM_Taylor3_Auto )  ELBDM_Taylor3_Coeff = ELBDM_SetTaylor3Coeff( dt, dh, ELBDM_Eta );

#  elif ( WAVE_SCHEME == WAVE_GRAMFE )

// set up GPU FFT if GPU is used for Gram Fourier extension FFT scheme
#  if   ( GRAMFE_SCHEME == GRAMFE_FFT )
// total size of shared memory required for storing FFT::ffts_per_block rows of data after Gram extension and the coefficients of the respective left and right extension polynomials
   auto size       = FFT::ffts_per_block*cufftdx::size_of<FFT>::value + 2*FFT::ffts_per_block*GRAMFE_NDELTA;
   auto size_bytes = size*sizeof(complex_type);

// shared memory must fit input data and must be big enough to run FFT
   cufftdx_shared_memory_size = std::max( (unsigned int)FFT::shared_memory_size, (unsigned int)size_bytes );

// increase max shared memory if needed
   CUDA_CHECK_ERROR(  cudaFuncSetAttribute( CUFLU_ELBDMSolver_GramFE_FFT, cudaFuncAttributeMaxDynamicSharedMemorySize,
                                            cufftdx_shared_memory_size )  );

#  elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
// time evolution matrix is copied to GPU in InvokeSolver()

#  else
#     error : ERROR : unsupported GRAMFE_SCHEME !!
#  endif // GRAMFE_SCHEME

#  else
#     error : ERROR : unsupported WAVE_SCHEME !!
#  endif // WAVE_SCHEME

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } // if ( UseWaveFlag )
#  endif

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   int *NPatch_per_Stream               = new int [GPU_NStream];
   int *UsedPatch                       = new int [GPU_NStream];
   int *Flu_MemSize_In                  = new int [GPU_NStream];
   int *Flu_MemSize_Out                 = new int [GPU_NStream];
   int *Flux_MemSize                    = new int [GPU_NStream];
#  ifdef MHD
   int *Mag_MemSize_In                  = new int [GPU_NStream];
   int *Mag_MemSize_Out                 = new int [GPU_NStream];
   int *Ele_MemSize                     = new int [GPU_NStream];
#  endif
#  ifdef UNSPLIT_GRAVITY
   int *USG_MemSize                     = new int [GPU_NStream];
   int *Corner_MemSize                  = new int [GPU_NStream];
#  endif
#  ifdef DUAL_ENERGY
   int *DE_MemSize_Out                  = new int [GPU_NStream];
#  endif
#  if ( MODEL == ELBDM )
   int *Flu_MemSize_IsCompletelyRefined = new int [GPU_NStream];
#  endif
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   int *Flu_MemSize_HasWaveCounterpart  = ( !UseWaveFlag ) ? new int [GPU_NStream] : NULL;
#  endif


// set the number of patches of each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatchGroup;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatchGroup / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatchGroup - UsedPatch[GPU_NStream-1];
   }


// set the size of data to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Flu_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN *CUBE(FLU_NXT);
      Flu_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NOUT*CUBE(PS2);
      Flux_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*NFLUX_TOTAL*9*SQR(PS2);
#     ifdef MHD
      Mag_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_MAG*FLU_NXT_P1*SQR(FLU_NXT);
      Mag_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_MAG*PS2P1*SQR(PS2);
      Ele_MemSize    [s] = sizeof(real  )*NPatch_per_Stream[s]*NCOMP_ELE*9*PS2P1*PS2;
#     endif
#     ifdef UNSPLIT_GRAVITY
      USG_MemSize    [s] = sizeof(real  )*NPatch_per_Stream[s]*CUBE(USG_NXT_F);
      Corner_MemSize [s] = sizeof(double)*NPatch_per_Stream[s]*3;
#     endif
#     ifdef DUAL_ENERGY
      DE_MemSize_Out [s] = sizeof(char  )*NPatch_per_Stream[s]*CUBE(PS2);
#     endif

//    optimization for phase scheme:
//    (a) transfer CUBE(HYB_NXT) instead of CUBE(FLU_NXT) cells to GPU
//    (b) when not in the debug mode, do not transfer STUB back from GPU (so only FLU_NIN instead of FLU_NOUT components)
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( !UseWaveFlag ) {
      Flu_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN*CUBE(HYB_NXT);
#     ifndef GAMER_DEBUG
      Flu_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN*CUBE(PS2);
#     endif
      }
#     endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )

#     if ( MODEL == ELBDM )
      Flu_MemSize_IsCompletelyRefined[s] = sizeof(bool)*NPatch_per_Stream[s];
#     endif
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( !UseWaveFlag )
      Flu_MemSize_HasWaveCounterpart [s] = sizeof(bool)*NPatch_per_Stream[s]*CUBE(HYB_NXT);
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( UseWaveFlag ) {
#     endif
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_F_In  + UsedPatch[s], h_Flu_Array_In  + UsedPatch[s],
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      } else {
      real (*smaller_d_Flu_Array_F_In)[FLU_NIN][CUBE(HYB_NXT)] = (real (*)[FLU_NIN][CUBE(HYB_NXT)]) d_Flu_Array_F_In;
      real (*smaller_h_Flu_Array_In  )[FLU_NIN][CUBE(HYB_NXT)] = (real (*)[FLU_NIN][CUBE(HYB_NXT)]) h_Flu_Array_In  ;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( smaller_d_Flu_Array_F_In + UsedPatch[s], smaller_h_Flu_Array_In + UsedPatch[s],
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
      }
#     endif
#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Mag_Array_F_In  + UsedPatch[s], h_Mag_Array_In  + UsedPatch[s],
                         Mag_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif

#     ifdef UNSPLIT_GRAVITY
      if ( UsePot )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_USG_F + UsedPatch[s], h_Pot_Array_USG + UsedPatch[s],
                         USG_MemSize   [s], cudaMemcpyHostToDevice, Stream[s] )  );

      if ( ExtAcc )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_F  + UsedPatch[s], h_Corner_Array  + UsedPatch[s],
                         Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif


#     if ( MODEL == ELBDM )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_IsCompletelyRefined + UsedPatch[s], h_IsCompletelyRefined + UsedPatch[s],
                         Flu_MemSize_IsCompletelyRefined[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif
#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( !UseWaveFlag )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_HasWaveCounterpart  + UsedPatch[s], h_HasWaveCounterpart  + UsedPatch[s],
                         Flu_MemSize_HasWaveCounterpart[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if   ( MODEL == HYDRO )

#        if   ( FLU_SCHEME == RTVD )

         CUFLU_FluidSolver_RTVD <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              dt, 1.0/dh, StoreFlux, XYZ, MinDens, MinPres, MinEint, EoS );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

         CUFLU_FluidSolver_MHM <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Mag_Array_F_In  + UsedPatch[s],
              d_Mag_Array_F_Out + UsedPatch[s],
              d_DE_Array_F_Out  + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Ele_Array       + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              d_PriVar          + UsedPatch[s],
              d_Slope_PPM       + UsedPatch[s],
              d_FC_Var          + UsedPatch[s],
              d_FC_Flux         + UsedPatch[s],
              d_FC_Mag_Half     + UsedPatch[s],
              d_EC_Ele          + UsedPatch[s],
              dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff, MinMod_MaxIter,
              Time, UsePot, ExtAcc, GPUExtAcc_Ptr, MinDens, MinPres, MinEint,
              DualEnergySwitch, NormPassive, NNorm, FracPassive, NFrac,
              JeansMinPres, JeansMinPres_Coeff, EoS, MicroPhy );

#        elif ( FLU_SCHEME == CTU )

         CUFLU_FluidSolver_CTU <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Mag_Array_F_In  + UsedPatch[s],
              d_Mag_Array_F_Out + UsedPatch[s],
              d_DE_Array_F_Out  + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Ele_Array       + UsedPatch[s],
              d_Corner_Array_F  + UsedPatch[s],
              d_Pot_Array_USG_F + UsedPatch[s],
              d_PriVar          + UsedPatch[s],
              d_Slope_PPM       + UsedPatch[s],
              d_FC_Var          + UsedPatch[s],
              d_FC_Flux         + UsedPatch[s],
              d_FC_Mag_Half     + UsedPatch[s],
              d_EC_Ele          + UsedPatch[s],
              dt, dh, StoreFlux, StoreElectric, LR_Limiter, MinMod_Coeff,
              Time, UsePot, ExtAcc, GPUExtAcc_Ptr, MinDens, MinPres, MinEint,
              DualEnergySwitch, NormPassive, NNorm, FracPassive, NFrac,
              JeansMinPres, JeansMinPres_Coeff, EoS );

#        else

#        error : unsupported GPU hydro scheme

#        endif // FLU_SCHEME

#     elif ( MODEL == ELBDM )

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      if ( UseWaveFlag ) {
#     endif

#     if   ( WAVE_SCHEME == WAVE_FD )

         CUFLU_ELBDMSolver_FD <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              dt, 1.0/dh, ELBDM_Eta, StoreFlux, ELBDM_Taylor3_Coeff, XYZ, MinDens );

#     elif ( WAVE_SCHEME == WAVE_GRAMFE )

#     if ( GRAMFE_SCHEME == GRAMFE_FFT )

//       create forward and backward cufftx workspaces
         cudaError_t error_code  = cudaSuccess;
         FFT::workspace_type cufftdx_workspace  = cufftdx::make_workspace<FFT>( error_code );
         CUDA_CHECK_ERROR(error_code);
         error_code              = cudaSuccess;
         IFFT::workspace_type cufftdx_iworkspace = cufftdx::make_workspace<IFFT>( error_code );
         CUDA_CHECK_ERROR(error_code);

         CUFLU_ELBDMSolver_GramFE_FFT <<< NPatch_per_Stream[s], FFT::block_dim, cufftdx_shared_memory_size, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              dt, 1.0/dh, ELBDM_Eta, StoreFlux, XYZ, MinDens, cufftdx_workspace, cufftdx_iworkspace );

#     elif ( GRAMFE_SCHEME == GRAMFE_MATMUL )
         CUFLU_ELBDMSolver_GramFE_MATMUL <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
            ( d_Flu_Array_F_In  + UsedPatch[s],
              d_Flu_Array_F_Out + UsedPatch[s],
              d_Flux_Array      + UsedPatch[s],
              d_Flu_TimeEvo,
              dt, dh, ELBDM_Eta, StoreFlux, XYZ, MinDens );
#     else
#        error : ERROR : unsupported GRAMFE_SCHEME !!
#     endif // GRAMFE_SCHEME

#     else // WAVE_SCHEME
#        error : ERROR : unsupported WAVE_SCHEME !!
#     endif // WAVE_SCHEME

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      } else {
         real (*smaller_d_Flu_Array_F_In) [FLU_NIN ][CUBE(HYB_NXT)] = (real (*)[FLU_NIN][CUBE(HYB_NXT)]) d_Flu_Array_F_In;
#        ifdef GAMER_DEBUG
         real (*smaller_d_Flu_Array_F_Out)[FLU_NOUT][CUBE(PS2)]     = d_Flu_Array_F_Out;
#        else
         real (*smaller_d_Flu_Array_F_Out)[FLU_NIN ][CUBE(PS2)]     = (real (*)[FLU_NIN][CUBE(PS2)]    ) d_Flu_Array_F_Out;
#        endif

         CUFLU_ELBDMSolver_HamiltonJacobi <<< NPatch_per_Stream[s], BlockDim_FluidSolver_HJ, 0, Stream[s] >>>
            (  smaller_d_Flu_Array_F_In  + UsedPatch[s],
               smaller_d_Flu_Array_F_Out + UsedPatch[s],
               d_Flux_Array              + UsedPatch[s],
               d_IsCompletelyRefined     + UsedPatch[s],
               d_HasWaveCounterpart      + UsedPatch[s],
               dt, 1.0/dh, ELBDM_Eta, StoreFlux, XYZ, MinDens );

      } // if ( UseWaveFlag ) ... else ...
#     endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )

#     else

#        error : unsupported MODEL !!

#     endif // MODEL


      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if ( ELBDM_SCHEME == ELBDM_HYBRID  &&  !defined(GAMER_DEBUG) )
      if ( UseWaveFlag ) {
#     endif
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array_Out + UsedPatch[s], d_Flu_Array_F_Out + UsedPatch[s],
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );
#     if ( ELBDM_SCHEME == ELBDM_HYBRID  &&  !defined(GAMER_DEBUG) )
      } else {
      real (*smaller_h_Flu_Array_Out  )[FLU_NIN][CUBE(PS2)] = (real (*)[FLU_NIN][CUBE(PS2)]) h_Flu_Array_Out;
      real (*smaller_d_Flu_Array_F_Out)[FLU_NIN][CUBE(PS2)] = (real (*)[FLU_NIN][CUBE(PS2)]) d_Flu_Array_F_Out;
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( smaller_h_Flu_Array_Out + UsedPatch[s], smaller_d_Flu_Array_F_Out + UsedPatch[s],
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );
      }
#     endif

      if ( StoreFlux )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flux_Array    + UsedPatch[s], d_Flux_Array      + UsedPatch[s],
                         Flux_MemSize[s],    cudaMemcpyDeviceToHost, Stream[s] )  );

#     ifdef MHD
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Mag_Array_Out + UsedPatch[s], d_Mag_Array_F_Out + UsedPatch[s],
                         Mag_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( StoreElectric )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Ele_Array     + UsedPatch[s], d_Ele_Array       + UsedPatch[s],
                         Ele_MemSize[s],    cudaMemcpyDeviceToHost, Stream[s] )  );
#     endif

#     ifdef DUAL_ENERGY
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_DE_Array_Out  + UsedPatch[s], d_DE_Array_F_Out  + UsedPatch[s],
                         DE_MemSize_Out[s],  cudaMemcpyDeviceToHost, Stream[s] )  );
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] Flu_MemSize_In;
   delete [] Flu_MemSize_Out;
   delete [] Flux_MemSize;
#  ifdef MHD
   delete [] Mag_MemSize_In;
   delete [] Mag_MemSize_Out;
   delete [] Ele_MemSize;
#  endif
#  ifdef UNSPLIT_GRAVITY
   delete [] USG_MemSize;
   delete [] Corner_MemSize;
#  endif
#  ifdef DUAL_ENERGY
   delete [] DE_MemSize_Out;
#  endif
#  if ( MODEL == ELBDM )
   delete [] Flu_MemSize_IsCompletelyRefined;
#  endif
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   delete [] Flu_MemSize_HasWaveCounterpart;
#  endif

} // FUNCTION : CUAPI_Asyn_FluidSolver



#endif // #ifdef GPU
