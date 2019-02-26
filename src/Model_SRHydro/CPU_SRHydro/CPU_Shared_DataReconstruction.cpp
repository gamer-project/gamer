#ifndef __CUFLU_DATARECONSTRUCTION__
#define __CUFLU_DATARECONSTRUCTION__



#include "CUFLU.h"

#if (  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else
# include "../../../include/SRHydroPrototypes.h"
#endif // #ifdef __CUDACC__ ... else ...


// internal functions (GPU_DEVICE is defined in CUFLU.h)
GPU_DEVICE
static void SRHydro_LimitSlope( const real L1[], const real C0[], const real R1[], const LR_Limiter_t LR_Limiter,
                                const real MinMod_Coeff, const real Gamma, const int XYZ, real Slope_Limiter[] );
#if ( FLU_SCHEME == MHM )
GPU_DEVICE
static void SRHydro_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                                    const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                                    const real MinDens, const real MinTemp, const real Gamma );
#endif



#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. Input data can be either conserved or primitive variables
//                   --> If the input data are conserved variables, one must provide g_ConVar[] and enable "Con2Pri"
//                       --> Primitive variables will be calculated by this function and stored in g_PriVar[]
//                       --> g_PriVar[] must be allocated in advance but it doesn't need to be initialized
//                       --> Adopted by MHM
//                   --> If the input data are primitive variables, one must provide g_PriVar[] and disable "Con2Pri"
//                       --> g_ConVar[] is useless here
//                           --> Although SRHydro_HancockPredict() still needs g_ConVar(), MHM provides conserved instead
//                               of primitive variables. So it's fine.
//                       --> Adopted by MHM_RP, where SRHydro_RiemannPredict() already returns primitive variables
//                3. Output data are always conserved variables
//                   --> Because SRHydro_HancockPredict() only works with conserved variables
//                4. PLM and PPM data reconstruction functions share the same function name
//                5. Face-centered variables will be advanced by half time-step for the MHM schemes
//                6. Data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION" in the header CUFLU.h
//                7. This function is shared by MHM, MHM_RP schemes
//
//
// Parameter   :  g_ConVar           : Array storing the input cell-centered conserved variables
//                g_PriVar           : Array storing/to store the cell-centered primitive variables
//                                     --> For MHM, g_ConVar[] and g_PriVar[] must point to different arrays since
//                                         SRHydro_HancockPredict() requires the original g_ConVar[]
//                g_FC_Var           : Array to store the output face-centered primitive variables
//                g_Slope_PPM        : Array to store the x/y/z slopes for the PPM reconstruction
//                                     --> Useless for PLM
//                NIn                : Size of g_PriVar[] along each direction
//                                     --> Can be smaller than FLU_NXT
//                NGhost             : Number of ghost zones
//                                      --> "NIn-2*NGhost" cells will be computed along each direction
//                                      --> Size of g_FC_Var[] is assumed to be "(NIn-2*NGhost)^3"
//                                      --> The reconstructed data at cell (i,j,k) will be stored in g_FC_Var[]
//                                          with the index "(i-NGhost,j-NGhost,k-NGhost)"
//                Gamma              : Ratio of specific heats
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP schemes
//                                     (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                dt                 : Time interval to advance solution (for the CTU scheme)
//                dh                 : Cell size
//                MinDens/Temp       : Minimum allowed density and temperature
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                       real g_PriVar   [][ CUBE(FLU_NXT) ],
                                       real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                       real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                                 const int NIn, const int NGhost, const real Gamma,
                                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                 const real dt, const real dh, const real MinDens, const real MinTemp )
{

   const int  didx_cc[3] = { 1, NIn, SQR(NIn) };
   const int  NOut       = NIn - 2*NGhost;      // number of output cells

// 0. conserved --> primitive variables
#     if ( FLU_SCHEME == MHM )
      real ConVar_1Cell[NCOMP_TOTAL], PriVar_1Cell[NCOMP_TOTAL];

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

#        ifdef CHECK_NEGATIVE_IN_FLUID
         SRHydro_CheckUnphysical( ConVar_1Cell, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true );
#        endif

         SRHydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, Gamma, MinTemp );

         for (int v=0; v<NCOMP_TOTAL; v++)   g_PriVar[v][idx] = PriVar_1Cell[v];
      }

#     ifdef __CUDACC__
      __syncthreads();
#     endif
#     endif

// data reconstruction
   const int NOut2 = SQR(NOut);
   CGPU_LOOP( idx_fc, CUBE(NOut) )
   {
      const int i_cc   = NGhost + idx_fc%NOut;
      const int j_cc   = NGhost + idx_fc%NOut2/NOut;
      const int k_cc   = NGhost + idx_fc/NOut2;
      const int idx_cc = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

      real cc_C[NCOMP_TOTAL], cc_L[NCOMP_TOTAL], cc_R[NCOMP_TOTAL];  // cell-centered variables of the Central/Left/Right cells
      real fc[6][NCOMP_TOTAL];                                       // face-centered variables of the central cell
      real Slope_Limiter[NCOMP_TOTAL];

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C[v] = g_PriVar[v][idx_cc];


//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
//       2. evaluate the monotonic slope
         const int faceL   = 2*d;      // left and right face indices
         const int faceR   = faceL+1;
         const int idx_ccL = idx_cc - didx_cc[d];
         const int idx_ccR = idx_cc + didx_cc[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            cc_L[v] = g_PriVar[v][idx_ccL];
            cc_R[v] = g_PriVar[v][idx_ccR];
         }

         SRHydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, Gamma, d, Slope_Limiter );


//       3. get the face-centered primitive variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fc[faceL][v] = cc_C[v] - (real)0.5*Slope_Limiter[v];
            fc[faceR][v] = cc_C[v] + (real)0.5*Slope_Limiter[v];
         }

//       ensure the face-centered variables lie between neighboring cell-centered values
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            real Min, Max;

            Min = ( cc_C[v] < cc_L[v] ) ? cc_C[v] : cc_L[v];
            Max = ( cc_C[v] > cc_L[v] ) ? cc_C[v] : cc_L[v];
            fc[faceL][v] = ( fc[faceL][v] > Min ) ? fc[faceL][v] : Min;
            fc[faceL][v] = ( fc[faceL][v] < Max ) ? fc[faceL][v] : Max;
            fc[faceR][v] = (real)2.0*cc_C[v] - fc[faceL][v];

            Min = ( cc_C[v] < cc_R[v] ) ? cc_C[v] : cc_R[v];
            Max = ( cc_C[v] > cc_R[v] ) ? cc_C[v] : cc_R[v];
            fc[faceR][v] = ( fc[faceR][v] > Min ) ? fc[faceR][v] : Min;
            fc[faceR][v] = ( fc[faceR][v] < Max ) ? fc[faceR][v] : Max;
            fc[faceL][v] = (real)2.0*cc_C[v] - fc[faceR][v];
         }

//       5. primitive variables --> conserved variables
         real tmp[NCOMP_TOTAL];  // input and output arrays must not overlap for Pri2Con()

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceL][v];
         SRHydro_Pri2Con( tmp, fc[faceL], Gamma );

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceR][v];
         SRHydro_Pri2Con( tmp, fc[faceR], Gamma );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    6. advance the face-centered variables by half time-step for the MHM integrator
      SRHydro_HancockPredict( fc, dt, dh, g_ConVar, idx_cc, MinDens, MinTemp, Gamma );
#     endif


//    7. store the face-centered values to the output array
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(NOut) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : SRHydro_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  See the PLM routine
//
// Parameter   :  See the PLM routine
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                       real g_PriVar   [][ CUBE(FLU_NXT) ],
                                       real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                       real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                                 const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                 const real dt, const real dh, const real MinDens, const real MinTemp )
{

   const int  didx_cc   [3] = { 1, NIn, SQR(NIn) };
   const int  didx_slope[3] = { 1, N_SLOPE_PPM, SQR(N_SLOPE_PPM) };
   const int  NOut       = NIn - 2*NGhost;      // number of output cells
   const int  NSlope     = N_SLOPE_PPM;         // size of g_Slope_PPM[] (which must be equal to NOut + 2)


// 0. conserved --> primitive variables
   if ( Con2Pri )
   {
      real ConVar_1Cell[NCOMP_TOTAL], PriVar_1Cell[NCOMP_TOTAL];

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

         SRHydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, Gamma, MinTemp );

         for (int v=0; v<NCOMP_TOTAL; v++)   g_PriVar[v][idx] = PriVar_1Cell[v];
      }

#     ifdef __CUDACC__
      __syncthreads();
#     endif
   } // if ( Con2Pri )


// 1. evaluate the monotonic slope of all cells
   const int NSlope2 = SQR(NSlope);
   CGPU_LOOP( idx_slope, CUBE(NSlope) )
   {
      const int i_cc   = NGhost - 1 + idx_slope%NSlope;
      const int j_cc   = NGhost - 1 + idx_slope%NSlope2/NSlope;
      const int k_cc   = NGhost - 1 + idx_slope/NSlope2;
      const int idx_cc = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

      real cc_C[NCOMP_TOTAL], cc_L[NCOMP_TOTAL], cc_R[NCOMP_TOTAL];  // cell-centered variables of the Central/Left/Right cells
      real Slope_Limiter[NCOMP_TOTAL];

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C[v] = g_PriVar[v][idx_cc];

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         const int idx_ccL = idx_cc - didx_cc[d];
         const int idx_ccR = idx_cc + didx_cc[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            cc_L[v] = g_PriVar[v][idx_ccL];
            cc_R[v] = g_PriVar[v][idx_ccR];
         }

         SRHydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, Gamma, d, Slope_Limiter );

//       store the results to g_Slope_PPM[]
         for (int v=0; v<NCOMP_TOTAL; v++)   g_Slope_PPM[d][v][idx_slope] = Slope_Limiter[v];

      } // for (int d=0; d<3; d++)
   } // CGPU_LOOP( idx_slope, CUBE(NSlope) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif


// data reconstruction
   const int NOut2 = SQR(NOut);
   CGPU_LOOP( idx_fc, CUBE(NOut) )
   {
      const int i_fc      = idx_fc%NOut;
      const int j_fc      = idx_fc%NOut2/NOut;
      const int k_fc      = idx_fc/NOut2;

      const int i_cc      = i_fc + NGhost;
      const int j_cc      = j_fc + NGhost;
      const int k_cc      = k_fc + NGhost;
      const int idx_cc    = IDX321( i_cc, j_cc, k_cc, NIn, NIn );

      const int i_slope   = i_fc + 1;   // because NSlope = NOut + 2
      const int j_slope   = j_fc + 1;
      const int k_slope   = k_fc + 1;
      const int idx_slope = IDX321( i_slope, j_slope, k_slope, NSlope, NSlope );

 //   cc/fc: cell/face-centered variables; _C_ncomp: central cell with all NCOMP_TOTAL variables
      real cc_C_ncomp[NCOMP_TOTAL], fc[6][NCOMP_TOTAL], dfc[NCOMP_TOTAL], dfc6[NCOMP_TOTAL];

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C_ncomp[v] = g_PriVar[v][idx_cc];


//    3. get the face-centered primitive variables
//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         const int faceL      = 2*d;      // left and right face indices
         const int faceR      = faceL+1;
         const int idx_ccL    = idx_cc - didx_cc[d];
         const int idx_ccR    = idx_cc + didx_cc[d];
         const int idx_slopeL = idx_slope - didx_slope[d];
         const int idx_slopeR = idx_slope + didx_slope[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          cc/fc: cell/face-centered variables; _C/L/R: Central/Left/Right cells
            real cc_C, cc_L, cc_R, dcc_L, dcc_R, dcc_C, fc_L, fc_R, Max, Min;

//          3-1. parabolic interpolation
            cc_L  = g_PriVar[v][idx_ccL];
            cc_R  = g_PriVar[v][idx_ccR];
            cc_C  = cc_C_ncomp[v];

            dcc_L = g_Slope_PPM[d][v][idx_slopeL];
            dcc_R = g_Slope_PPM[d][v][idx_slopeR];
            dcc_C = g_Slope_PPM[d][v][idx_slope ];

            fc_L  = (real)0.5*( cc_C + cc_L ) - (real)1.0/(real)6.0*( dcc_C - dcc_L );
            fc_R  = (real)0.5*( cc_C + cc_R ) - (real)1.0/(real)6.0*( dcc_R - dcc_C );


//          3-2. monotonicity constraint
            dfc [v] = fc_R - fc_L;
            dfc6[v] = (real)6.0*(  cc_C - (real)0.5*( fc_L + fc_R )  );

            if (  ( fc_R - cc_C )*( cc_C - fc_L ) <= (real)0.0  )
            {
               fc_L = cc_C;
               fc_R = cc_C;
            }
            else if ( dfc[v]*dfc6[v] > +dfc[v]*dfc[v] )
               fc_L = (real)3.0*cc_C - (real)2.0*fc_R;
            else if ( dfc[v]*dfc6[v] < -dfc[v]*dfc[v] )
               fc_R = (real)3.0*cc_C - (real)2.0*fc_L;


//          3-3. ensure the face-centered variables lie between neighboring cell-centered values
            Min  = ( cc_C < cc_L ) ? cc_C : cc_L;
            Max  = ( cc_C > cc_L ) ? cc_C : cc_L;
            fc_L = ( fc_L > Min  ) ? fc_L : Min;
            fc_L = ( fc_L < Max  ) ? fc_L : Max;

            Min  = ( cc_C < cc_R ) ? cc_C : cc_R;
            Max  = ( cc_C > cc_R ) ? cc_C : cc_R;
            fc_R = ( fc_R > Min  ) ? fc_R : Min;
            fc_R = ( fc_R < Max  ) ? fc_R : Max;

            fc[faceL][v] = fc_L;
            fc[faceR][v] = fc_R;

         } // for (int v=0; v<NCOMP_TOTAL; v++)


//       5. primitive variables --> conserved variables
         real tmp[NCOMP_TOTAL];  // input and output arrays must not overlap for Pri2Con()

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceL][v];
         SRHydro_Pri2Con( tmp, fc[faceL], Gamma );

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceR][v];
         SRHydro_Pri2Con( tmp, fc[faceR], Gamma );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    6. advance the face-centered variables by half time-step for the MHM integrator
      SRHydro_HancockPredict( fc, dt, dh, g_ConVar, idx_cc, MinDens, MinTemp );
#     endif


//    7. store the face-centered values to the output array
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(NOut) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : SRHydro_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//
// Parameter   :  L1            : Element x-1
//                C0            : Element x
//                R1            : Element x+1
//                LR_Limiter    : Slope limiter for the data reconstruction in the MHM/MHM_RP schemes
//                                (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff  : Coefficient of the generalized MinMod limiter
//                Gamma         : Ratio of specific heats
//                                --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ           : Target spatial direction : (0/1/2) --> (x/y/z)
//                                --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                Slope_Limiter : Array to store the output monotonic slope
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_LimitSlope( const real L1[], const real C0[], const real R1[], const LR_Limiter_t LR_Limiter,
                       const real MinMod_Coeff, const real Gamma, const int XYZ, real Slope_Limiter[] )
{
   real Slope_L[NCOMP_TOTAL], Slope_R[NCOMP_TOTAL], Slope_C[NCOMP_TOTAL], Slope_A[NCOMP_TOTAL], Slope_LR;


// evaluate different slopes
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_L[v] = C0[v] - L1[v];
      Slope_R[v] = R1[v] - C0[v];
      Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );
   }

   if ( LR_Limiter == VL_GMINMOD )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         if ( Slope_L[v]*Slope_R[v] > (real)0.0 )
            Slope_A[v] = (real)2.0*Slope_L[v]*Slope_R[v]/( Slope_L[v] + Slope_R[v] );
         else
            Slope_A[v] = (real)0.0;
      }
   }


// apply the slope limiter
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      Slope_LR = Slope_L[v]*Slope_R[v];

      if ( Slope_LR > (real)0.0 )
      {
         switch ( LR_Limiter )
         {
            case VANLEER:     // van-Leer
               Slope_Limiter[v] = (real)2.0*Slope_LR/( Slope_L[v] + Slope_R[v] );
               break;

            case GMINMOD:     // generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            case ALBADA:      // van-Albada
               Slope_Limiter[v] = Slope_LR*( Slope_L[v] + Slope_R[v] ) /
                                  ( Slope_L[v]*Slope_L[v] + Slope_R[v]*Slope_R[v] );
               break;

            case VL_GMINMOD:  // van-Leer + generalized MinMod
               Slope_L[v] *= MinMod_Coeff;
               Slope_R[v] *= MinMod_Coeff;
               Slope_Limiter[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_C[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v]  = FMIN(  FABS( Slope_A[v] ), Slope_Limiter[v]  );
               Slope_Limiter[v] *= SIGN( Slope_C[v] );
               break;

            default :
#              ifdef GAMER_DEBUG
               printf( "ERROR : incorrect parameter %s = %d !!\n", "LR_Limiter", LR_Limiter );
#              endif
               return;
         }
      } // if ( Slope_LR > (real)0.0 )

      else
      {
         Slope_Limiter[v] = (real)0.0;
      } // if ( Slope_LR > (real)0.0 ) ... else ...
   } // for (int v=0; v<NCOMP_TOTAL; v++)

} // FUNCTION : SRHydro_LimitSlope



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_HancockPredict
// Description :  Evolve the face-centered variables by half time-step by calculating the face-centered fluxes
//                (no Riemann solver is required)
//
// Note        :  1. Work for the MHM scheme
//                2. Do NOT require data in the neighboring cells
//                3. Input variables must be conserved variables
//
// Parameter   :  fc           : Face-centered conserved variables to be updated
//                dt           : Time interval to advance solution
//                dh           : Cell size
//                g_cc_array   : Array storing the cell-centered conserved variables for checking
//                               negative density and pressure
//                               --> It is just the input array Flu_Array_In[]
//                cc_idx       : Index for accessing g_cc_array[]
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                           const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                           const real MinDens, const real MinTemp, const real Gamma )
{

   const real dt_dh2 = (real)0.5*dt/dh;

   real Flux[6][NCOMP_TOTAL], dFlux;


// calculate flux
   for (int f=0; f<6; f++)    SRHydro_Con2Flux( f/2, Flux[f], fc[f], Gamma, MinTemp );

// update the face-centered variables
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      dFlux = dt_dh2*( Flux[1][v] - Flux[0][v] + Flux[3][v] - Flux[2][v] + Flux[5][v] - Flux[4][v] );

      for (int f=0; f<6; f++)    fc[f][v] -= dFlux;
   }

// check the negative density and energy
   for (int f=0; f<6; f++)
   {
      if ( SRHydro_CheckUnphysical( fc[f], NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, false ) )
      {
//       set to the cell-centered values before update
         for (int f=0; f<6; f++)
         for (int v=0; v<NCOMP_TOTAL; v++)
            fc[f][v] = g_cc_array[v][cc_idx];

         break;
      }
   }

// ensure positive density and pressure
#  ifdef CHECK_NEGATIVE_IN_FLUID
   for (int f=0; f<6; f++)
   {
    SRHydro_CheckUnphysical( fc[f], NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true );
   }
#  endif

} // FUNCTION : SRHydro_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



#endif // #if (  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



#endif // #ifndef __CUFLU_DATARECONSTRUCTION__
