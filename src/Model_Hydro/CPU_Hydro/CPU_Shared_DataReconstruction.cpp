#ifndef __CUFLU_DATARECONSTRUCTION__
#define __CUFLU_DATARECONSTRUCTION__



#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward );
real Hydro_CheckMinPres( const real InPres, const real MinPres );
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                    const bool NormPassive, const int NNorm, const int NormIdx[] );
#if ( FLU_SCHEME == MHM )
void Hydro_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres );
#endif

#endif // #ifdef __CUDACC__ ... else ...


// internal functions (GPU_DEVICE is defined in CUFLU.h)
GPU_DEVICE
static void Hydro_LimitSlope( const real L1[], const real C0[], const real R1[], const LR_Limiter_t LR_Limiter,
                              const real MinMod_Coeff, const real Gamma, const int XYZ, real Slope_Limiter[] );
#if ( FLU_SCHEME == CTU )
GPU_DEVICE
static void Hydro_GetEigenSystem( const real CC_Var[], real EigenVal[][NCOMP_FLUID], real LEigenVec[][NCOMP_FLUID],
                                  real REigenVec[][NCOMP_FLUID], const real Gamma );
#endif
#if ( FLU_SCHEME == MHM )
GPU_DEVICE
static void Hydro_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                                  const real Gamma_m1, const real _Gamma_m1,
                                  const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                                  const real MinDens, const real MinPres );
#endif
#ifdef CHAR_RECONSTRUCTION
GPU_DEVICE
static void Hydro_Pri2Char( real Var[], const real Gamma, const real Rho, const real Pres, const int XYZ );
GPU_DEVICE
static void Hydro_Char2Pri( real Var[], const real Gamma, const real Rho, const real Pres, const int XYZ );
#endif




#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. Input data can be either conserved or primitive variables
//                   --> If the input data are conserved variables, one must provide g_ConVar[] and enable "Con2Pri"
//                       --> Primitive variables will be calculated by this function and stored in g_PriVar[]
//                       --> g_PriVar[] must be allocated in advance but it doesn't need to be initialized
//                       --> Adopted by MHM and CTU
//                   --> If the input data are primitive variables, one must provide g_PriVar[] and disable "Con2Pri"
//                       --> g_ConVar[] is useless here
//                           --> Although Hydro_HancockPredict() still needs g_ConVar(), MHM provides conserved instead
//                               of primitive variables. So it's fine.
//                       --> Adopted by MHM_RP, where Hydro_RiemannPredict() already returns primitive variables
//                3. Output data are always conserved variables
//                   --> Because Hydro_HancockPredict() only works with conserved variables
//                4. PLM and PPM data reconstruction functions share the same function name
//                5. Face-centered variables will be advanced by half time-step for the MHM and CTU schemes
//                6. Data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION" in the header CUFLU.h
//                7. This function is shared by MHM, MHM_RP, and CTU schemes
//
//
// Parameter   :  g_ConVar           : Array storing the input cell-centered conserved variables
//                g_PriVar           : Array storing/to store the cell-centered primitive variables
//                                     --> For MHM, g_ConVar[] and g_PriVar[] must point to different arrays since
//                                         Hydro_HancockPredict() requires the original g_ConVar[]
//                g_FC_Var           : Array to store the output face-centered primitive variables
//                g_Slope_PPM        : Array to store the x/y/z slopes for the PPM reconstruction
//                                     --> Useless for PLM
//                Con2Pri            : Convert conserved variables in g_ConVar[] to primitive variables and
//                                     store the results in g_PriVar[]
//                NIn                : Size of g_PriVar[] along each direction
//                                     --> Can be smaller than FLU_NXT
//                NGhost             : Number of ghost zones
//                                      --> "NIn-2*NGhost" cells will be computed along each direction
//                                      --> Size of g_FC_Var[] is assumed to be "(NIn-2*NGhost)^3"
//                                      --> The reconstructed data at cell (i,j,k) will be stored in g_FC_Var[]
//                                          with the index "(i-NGhost,j-NGhost,k-NGhost)"
//                Gamma              : Ratio of specific heats
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                dt                 : Time interval to advance solution (for the CTU scheme)
//                dh                 : Cell size
//                MinDens/Pres       : Minimum allowed density and pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh, const real MinDens, const real MinPres,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const int  didx_cc[3] = { 1, NIn, SQR(NIn) };
   const int  NOut       = NIn - 2*NGhost;      // number of output cells
   const real  Gamma_m1  = Gamma - (real)1.0;
   const real _Gamma_m1  = (real)1.0 / Gamma_m1;

#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt/dh;

// constant components of the left and right eigenvector matrices must be initialized
   real LEigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #if ( FLU_SCHEME ==  CTU )


// 0. conserved --> primitive variables
   if ( Con2Pri )
   {
      real ConVar_1Cell[NCOMP_TOTAL], PriVar_1Cell[NCOMP_TOTAL];

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

         Hydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                        JeansMinPres, JeansMinPres_Coeff );

         for (int v=0; v<NCOMP_TOTAL; v++)   g_PriVar[v][idx] = PriVar_1Cell[v];
      }

#     ifdef __CUDACC__
      __syncthreads();
#     endif
   } // if ( Con2Pri )


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
#     if ( FLU_SCHEME == CTU )
      real EigenVal[3][NCOMP_FLUID], Correct_L[NCOMP_TOTAL], Correct_R[NCOMP_TOTAL], dfc[NCOMP_TOTAL];
      real Coeff_L, Coeff_R;
#     endif

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C[v] = g_PriVar[v][idx_cc];


//    1. evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Hydro_GetEigenSystem( cc_C, EigenVal, LEigenVec, REigenVec, Gamma );
#     endif


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

         Hydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, Gamma, d, Slope_Limiter );


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


//       4. advance the face-centered variables by half time-step for the CTU integrator
#        if ( FLU_SCHEME == CTU )

//       4-1. evaluate the slope (for passive scalars as well)
         for (int v=0; v<NCOMP_TOTAL; v++)   dfc[v] = fc[faceR][v] - fc[faceL][v];


//       4-2. re-order variables for the y/z directions
         Hydro_Rotate3D( dfc, d, true );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

//       4-2-a1. evaluate the corrections to the left and right face-centered variables

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = (real)0.0;
            Correct_R[v] = (real)0.0;
         }

#        ifdef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;

            for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[v];

            Coeff_L *= -dt_dh2*EigenVal[d][Mode];

            for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

         for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] = Correct_L[v];

#        else // #ifdef HLL_INCLUDE_ALL_WAVES

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[v];

               Coeff_L *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*dfc[v];

               Coeff_R *= -dt_dh2*EigenVal[d][Mode];

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // #ifdef HLL_INCLUDE_ALL_WAVES ... else ...


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT || ifndef HLL_NO_REF_STATE )

//       4-2-b1. evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][0], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][4], (real)0.0 );

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = Coeff_L*dfc[v];
            Correct_R[v] = Coeff_R*dfc[v];
         }


//       4-2-b2. evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*dfc[v];

               Coeff_L *= dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*dfc[v];

               Coeff_R *= dt_dh2*( EigenVal[d][4] - EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       4-2-b3. evaluate the corrections to the left and right face-centered passive scalars
//               --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*dfc[v];
            Correct_R[v] = Coeff_R*dfc[v];
         }
#        endif


//       4-2-b4. evaluate the face-centered variables at the half time-step
         Hydro_Rotate3D( Correct_L, d, false );
         Hydro_Rotate3D( Correct_R, d, false );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fc[faceL][v] += Correct_L[v];
            fc[faceR][v] += Correct_R[v];
         }

//       ensure positive density and pressure
         fc[faceL][0] = FMAX( fc[faceL][0], MinDens );
         fc[faceR][0] = FMAX( fc[faceR][0], MinDens );

         fc[faceL][4] = Hydro_CheckMinPres( fc[faceL][4], MinPres );
         fc[faceR][4] = Hydro_CheckMinPres( fc[faceR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         fc[faceL][v] = FMAX( fc[faceL][v], TINY_NUMBER );
         fc[faceR][v] = FMAX( fc[faceR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )


//       5. primitive variables --> conserved variables
         real tmp[NCOMP_TOTAL];  // input and output arrays must not overlap for Pri2Con()

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceL][v];
         Hydro_Pri2Con( tmp, fc[faceL], _Gamma_m1, NormPassive, NNorm, NormIdx );

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceR][v];
         Hydro_Pri2Con( tmp, fc[faceR], _Gamma_m1, NormPassive, NNorm, NormIdx );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    6. advance the face-centered variables by half time-step for the MHM integrator
      Hydro_HancockPredict( fc, dt, dh, Gamma_m1, _Gamma_m1, g_ConVar, idx_cc, MinDens, MinPres );
#     endif


//    7. store the face-centered values to the output array
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(NOut) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  See the PLM routine
//
// Parameter   :  See the PLM routine
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh, const real MinDens, const real MinPres,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const int  didx_cc   [3] = { 1, NIn, SQR(NIn) };
   const int  didx_slope[3] = { 1, N_SLOPE_PPM, SQR(N_SLOPE_PPM) };
   const int  NOut       = NIn - 2*NGhost;      // number of output cells
   const int  NSlope     = N_SLOPE_PPM;         // size of g_Slope_PPM[] (which must be equal to NOut + 2)
   const real  Gamma_m1  = Gamma - (real)1.0;
   const real _Gamma_m1  = (real)1.0 / Gamma_m1;


#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt/dh;

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )
#  ifdef HLL_INCLUDE_ALL_WAVES
   const bool HLL_Include_All_Waves = true;
#  else
   const bool HLL_Include_All_Waves = false;
#  endif
#  endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

// constant components of the left and right eigenvector matrices must be initialized
   real LEigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   real REigenVec[NCOMP_FLUID][NCOMP_FLUID] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                                { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                                { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                                { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
#  endif // #if ( FLU_SCHEME == CTU )


// 0. conserved --> primitive variables
   if ( Con2Pri )
   {
      real ConVar_1Cell[NCOMP_TOTAL], PriVar_1Cell[NCOMP_TOTAL];

      CGPU_LOOP( idx, CUBE(NIn) )
      {
         for (int v=0; v<NCOMP_TOTAL; v++)   ConVar_1Cell[v] = g_ConVar[v][idx];

         Hydro_Con2Pri( ConVar_1Cell, PriVar_1Cell, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                        JeansMinPres, JeansMinPres_Coeff );

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

         Hydro_LimitSlope( cc_L, cc_C, cc_R, LR_Limiter, MinMod_Coeff, Gamma, d, Slope_Limiter );

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
#     if ( FLU_SCHEME == CTU )
      real EigenVal[3][NCOMP_FLUID], Correct_L[NCOMP_TOTAL], Correct_R[NCOMP_TOTAL];
      real Coeff_L, Coeff_R;
#     endif

      for (int v=0; v<NCOMP_TOTAL; v++)   cc_C_ncomp[v] = g_PriVar[v][idx_cc];


//    2. evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Hydro_GetEigenSystem( cc_C_ncomp, EigenVal, LEigenVec, REigenVec, Gamma );
#     endif


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


//       4. advance the face-centered variables by half time-step for the CTU integrator
#        if ( FLU_SCHEME == CTU )

//       4-1. compute the PPM coefficient (for the passive scalars as well)
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            dfc [v] = fc[faceR][v] - fc[faceL][v];
            dfc6[v] = (real)6.0*(  cc_C_ncomp[v] - (real)0.5*( fc[faceL][v] + fc[faceR][v] )  );
         }

//       4-2. re-order variables for the y/z directions
         Hydro_Rotate3D( dfc,  d, true );
         Hydro_Rotate3D( dfc6, d, true );


//       =====================================================================================
//       a. for the HLL solvers (HLLE/HLLC)
//       =====================================================================================
#        if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

//       4-2-a1. evaluate the corrections to the left and right face-centered variables
         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = (real)0.0;
            Correct_R[v] = (real)0.0;
         }

         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] <= (real)0.0 )
            {
               const real Coeff_C = -dt_dh2*EigenVal[d][Mode];
               const real Coeff_D = real(-4.0/3.0)*SQR(Coeff_C);

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[v] + dfc6[v] ) +
                                                                                    Coeff_D*( dfc6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] >= (real)0.0 )
            {
               const real Coeff_A = -dt_dh2*EigenVal[d][Mode];
               const real Coeff_B = real(-4.0/3.0)*SQR(Coeff_A);

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[v] - dfc6[v] ) +
                                                                                    Coeff_B*( dfc6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)


//       =====================================================================================
//       b. for the Roe's and exact solvers
//       =====================================================================================
#        else // ( RSOLVER == ROE/EXACT || ifndef HLL_NO_REF_STATE )

//       4-2-b1. evaluate the reference states
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][0], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][4], (real)0.0 );

         for (int v=0; v<NCOMP_FLUID; v++)
         {
            Correct_L[v] = Coeff_L*(  dfc[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[v]  );
            Correct_R[v] = Coeff_R*(  dfc[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[v]  );
         }


//       4-2-b2. evaluate the corrections to the left and right face-centered variables
         for (int Mode=0; Mode<NCOMP_FLUID; Mode++)
         {
            Coeff_L = (real)0.0;
            Coeff_R = (real)0.0;

            if ( EigenVal[d][Mode] <= (real)0.0 )
            {
               const real Coeff_C = dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_D=0 when Coeff_C=0
//             Coeff_D = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][   0]*EigenVal[d][   0] -
//                                                      EigenVal[d][Mode]*EigenVal[d][Mode]   );
               const real Coeff_D = real(4.0/3.0)*dt_dh2*Coeff_C*( EigenVal[d][0] + EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[v] + dfc6[v] ) +
                                                                                    Coeff_D*( dfc6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_L[v] += Coeff_L*REigenVec[Mode][v];
            }

            if ( EigenVal[d][Mode] >= (real)0.0 )
            {
               const real Coeff_A = dt_dh2*( EigenVal[d][4] - EigenVal[d][Mode] );
//             write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_B=0 when Coeff_A=0
//             Coeff_B = real(4.0/3.0)*dt_dh2*dt_dh2* ( EigenVal[d][   4]*EigenVal[d][   4] -
//                                                      EigenVal[d][Mode]*EigenVal[d][Mode]   );
               const real Coeff_B = real(4.0/3.0)*dt_dh2*Coeff_A*( EigenVal[d][4] + EigenVal[d][Mode] );

               for (int v=0; v<NCOMP_FLUID; v++)   Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[v] - dfc6[v] ) +
                                                                                    Coeff_B*( dfc6[v]          )   );

               for (int v=0; v<NCOMP_FLUID; v++)   Correct_R[v] += Coeff_R*REigenVec[Mode][v];
            }
         } // for (int Mode=0; Mode<NCOMP_FLUID; Mode++)

#        endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//       4-3. evaluate the corrections to the left and right face-centered passive scalars
//            --> passive scalars travel with fluid velocity (i.e., entropy mode)
#        if ( NCOMP_PASSIVE > 0 )
         Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
         Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Correct_L[v] = Coeff_L*(  dfc[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[v]  );
            Correct_R[v] = Coeff_R*(  dfc[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[v]  );
         }
#        endif


//       4-4. evaluate the face-centered variables at the half time-step
         Hydro_Rotate3D( Correct_L, d, false );
         Hydro_Rotate3D( Correct_R, d, false );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            fc[faceL][v] += Correct_L[v];
            fc[faceR][v] += Correct_R[v];
         }

//       ensure positive density and pressure
         fc[faceL][0] = FMAX( fc[faceL][0], MinDens );
         fc[faceR][0] = FMAX( fc[faceR][0], MinDens );

         fc[faceL][4] = Hydro_CheckMinPres( fc[faceL][4], MinPres );
         fc[faceR][4] = Hydro_CheckMinPres( fc[faceR][4], MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         fc[faceL][v] = FMAX( fc[faceL][v], TINY_NUMBER );
         fc[faceR][v] = FMAX( fc[faceR][v], TINY_NUMBER ); }
#        endif

#        endif // #if ( FLU_SCHEME == CTU )


//       5. primitive variables --> conserved variables
         real tmp[NCOMP_TOTAL];  // input and output arrays must not overlap for Pri2Con()

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceL][v];
         Hydro_Pri2Con( tmp, fc[faceL], _Gamma_m1, NormPassive, NNorm, NormIdx );

         for (int v=0; v<NCOMP_TOTAL; v++)   tmp[v] = fc[faceR][v];
         Hydro_Pri2Con( tmp, fc[faceR], _Gamma_m1, NormPassive, NNorm, NormIdx );

      } // for (int d=0; d<3; d++)


#     if ( FLU_SCHEME == MHM )
//    6. advance the face-centered variables by half time-step for the MHM integrator
      Hydro_HancockPredict( fc, dt, dh, Gamma_m1, _Gamma_m1, g_ConVar, idx_cc, MinDens, MinPres );
#     endif


//    7. store the face-centered values to the output array
      for (int f=0; f<6; f++)
      for (int v=0; v<NCOMP_TOTAL; v++)
         g_FC_Var[f][v][idx_fc] = fc[f][v];

   } // CGPU_LOOP( idx_fc, CUBE(NOut) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )



#ifdef CHAR_RECONSTRUCTION
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Pri2Char
// Description :  Primitive variables --> characteristic variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//                2. Input and output share the same array
//
// Parameter   :  InOut : Array storing both the input primitive variables and output characteristic variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Pri2Char( real InOut[], const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Cs2 = (real)1.0 / ( Gamma*Pres/Rho );
   const real _Cs  = SQRT( _Cs2 );
   real Temp[NCOMP_FLUID];

   for (int v=0; v<NCOMP_FLUID; v++)   Temp[v] = InOut[v];

   Hydro_Rotate3D( Temp, XYZ, true );

   InOut[0] = -(real)0.5*Rho*_Cs*Temp[1] + (real)0.5*_Cs2*Temp[4];
   InOut[1] = Temp[0] - _Cs2*Temp[4];
   InOut[2] = Temp[2];
   InOut[3] = Temp[3];
   InOut[4] = +(real)0.5*Rho*_Cs*Temp[1] + (real)0.5*_Cs2*Temp[4];

} // FUNCTION : Hydro_Pri2Char



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Char2Pri
// Description :  Characteristic variables --> primitive variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//                2. Input and output share the same array
//
// Parameter   :  InOut : Array storing both the input characteristic variables and output primitive variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Char2Pri( real InOut[], const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(Pres) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Rho = (real)1.0 / Rho;
   const real Cs2  = Gamma*Pres*_Rho;
   const real Cs   = SQRT( Cs2 );
   real Temp[NCOMP_FLUID];

   Temp[0] = InOut[0] + InOut[1] + InOut[4];
   Temp[1] = Cs*_Rho*( -InOut[0] + InOut[4] );
   Temp[2] = InOut[2];
   Temp[3] = InOut[3];
   Temp[4] = Cs2*( InOut[0] + InOut[4] );

   Hydro_Rotate3D( Temp, XYZ, false );

   for (int v=0; v<NCOMP_FLUID; v++)   InOut[v] = Temp[v];

} // FUNCTION : Hydro_Char2Pri
#endif



#if ( FLU_SCHEME == CTU )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_GetEigenSystem
// Description :  Evaluate the eigenvalues and left/right eigenvectors
//
// Note        :  1. Input data must be primitive variables
//                2. Constant components of eigenvectors must be set in advance
//                3. Work for the CTU scheme
//                4. Do not need to consider passive scalars
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  CC_Var      : Array storing the input cell-centered primitive variables
//                EigenVal    : Array to store the output eigenvalues (in three spatial directions)
//                L/REigenVec : Array to store the output left/right eigenvectors
//                Gamma       : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_GetEigenSystem( const real CC_Var[], real EigenVal[][NCOMP_FLUID], real LEigenVec[][NCOMP_FLUID],
                           real REigenVec[][NCOMP_FLUID], const real Gamma )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(CC_Var[4]) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              CC_Var[4], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(CC_Var[0]) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              CC_Var[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real  Rho = CC_Var[0];
   const real _Rho = (real)1.0/Rho;
   const real  Vx  = CC_Var[1];
   const real  Vy  = CC_Var[2];
   const real  Vz  = CC_Var[3];
   const real Pres = CC_Var[4];
   const real  Cs2 = Gamma*Pres*_Rho;
   const real  Cs  = SQRT( Cs2 );
   const real _Cs  = (real)1.0/Cs;
   const real _Cs2 = _Cs*_Cs;

// a. eigenvalues in three spatial directions
   EigenVal[0][0] = Vx - Cs;
   EigenVal[0][1] = Vx;
   EigenVal[0][2] = Vx;
   EigenVal[0][3] = Vx;
   EigenVal[0][4] = Vx + Cs;

   EigenVal[1][0] = Vy - Cs;
   EigenVal[1][1] = Vy;
   EigenVal[1][2] = Vy;
   EigenVal[1][3] = Vy;
   EigenVal[1][4] = Vy + Cs;

   EigenVal[2][0] = Vz - Cs;
   EigenVal[2][1] = Vz;
   EigenVal[2][2] = Vz;
   EigenVal[2][3] = Vz;
   EigenVal[2][4] = Vz + Cs;


// NOTE : the left and right eigenvectors have the same form in different spatial directions
// b. left eigenvectors (rows of the matrix LEigenVec)
   LEigenVec[0][1] = -(real)0.5*Rho*_Cs;
   LEigenVec[0][4] = (real)0.5*_Cs2;
   LEigenVec[1][4] = -_Cs2;
   LEigenVec[4][1] = -LEigenVec[0][1];
   LEigenVec[4][4] = +LEigenVec[0][4];


// c. right eigenvectors (rows of the matrix REigenVec)
   REigenVec[0][1] = -Cs*_Rho;
   REigenVec[0][4] = Cs2;
   REigenVec[4][1] = -REigenVec[0][1];
   REigenVec[4][4] = Cs2;

} // FUNCTION : Hydro_GetEigenSystem
#endif // #if ( FLU_SCHEME == CTU )



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//
// Parameter   :  L1            : Element x-1
//                C0            : Element x
//                R1            : Element x+1
//                LR_Limiter    : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                (0/1/2/3) = (vanLeer/generalized MinMod/vanAlbada/vanLeer+generalized MinMod) limiter
//                MinMod_Coeff  : Coefficient of the generalized MinMod limiter
//                Gamma         : Ratio of specific heats
//                                --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ           : Target spatial direction : (0/1/2) --> (x/y/z)
//                                --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                Slope_Limiter : Array to store the output monotonic slope
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_LimitSlope( const real L1[], const real C0[], const real R1[], const LR_Limiter_t LR_Limiter,
                       const real MinMod_Coeff, const real Gamma, const int XYZ, real Slope_Limiter[] )
{

#  ifdef CHAR_RECONSTRUCTION
   const real Rho  = C0[0];
   const real Pres = C0[4];
#  endif

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


// primitive variables --> characteristic variables
#  ifdef CHAR_RECONSTRUCTION
   Hydro_Pri2Char( Slope_L, Gamma, Rho, Pres, XYZ );
   Hydro_Pri2Char( Slope_R, Gamma, Rho, Pres, XYZ );
   Hydro_Pri2Char( Slope_C, Gamma, Rho, Pres, XYZ );

   if ( LR_Limiter == VL_GMINMOD )
      Hydro_Pri2Char( Slope_A, Gamma, Rho, Pres, XYZ );
#  endif


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


// characteristic variables --> primitive variables
#  ifdef CHAR_RECONSTRUCTION
   Hydro_Char2Pri( Slope_Limiter, Gamma, Rho, Pres, XYZ );
#  endif

} // FUNCTION : Hydro_LimitSlope



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_HancockPredict
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
//                Gamma_m1     : Gamma - 1
//                _Gamma_m1    : 1 / (Gamma - 1)
//                g_cc_array   : Array storing the cell-centered conserved variables for checking
//                               negative density and pressure
//                               --> It is just the input array Flu_Array_In[]
//                cc_idx       : Index for accessing g_cc_array[]
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_HancockPredict( real fc[][NCOMP_TOTAL], const real dt, const real dh,
                           const real Gamma_m1, const real _Gamma_m1,
                           const real g_cc_array[][ CUBE(FLU_NXT) ], const int cc_idx,
                           const real MinDens, const real MinPres )
{

   const real dt_dh2 = (real)0.5*dt/dh;

   real Flux[6][NCOMP_TOTAL], dFlux;


// calculate flux
   for (int f=0; f<6; f++)    Hydro_Con2Flux( f/2, Flux[f], fc[f], Gamma_m1, MinPres );

// update the face-centered variables
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      dFlux = dt_dh2*( Flux[1][v] - Flux[0][v] + Flux[3][v] - Flux[2][v] + Flux[5][v] - Flux[4][v] );

      for (int f=0; f<6; f++)    fc[f][v] -= dFlux;
   }

// check the negative density and energy
   for (int f=0; f<6; f++)
   {
      if ( fc[f][0] <= (real)0.0  ||  fc[f][4] <= (real)0.0 )
      {
//       set to the cell-centered values before update
         for (int f=0; f<6; f++)
         for (int v=0; v<NCOMP_TOTAL; v++)
            fc[f][v] = g_cc_array[v][cc_idx];

         break;
      }
   }

// ensure positive density and pressure
   for (int f=0; f<6; f++)
   {
      fc[f][0] = FMAX( fc[f][0], MinDens );
      fc[f][4] = Hydro_CheckMinPresInEngy( fc[f][0], fc[f][1], fc[f][2], fc[f][3], fc[f][4],
                                           Gamma_m1, _Gamma_m1, MinPres );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      fc[f][v] = FMAX( fc[f][v], TINY_NUMBER );
#     endif
   }

} // FUNCTION : Hydro_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )



#endif // #ifndef __CUFLU_DATARECONSTRUCTION__
