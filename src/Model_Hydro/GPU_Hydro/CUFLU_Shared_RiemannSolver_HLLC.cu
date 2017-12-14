#ifndef __CUFLU_RIEMANNSOLVER_HLLC_CU__
#define __CUFLU_RIEMANNSOLVER_HLLC_CU__



#include "CUFLU.h"
#include "CUFLU_Shared_FluUtility.cu"

static __device__ FluVar CUFLU_RiemannSolver_HLLC( const int XYZ, const FluVar L_In, const FluVar R_In,
                                                   const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                The wave speed is estimated by the same formula in HLLE solver
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput.,
//                            18, 1553
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//                4. The "__noinline__" qualifier is added in Fermi GPUs when "CHECK_INTERMEDIATE == HLLC"
//                   is adopted for higher performance
//
// Parameter   :  XYZ     : Target spatial direction : (0/1/2) --> (x/y/z)
//                L_In    : Input left  state (conserved variables)
//                R_In    : Input right state (conserved variables)
//                Gamma   : Ratio of specific heats
//                MinPres : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
#if ( __CUDA_ARCH__ >= 200  &&  RSOLVER != HLLC ) // for CHECK_INTERMEDIATE == HLLC
__noinline__
#endif
__device__ FluVar CUFLU_RiemannSolver_HLLC( const int XYZ, const FluVar L_In, const FluVar R_In,
                                            const real Gamma, const real MinPres )
{

// 1. reorder the input variables for different spatial directions
   FluVar L = CUFLU_Rotate3D( L_In, XYZ, true );
   FluVar R = CUFLU_Rotate3D( R_In, XYZ, true );


// 2. evaluate the Roe's average values
   const real Gamma_m1 = Gamma - (real)1.0;
   const real  TempRho = (real)0.5*( L.Rho + R.Rho );
   const real _TempRho = (real)1.0/TempRho;

   real _RhoL, _RhoR, P_L, P_R, H_L, H_R, u, v, w, V2, H, Cs;
   real RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum, GammaP_Rho, TempPres;

   _RhoL = (real)1.0 / L.Rho;
   _RhoR = (real)1.0 / R.Rho;
   P_L   = Gamma_m1*(  L.Egy - (real)0.5*( L.Px*L.Px + L.Py*L.Py + L.Pz*L.Pz )*_RhoL  );
   P_R   = Gamma_m1*(  R.Egy - (real)0.5*( R.Px*R.Px + R.Py*R.Py + R.Pz*R.Pz )*_RhoR  );
   P_L   = CUFLU_CheckMinPres( P_L, MinPres );
   P_R   = CUFLU_CheckMinPres( P_R, MinPres );
   H_L   = ( L.Egy + P_L )*_RhoL;
   H_R   = ( R.Egy + P_R )*_RhoR;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(L.Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L.Rho, __FILE__, __LINE__, __FUNCTION__ );

   if ( CUFLU_CheckNegative(R.Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R.Rho, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   RhoL_sqrt       = SQRT( L.Rho );
   RhoR_sqrt       = SQRT( R.Rho );
   _RhoL_sqrt      = (real)1.0 / RhoL_sqrt;
   _RhoR_sqrt      = (real)1.0 / RhoR_sqrt;
   _RhoLR_sqrt_sum = (real)1.0 / (RhoL_sqrt + RhoR_sqrt);

   u  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L.Px + _RhoR_sqrt*R.Px );
   v  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L.Py + _RhoR_sqrt*R.Py );
   w  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L.Pz + _RhoR_sqrt*R.Pz );
   V2 = u*u + v*v + w*w;
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

   GammaP_Rho = Gamma_m1*( H - (real)0.5*V2 );
   TempPres   = GammaP_Rho*TempRho/Gamma;
   TempPres   = CUFLU_CheckMinPres( TempPres, MinPres );
   GammaP_Rho = Gamma*TempPres*_TempRho;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(GammaP_Rho) )
      printf( "ERROR : negative GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Cs = SQRT( GammaP_Rho );


// 3. estimate the maximum wave speeds
   FluVar5 EVal;
   real u_L, u_R, Cs_L, Cs_R, W_L, W_R, MaxV_L, MaxV_R;

   EVal.Rho = u - Cs;
   EVal.Px  = u;
   EVal.Py  = u;
   EVal.Pz  = u;
   EVal.Egy = u + Cs;

   u_L    = _RhoL*L.Px;
   u_R    = _RhoR*R.Px;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(P_L) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_L, __FILE__, __LINE__, __FUNCTION__ );

   if ( CUFLU_CheckNegative(P_R) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_R, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Cs_L   = SQRT( Gamma*P_L*_RhoL );
   Cs_R   = SQRT( Gamma*P_R*_RhoR );
   W_L    = FMIN( EVal.Rho, u_L-Cs_L );
   W_R    = FMAX( EVal.Egy, u_R+Cs_R );
   MaxV_L = FMIN( W_L, (real)0.0 );
   MaxV_R = FMAX( W_R, (real)0.0 );


// 4. evaluate the star-region velocity (V_S) and pressure (P_S)
   real V_S, P_S, temp1_L, temp1_R, temp2_L, temp2_R, temp3;

// take care with large round-off errors when Cs_L<<u_L && Cs_R<<u_R && EVal.Rho>u_L && EVal.Egy<u_R
// ==> temp1_L ~ temp1_R ~ 0.0 ==> temp3 = inf
   /*
   temp1_L = L.Rho*( u_L - W_L );
   temp1_R = R.Rho*( u_R - W_R );
   */
   temp1_L = L.Rho*(  (EVal.Rho<u_L-Cs_L) ? (u_L-EVal.Rho) : (+Cs_L)  );
   temp1_R = R.Rho*(  (EVal.Egy>u_R+Cs_R) ? (u_R-EVal.Egy) : (-Cs_R)  );

   temp2_L = P_L + temp1_L*u_L;
   temp2_R = P_R + temp1_R*u_R;
   temp3   = real(1.0) / ( temp1_L - temp1_R );

   V_S = temp3*( P_L - P_R + temp1_L*u_L - temp1_R*u_R );
   P_S = temp3*( temp1_L*temp2_R - temp1_R*temp2_L );
   P_S = CUFLU_CheckMinPres( P_S, MinPres );


// 5. evaluate the weightings of the left(right) fluxes and contact wave
   FluVar Flux_LR;
   real   temp4, Coeff_S;

   if ( V_S >= (real)0.0 )
   {
      Flux_LR = CUFLU_Con2Flux( L, Gamma_m1, 0, MinPres );

//    fluxes along the maximum wave speed
      Flux_LR.Rho -= MaxV_L*L.Rho;
      Flux_LR.Px  -= MaxV_L*L.Px;
      Flux_LR.Py  -= MaxV_L*L.Py;
      Flux_LR.Pz  -= MaxV_L*L.Pz;
      Flux_LR.Egy -= MaxV_L*L.Egy;

      temp4   = (real)1.0 / ( V_S - MaxV_L );
      Coeff_S = -temp4*MaxV_L*P_S;
   }

   else // V_S < 0.0
   {
      Flux_LR = CUFLU_Con2Flux( R, Gamma_m1, 0, MinPres );

//    fluxes along the maximum wave speed
      Flux_LR.Rho -= MaxV_R*R.Rho;
      Flux_LR.Px  -= MaxV_R*R.Px;
      Flux_LR.Py  -= MaxV_R*R.Py;
      Flux_LR.Pz  -= MaxV_R*R.Pz;
      Flux_LR.Egy -= MaxV_R*R.Egy;

      temp4   = (real)1.0 / ( V_S - MaxV_R );
      Coeff_S = -temp4*MaxV_R*P_S;
   }


// 6. evaluate the HLLC fluxes
   const real Coeff_LR = temp4*V_S;

   Flux_LR.Rho *= Coeff_LR;
   Flux_LR.Px  *= Coeff_LR;
   Flux_LR.Py  *= Coeff_LR;
   Flux_LR.Pz  *= Coeff_LR;
   Flux_LR.Egy *= Coeff_LR;

   Flux_LR.Px  += Coeff_S;
   Flux_LR.Egy += Coeff_S*V_S;


// 7. evaluate the fluxes for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_LR.Rho >= (real)0.0 )
   {
      const real vx = Flux_LR.Rho*_RhoL;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_LR.Passive[v] = L_In.Passive[v]*vx;
   }

   else
   {
      const real vx = Flux_LR.Rho*_RhoR;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_LR.Passive[v] = R_In.Passive[v]*vx;
   }
#  endif


// 8. restore the correct order
   Flux_LR = CUFLU_Rotate3D( Flux_LR, XYZ, false );

   return Flux_LR;

} // FUNCTION : CUFLU_RiemannSolver_HLLC



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLC_CU__
