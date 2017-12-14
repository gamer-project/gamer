#ifndef __CUFLU_RIEMANNSOLVER_HLLE_CU__
#define __CUFLU_RIEMANNSOLVER_HLLE_CU__



#include "CUFLU.h"
#include "CUFLU_Shared_FluUtility.cu"

static __device__ FluVar CUFLU_RiemannSolver_HLLE( const int XYZ, const FluVar L_In, const FluVar R_In,
                                                   const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemannSolver_HLLE
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                Estimate the wave speed by Einfeldt et al. (1991).
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         b. Einfeldt, B., et al. J. 1991, J. Comput. Phys., 92, 273
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//                4. The "__noinline__" qualifier is added in Fermi GPUs when "CHECK_INTERMEDIATE == HLLE"
//                   is adopted for higher performance
//
// Parameter   :  XYZ     : Target spatial direction : (0/1/2) --> (x/y/z)
//                L_In    : Input left  state (conserved variables)
//                R_In    : Input right state (conserved variables)
//                Gamma   : Ratio of specific heats
//                MinPres : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
#if ( __CUDA_ARCH__ >= 200  &&  RSOLVER != HLLE ) // for CHECK_INTERMEDIATE == HLLE
__noinline__
#endif
__device__ FluVar CUFLU_RiemannSolver_HLLE( const int XYZ, const FluVar L_In, const FluVar R_In,
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
   real    u_L, u_R, Cs_L, Cs_R, MaxV_L, MaxV_R;

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
   MaxV_L = FMIN( EVal.Rho, u_L-Cs_L );
   MaxV_R = FMAX( EVal.Egy, u_R+Cs_R );
   MaxV_L = FMIN( MaxV_L, (real)0.0 );
   MaxV_R = FMAX( MaxV_R, (real)0.0 );


// 4. evaluate the left and right fluxes along the maximum wave speeds
   FluVar Flux_L = CUFLU_Con2Flux( L, Gamma_m1, 0, MinPres );
   FluVar Flux_R = CUFLU_Con2Flux( R, Gamma_m1, 0, MinPres );

   Flux_L.Rho -= MaxV_L*L.Rho;
   Flux_L.Px  -= MaxV_L*L.Px;
   Flux_L.Py  -= MaxV_L*L.Py;
   Flux_L.Pz  -= MaxV_L*L.Pz;
   Flux_L.Egy -= MaxV_L*L.Egy;

   Flux_R.Rho -= MaxV_R*R.Rho;
   Flux_R.Px  -= MaxV_R*R.Px;
   Flux_R.Py  -= MaxV_R*R.Py;
   Flux_R.Pz  -= MaxV_R*R.Pz;
   Flux_R.Egy -= MaxV_R*R.Egy;


// 5. evaluate the HLLE fluxes
   const real _MaxV_R_minus_L = (real)1.0 / ( MaxV_R - MaxV_L );
   FluVar Flux_Out;

   Flux_Out.Rho = _MaxV_R_minus_L*( MaxV_R*Flux_L.Rho - MaxV_L*Flux_R.Rho );
   Flux_Out.Px  = _MaxV_R_minus_L*( MaxV_R*Flux_L.Px  - MaxV_L*Flux_R.Px  );
   Flux_Out.Py  = _MaxV_R_minus_L*( MaxV_R*Flux_L.Py  - MaxV_L*Flux_R.Py  );
   Flux_Out.Pz  = _MaxV_R_minus_L*( MaxV_R*Flux_L.Pz  - MaxV_L*Flux_R.Pz  );
   Flux_Out.Egy = _MaxV_R_minus_L*( MaxV_R*Flux_L.Egy - MaxV_L*Flux_R.Egy );


// 6. evaluate the fluxes for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out.Rho >= (real)0.0 )
   {
      const real vx = Flux_Out.Rho*_RhoL;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_Out.Passive[v] = L_In.Passive[v]*vx;
   }

   else
   {
      const real vx = Flux_Out.Rho*_RhoR;

      for (int v=0; v<NCOMP_PASSIVE; v++)    Flux_Out.Passive[v] = R_In.Passive[v]*vx;
   }
#  endif


// 7. restore the correct order
   Flux_Out = CUFLU_Rotate3D( Flux_Out, XYZ, false );

   return Flux_Out;

} // FUNCTION : CUFLU_RiemannSolver_HLLE



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLE_CU__
