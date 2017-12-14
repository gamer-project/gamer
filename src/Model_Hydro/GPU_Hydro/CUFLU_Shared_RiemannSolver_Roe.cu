#ifndef __CUFLU_RIEMANNSOLVER_ROE_CU__
#define __CUFLU_RIEMANNSOLVER_ROE_CU__



#include "CUFLU.h"
#include "CUFLU_Shared_FluUtility.cu"
#if   ( CHECK_INTERMEDIATE == EXACT )
#include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( CHECK_INTERMEDIATE == HLLE )
#include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( CHECK_INTERMEDIATE == HLLC )
#include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

static __device__ FluVar CUFLU_RiemannSolver_Roe( const int XYZ, const FluVar L_In, const FluVar R_In,
                                                  const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemannSolver_Roe
// Description :  Approximate Riemann solver of Roe
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : "Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                         ~ by Eleuterio F. Toro"
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//                4. The "__forceinline__" qualifier is added for higher performance
//
// Parameter   :  XYZ         : Target spatial direction : (0/1/2) --> (x/y/z)
//                L_In        : Input left  state (conserved variables)
//                R_In        : Input right state (conserved variables)
//                Gamma       : Gamma
//                MinPres     : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
__forceinline__
__device__ FluVar CUFLU_RiemannSolver_Roe( const int XYZ, const FluVar L_In, const FluVar R_In,
                                           const real Gamma, const real MinPres )
{

// 1. reorder the input variables for different spatial directions
   FluVar L = CUFLU_Rotate3D( L_In, XYZ, true );
   FluVar R = CUFLU_Rotate3D( R_In, XYZ, true );


// 2. evaluate the average values
   const real Gamma_m1 = Gamma - (real)1.0;
   const real  TempRho = (real)0.5*( L.Rho + R.Rho );
   const real _TempRho = (real)1.0/TempRho;
   real _RhoL, _RhoR, HL, HR, u, v, w, V2, H, Cs, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real GammaP_Rho, TempPres;

   _RhoL = (real)1.0 / L.Rho;
   _RhoR = (real)1.0 / R.Rho;
   HL    = (  L.Egy + Gamma_m1*( L.Egy - (real)0.5*( L.Px*L.Px + L.Py*L.Py + L.Pz*L.Pz )*_RhoL )  )*_RhoL;
   HR    = (  R.Egy + Gamma_m1*( R.Egy - (real)0.5*( R.Px*R.Px + R.Py*R.Py + R.Pz*R.Pz )*_RhoR )  )*_RhoR;

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
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*HL   +  RhoR_sqrt*HR   );

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


// 3. evaluate the eigenvalues
   FluVar5 EVal;

   EVal.Rho = u - Cs;
   EVal.Px  = u;
   EVal.Py  = u;
   EVal.Pz  = u;
   EVal.Egy = u + Cs;


// 4. evaluate the left and right fluxes
   const FluVar Flux_L = CUFLU_Con2Flux( L, Gamma_m1, 0, MinPres );
   const FluVar Flux_R = CUFLU_Con2Flux( R, Gamma_m1, 0, MinPres );


// 5. return the upwind fluxes if flow is supersonic
   FluVar Flux_Out;

   if ( EVal.Rho >= (real)0.0 )
   {
      Flux_Out = CUFLU_Rotate3D( Flux_L, XYZ, false );

      return Flux_Out;
   }

   if ( EVal.Egy <= (real)0.0 )
   {
      Flux_Out = CUFLU_Rotate3D( Flux_R, XYZ, false );

      return Flux_Out;
   }


// 6. evaluate the amplitudes along different characteristics (eigenvectors)
   const real _Cs = (real)1.0/Cs;
   FluVar5 Jump, Amp;

   Jump.Rho = R.Rho - L.Rho;
   Jump.Px  = R.Px  - L.Px;
   Jump.Py  = R.Py  - L.Py;
   Jump.Pz  = R.Pz  - L.Pz;
   Jump.Egy = R.Egy - L.Egy;

   Amp.Py  = Jump.Py - v*Jump.Rho;
   Amp.Pz  = Jump.Pz - w*Jump.Rho;
   Amp.Px  = Gamma_m1*_Cs*_Cs*( Jump.Rho*(H-u*u) + u*Jump.Px - Jump.Egy + v*Amp.Py + w*Amp.Pz );
   Amp.Rho = (real)0.5*_Cs*( Jump.Rho*(u+Cs) - Jump.Px - Cs*Amp.Px );
   Amp.Egy = Jump.Rho - Amp.Rho - Amp.Px;


// 7. evaluate the eigenvectors
   const real ONE  = (real)1.0;
   const real ZERO = (real)0.0;
   FluVar5 EVec1, EVec2, EVec3, EVec4, EVec5;

   EVec1.Rho = ONE;    EVec1.Px = u - Cs;   EVec1.Py = v;      EVec1.Pz = w;      EVec1.Egy = H - u*Cs;
   EVec2.Rho = ONE;    EVec2.Px = u;        EVec2.Py = v;      EVec2.Pz = w;      EVec2.Egy = (real)0.5*V2;
   EVec3.Rho = ZERO;   EVec3.Px = ZERO;     EVec3.Py = ONE;    EVec3.Pz = ZERO;   EVec3.Egy = v;
   EVec4.Rho = ZERO;   EVec4.Px = ZERO;     EVec4.Py = ZERO;   EVec4.Pz = ONE;    EVec4.Egy = w;
   EVec5.Rho = ONE;    EVec5.Px = u + Cs;   EVec5.Py = v;      EVec5.Pz = w;      EVec5.Egy = H + u*Cs;


// 8. verify that the density and pressure in the intermediate states are positive
#  ifdef CHECK_INTERMEDIATE
   FluVar5 I_States;
   real    I_Pres;

#  if   ( CHECK_INTERMEDIATE == EXACT )   // recalculate fluxes by exact solver

#     define Recalculate_Flux( L, R, Flux_Out )                                                    \
      {                                                                                            \
         /* do NOT convert any passive variable to mass fraction for the Riemann solvers */        \
         const bool NormPassive_No  = false;                                                       \
         const bool JeansMinPres_No = false;                                                       \
                                                                                                   \
         L = CUFLU_Con2Pri( L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL ); \
         R = CUFLU_Con2Pri( R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL ); \
                                                                                                   \
         Flux_Out = CUFLU_RiemannSolver_Exact( 0, NULL, NULL, NULL, L, R, Gamma );                 \
      } // Recalculate_Flux

#  elif ( CHECK_INTERMEDIATE == HLLE )    // recalculate fluxes by HLLE solver

#     define Recalculate_Flux( L, R, Flux_Out )                                                    \
      {                                                                                            \
         Flux_Out = CUFLU_RiemannSolver_HLLE( 0, L, R, Gamma, MinPres );                           \
      } // Recalculate_Flux

#  elif ( CHECK_INTERMEDIATE == HLLC )    // recalculate fluxes by HLLC solver

#     define Recalculate_Flux( L, R, Flux_Out )                                                    \
      {                                                                                            \
         Flux_Out = CUFLU_RiemannSolver_HLLC( 0, L, R, Gamma, MinPres );                           \
      } // Recalculate_Flux

#  else

#  error : ERROR : unsupported CHECK_INTERMEDIATE (EXACT/HLLE/HLLC) !!

#  endif // CHECK_INTERMEDIATE == EXACT/HLLE/HLLC


#  define Get_I_States( comp, comp_next, EVec )                                                    \
   {                                                                                               \
      I_States.Rho += Amp.comp*EVec.Rho;                                                           \
      I_States.Px  += Amp.comp*EVec.Px;                                                            \
      I_States.Py  += Amp.comp*EVec.Py;                                                            \
      I_States.Pz  += Amp.comp*EVec.Pz;                                                            \
      I_States.Egy += Amp.comp*EVec.Egy;                                                           \
                                                                                                   \
      /* skip the degenerate states */                                                             \
      if ( EVal.comp_next > EVal.comp )                                                            \
      {                                                                                            \
         I_Pres = I_States.Egy - (real)0.5*( I_States.Px*I_States.Px + I_States.Py*I_States.Py +   \
                                             I_States.Pz*I_States.Pz ) / I_States.Rho;             \
                                                                                                   \
         if ( I_States.Rho <= (real)0.0  ||  I_Pres <= (real)0.0 )                                 \
         {                                                                                         \
            Recalculate_Flux( L, R, Flux_Out );                                                    \
                                                                                                   \
            Flux_Out = CUFLU_Rotate3D( Flux_Out, XYZ, false );                                     \
                                                                                                   \
            return Flux_Out;                                                                       \
         }                                                                                         \
      }                                                                                            \
   } // Get_I_States

   I_States.Rho = L.Rho;
   I_States.Px  = L.Px;
   I_States.Py  = L.Py;
   I_States.Pz  = L.Pz;
   I_States.Egy = L.Egy;

   Get_I_States( Rho, Px,  EVec1 );
   Get_I_States( Px,  Py,  EVec2 );
   Get_I_States( Py,  Pz,  EVec3 );
   Get_I_States( Pz,  Egy, EVec4 );

#  undef Recalculate_Flux
#  undef Get_I_States

#  endif // #ifdef CHECK_INTERMEDIATE


// 9. evaluate the ROE fluxes
   Amp.Rho *= FABS( EVal.Rho );
   Amp.Px  *= FABS( EVal.Px  );
   Amp.Py  *= FABS( EVal.Py  );
   Amp.Pz  *= FABS( EVal.Pz  );
   Amp.Egy *= FABS( EVal.Egy );

#  define GetFlux( comp )                                                                    \
   {                                                                                         \
      Flux_Out.comp = (real)0.5*(  Flux_L.comp + Flux_R.comp - ( Amp.Rho*EVec1.comp +        \
                                                                 Amp.Px *EVec2.comp +        \
                                                                 Amp.Py *EVec3.comp +        \
                                                                 Amp.Pz *EVec4.comp +        \
                                                                 Amp.Egy*EVec5.comp )  );    \
   } // GetFlux

   GetFlux( Rho );
   GetFlux( Px  );
   GetFlux( Py  );
   GetFlux( Pz  );
   GetFlux( Egy );

#  undef GetFlux


// 10. evaluate the fluxes for passive scalars
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


// 11. restore the correct order
   Flux_Out = CUFLU_Rotate3D( Flux_Out, XYZ, false );

   return Flux_Out;

} // FUNCTION : CUFLU_RiemannSolver_Roe



#endif // #ifndef __CUFLU_RIEMANNSOLVER_ROE_CU__
