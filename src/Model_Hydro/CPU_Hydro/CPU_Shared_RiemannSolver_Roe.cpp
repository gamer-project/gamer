#ifndef __CUFLU_RIEMANNSOLVER_ROE__
#define __CUFLU_RIEMANNSOLVER_ROE__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#if   ( CHECK_INTERMEDIATE == EXACT )
# include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( CHECK_INTERMEDIATE == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( CHECK_INTERMEDIATE == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres );
#if   ( CHECK_INTERMEDIATE == EXACT )
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_RiemannSolver_Exact( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( CHECK_INTERMEDIATE == HLLE )
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres );
#elif ( CHECK_INTERMEDIATE == HLLC )
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres );
#endif
real Hydro_CheckMinPres( const real InPres, const real MinPres );

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_Roe
// Description :  Approximate Riemann solver of Roe
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : "Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                         ~ by Eleuterio F. Toro"
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ         : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out    : Array to store the output flux
//                L_In        : Input left  state (conserved variables)
//                R_In        : Input right state (conserved variables)
//                Gamma       : Ratio of specific heats
//                MinPres     : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                              const real Gamma, const real MinPres )
{

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL], R[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true );
   Hydro_Rotate3D( R, XYZ, true );


// 2. evaluate the average values
   const real Gamma_m1 = Gamma - (real)1.0;
   const real  TempRho = (real)0.5*( L[0] + R[0] );
   const real _TempRho = (real)1.0/TempRho;
   real _RhoL, _RhoR, HL, HR, u, v, w, V2, H, Cs, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real GammaP_Rho, TempPres;

   _RhoL = (real)1.0 / L[0];
   _RhoR = (real)1.0 / R[0];
   HL    = (  L[4] + Gamma_m1*( L[4] - (real)0.5*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] )*_RhoL )  )*_RhoL;
   HR    = (  R[4] + Gamma_m1*( R[4] - (real)0.5*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] )*_RhoR )  )*_RhoR;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(L[0]) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[0], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(R[0]) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   _RhoL_sqrt      = (real)1.0 / RhoL_sqrt;
   _RhoR_sqrt      = (real)1.0 / RhoR_sqrt;
   _RhoLR_sqrt_sum = (real)1.0 / (RhoL_sqrt + RhoR_sqrt);

   u  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   V2 = u*u + v*v + w*w;
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*HL   +  RhoR_sqrt*HR   );

   GammaP_Rho = Gamma_m1*( H - (real)0.5*V2 );
   TempPres   = GammaP_Rho*TempRho/Gamma;
   TempPres   = Hydro_CheckMinPres( TempPres, MinPres );
   GammaP_Rho = Gamma*TempPres*_TempRho;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(GammaP_Rho) )
      printf( "ERROR : negative GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Cs = SQRT( GammaP_Rho );


// 3. evaluate the eigenvalues
   const real EigenVal[NCOMP_FLUID] = { u-Cs, u, u, u, u+Cs };


// 4. evaluate the left and right fluxes
   real Flux_L[NCOMP_TOTAL], Flux_R[NCOMP_TOTAL];

   Hydro_Con2Flux( 0, Flux_L, L, Gamma_m1, MinPres );
   Hydro_Con2Flux( 0, Flux_R, R, Gamma_m1, MinPres );


// 5. return the upwind fluxes if flow is supersonic
   if ( EigenVal[0] >= (real)0.0 )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Out[v] = Flux_L[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false );

      return;
   }

   if ( EigenVal[4] <= (real)0.0 )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Out[v] = Flux_R[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false );

      return;
   }


// 6. evaluate the amplitudes along different characteristics (eigenvectors)
   real Jump[NCOMP_FLUID], Amp[NCOMP_FLUID];

   for (int v=0; v<NCOMP_FLUID; v++)   Jump[v] = R[v] - L[v];

   Amp[2] = Jump[2] - v*Jump[0];
   Amp[3] = Jump[3] - w*Jump[0];
   Amp[1] = Gamma_m1/(Cs*Cs)*( Jump[0]*(H-u*u) + u*Jump[1] - Jump[4] + v*Amp[2] + w*Amp[3] );
   Amp[0] = (real)0.5/Cs*( Jump[0]*(u+Cs) - Jump[1] - Cs*Amp[1] );
   Amp[4] = Jump[0] - Amp[0] - Amp[1];


// 7. evaluate the eigenvectors
   const real EigenVec[NCOMP_FLUID][NCOMP_FLUID] = {  { (real)1.0,       u-Cs,         v,         w,       H-u*Cs },
                                                      { (real)1.0,          u,         v,         w, (real)0.5*V2 },
                                                      { (real)0.0,  (real)0.0, (real)1.0, (real)0.0,            v },
                                                      { (real)0.0,  (real)0.0, (real)0.0, (real)1.0,            w },
                                                      { (real)1.0,       u+Cs,         v,         w,       H+u*Cs }  };


// 8. verify that the density and pressure in the intermediate states are positive
#  ifdef CHECK_INTERMEDIATE
   real I_Pres, I_States[NCOMP_FLUID];

#  if ( CHECK_INTERMEDIATE == EXACT )
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif


   for (int v=0; v<NCOMP_FLUID; v++)   I_States[v] = L[v];

   for (int t=0; t<4; t++)
   {
      for (int v=0; v<NCOMP_FLUID; v++)   I_States[v] += Amp[t]*EigenVec[t][v];

      if ( EigenVal[t+1] > EigenVal[t] )  // skip the degenerate states
      {
         I_Pres = I_States[4] - (real)0.5*( I_States[1]*I_States[1] + I_States[2]*I_States[2] +
                                            I_States[3]*I_States[3] ) / I_States[0];

         if ( I_States[0] <= (real)0.0  ||  I_Pres <= (real)0.0 )
         {
#           ifdef GAMER_DEBUG
            printf( "WARNING : intermediate states check failed !!\n" );
#           endif

#           if   ( CHECK_INTERMEDIATE == EXACT )   // recalculate fluxes by exact solver

            const bool NormPassive_No  = false; // do NOT convert any passive variable to mass fraction for the Riemann solvers
            const bool JeansMinPres_No = false;

            Hydro_Con2Pri( L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
            Hydro_Con2Pri( R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

            Hydro_RiemannSolver_Exact( 0, Flux_Out, PriVar_L, PriVar_R, Gamma );

#           elif ( CHECK_INTERMEDIATE == HLLE )    // recalculate fluxes by HLLE solver

            Hydro_RiemannSolver_HLLE ( 0, Flux_Out, L, R, Gamma, MinPres );

#           elif ( CHECK_INTERMEDIATE == HLLC )    // recalculate fluxes by HLLC solver

            Hydro_RiemannSolver_HLLC ( 0, Flux_Out, L, R, Gamma, MinPres );

#           else

#           error : ERROR : unsupported CHECK_INTERMEDIATE (EXACT/HLLE/HLLC) !!

#           endif

            Hydro_Rotate3D( Flux_Out, XYZ, false );
            return;

         } // if ( I_States[0] <= (real)0.0  ||  I_Pres <= (real)0.0 )
      } // if ( EigenVal[t+1] > EigenVal[t] )
   } // for (int t=0; t<4; t++)
#  endif // #ifdef CHECK_INTERMEDIATE


// 9. evaluate the Roe fluxes
   for (int v=0; v<NCOMP_FLUID; v++)   Amp[v] *= FABS( EigenVal[v] );

   for (int v=0; v<NCOMP_FLUID; v++)
      Flux_Out[v] = (real)0.5*( Flux_L[v] + Flux_R[v] ) - (real)0.5*(   Amp[0]*EigenVec[0][v]
                                                                      + Amp[1]*EigenVec[1][v]
                                                                      + Amp[2]*EigenVec[2][v]
                                                                      + Amp[3]*EigenVec[3][v]
                                                                      + Amp[4]*EigenVec[4][v] );

// 10. evaluate the fluxes for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoL;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
   }

   else
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoR;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
   }
#  endif


// 11. restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false );

} // FUNCTION : Hydro_RiemannSolver_Roe



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_ROE__
