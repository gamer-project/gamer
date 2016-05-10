#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

#if ( MODEL == HYDRO )



extern void CPU_Rotate3D( real InOut[], const int XYZ, const bool Forward );
extern void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma );
#if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
extern real CPU_PositivePres( const real Pres_In, const real Dens, const real _Dens );
extern real CPU_PositivePres_In_Engy( const real ConVar[], const real Gamma_m1, const real _Gamma_m1 );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannSolver_HLLE
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                Estimate the wave speed by Einfeldt et al. (1991).
//
// Note        :  1. The input data should be conserved variables 
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro 
//                         b. Einfeldt, B., et al. J. 1991, J. Comput. Phys., 92, 273
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ      : Targeted spatial direction : (0/1/2) --> (x/y/z)  
//                Flux_Out : Array to store the output flux
//                L_In     : Input left  state (conserved variables)
//                R_In     : Input right state (conserved variables)
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                             const real Gamma )
{

// 1. reorder the input variables for different spatial directions
   real L[5], R[5];

   for (int v=0; v<5; v++)    
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   CPU_Rotate3D( L, XYZ, true );
   CPU_Rotate3D( R, XYZ, true );


// 2. evaluate the Roe's average values
   const real Gamma_m1 = Gamma - (real)1.0;
   real _RhoL, _RhoR, P_L, P_R, H_L, H_R, u, v, w, V2, H, Cs; 
   real RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum, GammaP_Rho;

#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   const real  TempRho  = (real)0.5*( L[0] + R[0] );
   const real _TempRho  = (real)1.0/TempRho;
   real TempPres;
#  endif
   
   _RhoL = (real)1.0 / L[0];
   _RhoR = (real)1.0 / R[0];
   P_L   = Gamma_m1*(  L[4] - (real)0.5*( L[1]*L[1] + L[2]*L[2] + L[3]*L[3] )*_RhoL  );
   P_R   = Gamma_m1*(  R[4] - (real)0.5*( R[1]*R[1] + R[2]*R[2] + R[3]*R[3] )*_RhoR  );
#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   P_L   = CPU_PositivePres( P_L, L[0], _RhoL );
   P_R   = CPU_PositivePres( P_R, R[0], _RhoR );
#  endif
   H_L   = ( L[4] + P_L )*_RhoL;  
   H_R   = ( R[4] + P_R )*_RhoR;  

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(L[0]) )
      Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   L[0], __FILE__, __LINE__, __FUNCTION__ );
   if ( CPU_CheckNegative(R[0]) )
      Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
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
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

   GammaP_Rho = Gamma_m1*( H - (real)0.5*V2 );
#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   TempPres   = GammaP_Rho*TempRho/Gamma;
   TempPres   = CPU_PositivePres( TempPres, TempRho, _TempRho );
   GammaP_Rho = Gamma*TempPres*_TempRho;
#  endif
#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(GammaP_Rho) )
      Aux_Message( stderr, "ERROR : negative GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n", 
                   GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Cs = SQRT( GammaP_Rho );


// 3. estimate the maximum wave speeds
   const real EVal[NCOMP] = { u-Cs, u, u, u, u+Cs };
   real u_L, u_R, Cs_L, Cs_R, MaxV_L, MaxV_R;

   u_L    = _RhoL*L[1];
   u_R    = _RhoR*R[1];

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CPU_CheckNegative(P_L) )
      Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   P_L, __FILE__, __LINE__, __FUNCTION__ );
   if ( CPU_CheckNegative(P_R) )
      Aux_Message( stderr, "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                   P_R, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Cs_L   = SQRT( Gamma*P_L*_RhoL );
   Cs_R   = SQRT( Gamma*P_R*_RhoR );
   MaxV_L = FMIN( EVal[      0], u_L-Cs_L );
   MaxV_R = FMAX( EVal[NCOMP-1], u_R+Cs_R );
   MaxV_L = FMIN( MaxV_L, (real)0.0 );
   MaxV_R = FMAX( MaxV_R, (real)0.0 );


// 4. evaluate the left and right fluxes along the maximum wave speeds
   real Flux_L[5], Flux_R[5];

   CPU_Con2Flux( 0, Flux_L, L, Gamma );
   CPU_Con2Flux( 0, Flux_R, R, Gamma );

   for (int v=0; v<5; v++)
   {
      Flux_L[v] -= MaxV_L*L[v];
      Flux_R[v] -= MaxV_R*R[v];
   }


// 5. evaluate the HLLE fluxes
   const real _MaxV_R_minus_L = (real)1.0 / ( MaxV_R - MaxV_L );

   for (int v=0; v<5; v++)    Flux_Out[v] = _MaxV_R_minus_L*( MaxV_R*Flux_L[v] - MaxV_L*Flux_R[v] );


// 6. restore the correct order
   CPU_Rotate3D( Flux_Out, XYZ, false );

} // FUNCTION : CPU_RiemannSolver_HLLE



#endif // #if ( MODEL == HYDRO )
