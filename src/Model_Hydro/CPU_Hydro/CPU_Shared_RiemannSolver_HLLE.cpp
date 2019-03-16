#ifndef __CUFLU_RIEMANNSOLVER_HLLE__
#define __CUFLU_RIEMANNSOLVER_HLLE__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real Gamma_m1, const real MinPres );
real Hydro_CheckMinPres( const real InPres, const real MinPres );

#endif // #ifdef __CUDACC__ ... else ...



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLE
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                Estimate the wave speed by Einfeldt et al. (1991).
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : (a) Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                         (c) Einfeldt et al., J. Comput. Phys., 92, 273 (1991)
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out : Array to store the output flux
//                L_In     : Input left  state (conserved variables)
//                R_In     : Input right state (conserved variables)
//                Gamma    : Ratio of specific heats
//                MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres )
{

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL_PLUS_MAG], R[NCOMP_TOTAL_PLUS_MAG];

   for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );


// 2. evaluate the Roe's average values
   const real ZERO     = (real)0.0;
   const real ONE      = (real)1.0;
   const real TWO      = (real)2.0;
   const real _TWO     = (real)0.5;
   const real Gamma_m1 = Gamma - ONE;
#  ifdef MHD
   const real Gamma_m2 = Gamma - TWO;
#  endif

   real Rho, _Rho, _RhoL, _RhoR, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real P_L, P_R, H_L, H_R, u, v, w, V2, H, Cf, Cf2, GammaP_Rho;
#  ifdef MHD
   real BxL, ByL, BzL, BxR, ByR, BzR;           // magnetic field from left and right states
   real BxL2, BtL2, BxR2, BtR2, BL2, BR2;       // ...
   real Bx, By, Bz, B2, Bt2, B2_Rho;            // Roe-average magnetic field
   real a2, Cax2, Cat2;                         // wave speeds
   real Ca2_plus_a2, Ca2_min_a2, Cf2_min_Cs2;   // Ca^2+a^2, Ca^2-a^2, Cf^2-Cs^2
   real X, Y;                                   // Eqs. (B15) and (B16) in ref-b
#  endif

   _RhoL = ONE / L[0];
   _RhoR = ONE / R[0];
   P_L   = Gamma_m1*(  L[4] - _TWO*( SQR(L[1]) + SQR(L[2]) + SQR(L[3]) )*_RhoL  );
   P_R   = Gamma_m1*(  R[4] - _TWO*( SQR(R[1]) + SQR(R[2]) + SQR(R[3]) )*_RhoR  );
#  ifdef MHD
   BxL   = L[ MAG_OFFSET + 0 ];
   ByL   = L[ MAG_OFFSET + 1 ];
   BzL   = L[ MAG_OFFSET + 2 ];
   BxR   = R[ MAG_OFFSET + 0 ];
   ByR   = R[ MAG_OFFSET + 1 ];
   BzR   = R[ MAG_OFFSET + 2 ];
   BxL2  = SQR( BxL );
   BxR2  = SQR( BxR );
   BtL2  = SQR( ByL ) + SQR( BzL );
   BtR2  = SQR( ByR ) + SQR( BzR );
   BL2   = BxL2 + BtL2;
   BR2   = BxR2 + BtR2;
   P_L  -= _TWO*Gamma_m1*BL2;
   P_R  -= _TWO*Gamma_m1*BR2;
#  endif // #ifdef MHD
   P_L   = Hydro_CheckMinPres( P_L, MinPres );
   P_R   = Hydro_CheckMinPres( P_R, MinPres );
   H_L   = ( L[4] + P_L )*_RhoL;
   H_R   = ( R[4] + P_R )*_RhoR;
#  ifdef MHD
   H_L  += _TWO*BL2*_RhoL;
   H_R  += _TWO*BR2*_RhoR;
#  endif

#  ifdef GAMER_DEBUG
// longitudinal B field in the left and right states should be the same
#  ifdef MHD
   if ( BxL != BxR )
      printf( "ERROR : BxL (%24.17e) != BxR (%24.17e) for XYZ %d at file <%s>, line <%d>, function <%s>!!\n",
              BxL, BxR, XYZ, __FILE__, __LINE__, __FUNCTION__ );
#  endif

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(L[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              L[0], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(R[0]) )
      printf( "ERROR : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              R[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif
#  endif // #ifdef GAMER_DEBUG

   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   Rho             = RhoL_sqrt*RhoR_sqrt;
   _Rho            = ONE/Rho;
   _RhoL_sqrt      = ONE/RhoL_sqrt;
   _RhoR_sqrt      = ONE/RhoR_sqrt;
   _RhoLR_sqrt_sum = ONE/(RhoL_sqrt + RhoR_sqrt);

   u  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w  = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   V2 = u*u + v*v + w*w;
   H  = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

#  ifdef MHD
   Bx     = BxL;
   By     = _RhoLR_sqrt_sum*( RhoL_sqrt*ByR + RhoR_sqrt*ByL );
   Bz     = _RhoLR_sqrt_sum*( RhoL_sqrt*BzR + RhoR_sqrt*BzL );
   Bt2    = SQR( By ) + SQR( Bz );
   B2     = SQR( Bx ) + Bt2;
   B2_Rho = B2*_Rho;
   X      = _TWO*( SQR(ByR-ByL) + SQR(BzR-BzL) )*SQR( _RhoLR_sqrt_sum );
   X     *= Gamma_m2;
#  ifdef EULERY
   Y      = _TWO*( L[0] + R[0] )*_Rho ;
#  else
   Y      = ONE;
#  endif
   Y     *= Gamma_m2;
#  endif // #ifdef MHD

   GammaP_Rho  = Gamma_m1*( H - _TWO*V2 );
#  ifdef MHD
   GammaP_Rho -= Gamma_m1*B2_Rho;   // H = 0.5*v^2 + B^2/rho + gamma/(gamma-1)*P/rho
#  endif
   GammaP_Rho  = Gamma*_Rho*Hydro_CheckMinPres( GammaP_Rho*Rho/Gamma, MinPres );  // apply pressure floor

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(GammaP_Rho) )
      printf( "ERROR : invalid GammaP_Rho (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              GammaP_Rho, __FILE__, __LINE__, __FUNCTION__ );
#  ifdef MHD
   if ( Hydro_CheckNegative(Gamma_m1-Y) )
      printf( "ERROR : invalid Gamma_m1-Y (%14.7e, Gamma_m1 %14.7e, Y %14.7e) at file <%s>, line <%d>, function <%s>\n",
              Gamma_m1-Y, Gamma_m1, Y, __FILE__, __LINE__, __FUNCTION__ );
#  endif
#  endif // #ifdef CHECK_NEGATIVE_IN_FLUID

#  ifdef MHD
   a2          = GammaP_Rho - X;
#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(a2) )
      printf( "ERROR : invalid a2 (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              a2, __FILE__, __LINE__, __FUNCTION__ );
#  endif
   Cax2        = SQR(Bx)*_Rho;
   Cat2        = ( Gamma_m1 - Y )*Bt2*_Rho;
   Ca2_plus_a2 = Cat2 + Cax2 + a2;
   Ca2_min_a2  = Cat2 + Cax2 - a2;
   Cf2_min_Cs2 = SQRT( SQR(Ca2_min_a2) + (real)4.0*a2*Cat2 );

// evaluate the fast wave speed Cf for the Roe's average values
   if ( Cat2 == ZERO )
   {
      if      ( Cax2 == a2 )  Cf2 = a2;
      else if ( Cax2 >  a2 )  Cf2 = Cax2;
      else                    Cf2 = a2;
   }

   else
   {
      if ( Cax2 == ZERO )     Cf2 = a2 + Cat2;
      else                    Cf2 = _TWO*( Ca2_plus_a2 + Cf2_min_Cs2 );
   } // if ( Cat2 == ZERO ) ... else ...

#  else  // #ifdef MHD
   Cf2 = GammaP_Rho;
#  endif // #ifdef MHD ... else ...
   Cf  = SQRT( Cf2 );


// 3. estimate the maximum wave speeds
   const real EVal_min = u - Cf;
   const real EVal_max = u + Cf;
   real u_L, u_R, Cf_L, Cf_R, MaxV_L, MaxV_R;

   u_L = _RhoL*L[1];
   u_R = _RhoR*R[1];

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(P_L) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_L, __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(P_R) )
      printf( "ERROR : invalid pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              P_R, __FILE__, __LINE__, __FUNCTION__ );
#  endif

#  ifdef MHD
// evaluate the fast wave speed Cf for the left state
   a2          = Gamma*P_L*_RhoL;
   Cax2        = BxL2*_RhoL;
   Cat2        = BtL2*_RhoL;
   Ca2_plus_a2 = Cat2 + Cax2 + a2;
   Ca2_min_a2  = Cat2 + Cax2 - a2;
   Cf2_min_Cs2 = SQRT( SQR(Ca2_min_a2) + (real)4.0*a2*Cat2 );

   if ( Cat2 == ZERO )
   {
      if      ( Cax2 == a2 )  Cf2 = a2;
      else if ( Cax2 >  a2 )  Cf2 = Cax2;
      else                    Cf2 = a2;
   }

   else
   {
      if ( Cax2 == ZERO )     Cf2 = a2 + Cat2;
      else                    Cf2 = _TWO*( Ca2_plus_a2 + Cf2_min_Cs2 );
   } // if ( Cat2 == ZERO ) ... else ...

   Cf_L = SQRT( Cf2 );

// evaluate the fast wave speed Cf for the right state
   a2          = Gamma*P_R*_RhoR;
   Cax2        = BxR2*_RhoR;
   Cat2        = BtR2*_RhoR;
   Ca2_plus_a2 = Cat2 + Cax2 + a2;
   Ca2_min_a2  = Cat2 + Cax2 - a2;
   Cf2_min_Cs2 = SQRT( SQR(Ca2_min_a2) + (real)4.0*a2*Cat2 );

   if ( Cat2 == ZERO )
   {
      if      ( Cax2 == a2 )  Cf2 = a2;
      else if ( Cax2 >  a2 )  Cf2 = Cax2;
      else                    Cf2 = a2;
   }

   else
   {
      if ( Cax2 == ZERO )     Cf2 = a2 + Cat2;
      else                    Cf2 = _TWO*( Ca2_plus_a2 + Cf2_min_Cs2 );
   } // if ( Cat2 == ZERO ) ... else ...

   Cf_R = SQRT( Cf2 );

#  else // #ifdef MHD
   Cf_L   = SQRT( Gamma*P_L*_RhoL );
   Cf_R   = SQRT( Gamma*P_R*_RhoR );
#  endif // #ifdef MHD ... else ...

   MaxV_L = FMIN( EVal_min, u_L-Cf_L );
   MaxV_R = FMAX( EVal_max, u_R+Cf_R );
   MaxV_L = FMIN( MaxV_L, ZERO );
   MaxV_R = FMAX( MaxV_R, ZERO );


// 4. evaluate the left and right fluxes along the maximum wave speeds
#  ifdef MHD
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4, MAG_OFFSET+1, MAG_OFFSET+2 };
#  else
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };
#  endif
   real Flux_L[NCOMP_TOTAL_PLUS_MAG], Flux_R[NCOMP_TOTAL_PLUS_MAG];  // use NCOMP_TOTAL_PLUS_MAG for Hydro_Con2Flux()

   Hydro_Con2Flux( 0, Flux_L, L, Gamma_m1, MinPres );
   Hydro_Con2Flux( 0, Flux_R, R, Gamma_m1, MinPres );

   for (int v=0; v<NWAVE; v++)
   {
      Flux_L[ idx_wave[v] ] -= MaxV_L*L[ idx_wave[v] ];
      Flux_R[ idx_wave[v] ] -= MaxV_R*R[ idx_wave[v] ];
   }


// 5. evaluate the HLLE fluxes
   const real _MaxV_R_minus_L = ONE / ( MaxV_R - MaxV_L );

   for (int v=0; v<NWAVE; v++)
      Flux_Out[ idx_wave[v] ] = _MaxV_R_minus_L*( MaxV_R*Flux_L[ idx_wave[v] ] - MaxV_L*Flux_R[ idx_wave[v] ] );

// longitudinal magnetic flux is always zero
#  ifdef MHD
   Flux_Out[MAG_OFFSET] = ZERO;
#  endif


// 6. evaluate the fluxes for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out[FLUX_DENS] >= ZERO )
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


// 7. restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );

} // FUNCTION : Hydro_RiemannSolver_HLLE



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLE__
