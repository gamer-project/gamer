#ifndef __CUFLU_RIEMANNSOLVER_HLLC__
#define __CUFLU_RIEMANNSOLVER_HLLC__



#include "CUFLU.h"

#if ( MODEL == HYDRO )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real MinPres, const long PassiveFloor,
                     const EoS_DE2P_t EoS_DensEint2Pres, const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                     const real *const EoS_Table[EOS_NTABLE_MAX], const real* const PresIn );

#endif // #ifdef __CUDACC__ ... else ...




// ##################
// ## SRHD version ##
// ##################
#ifdef SRHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer extended to include the contact wave
//                (SRHD version)
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro (1999)
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput., 18, 1553
//                         c. Coleman, M. S. B. 2020, ApJS, 248, 7
//                3. Wave-speed estimator is set by HLLC_WAVESPEED in CUFLU.h
//                4. Support general EoS
//                5. Shared by the MHM, MHM_RP, and CTU schemes
//
// Parameter   :  XYZ               : Target spatial direction : (0/1/2) --> (x/y/z)
//                Flux_Out          : Array to store the output flux
//                L/R_In            : Input left/right states (conserved variables)
//                MinDens/Pres      : Density and pressure floors
//                PassiveFloor      : Bitwise flag to specify the passive scalars to be floored
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed squared
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//                FreezeHydro       : Freeze hydrodynamic fluxes
//
// Return      :  Flux_Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const long PassiveFloor, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                               const real* const EoS_Table[EOS_NTABLE_MAX], bool FreezeHydro )
{

   if ( FreezeHydro )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Out[v] = 0.0;
      return;
   }

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL], R[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical( UNPHY_MODE_CONS, L, NULL_REAL,
                       EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table,
                       PassiveFloor, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical( UNPHY_MODE_CONS, R, NULL_REAL,
                       EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table,
                       PassiveFloor, ERROR_INFO, UNPHY_VERBOSE );
#  endif


// 2. compute primitive variables from conserved variables
   real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL];
   real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
   real Usl[NCOMP_TOTAL], Usr[NCOMP_TOTAL];
   real cslsq, csrsq, gammasql, gammasqr;
   real ssl, ssr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
   real lmdal,lmdar;
   real lmdas;
   real a,b,c;
   real den,ps;
   real lV1, rV1, lV2, rV2, lV3, rV3;
   real lFactor,rFactor;

   Hydro_Con2Pri( L, PL, MinPres, PassiveFloor, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                  (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, &lFactor );

   Hydro_Con2Pri( R, PR, MinPres, PassiveFloor, NULL_BOOL, NULL_INT, NULL, NULL_BOOL,
                  (real)NULL_REAL, NULL, NULL, EoS_GuessHTilde, EoS_HTilde2Temp,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL, &rFactor );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical( UNPHY_MODE_PRIM, PL, NULL_REAL,
                       EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table,
                       PassiveFloor, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical( UNPHY_MODE_PRIM, PR, NULL_REAL,
                       EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                       EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table,
                       PassiveFloor, ERROR_INFO, UNPHY_VERBOSE );
#  endif


// 3. transform 4-velocity to 3-velocity
   const real _lFactor = (real)1.0 / lFactor;
   lV1 = PL[1]*_lFactor;
   lV2 = PL[2]*_lFactor;
   lV3 = PL[3]*_lFactor;

   const real _rFactor = (real)1.0 / rFactor;
   rV1 = PR[1]*_rFactor;
   rV2 = PR[2]*_rFactor;
   rV3 = PR[3]*_rFactor;


// 4. compute the max and min wave speeds used in Mignone
   cslsq = EoS_DensPres2CSqr( PL[0], PL[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );
   csrsq = EoS_DensPres2CSqr( PR[0], PR[4], NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   if ( cslsq >= (real)1.0  ||  csrsq >= (real)1.0  ||  cslsq < (real)0.0  ||  csrsq < (real)0.0 )
      printf( "cslsq = %14.7e, cslrq = %14.7e\n", cslsq, csrsq );
#  endif

// square of Lorentz factor
   gammasql = SQR( lFactor );
   gammasqr = SQR( rFactor );
   ssl      = cslsq / ( -gammasql*cslsq + gammasql ); // Mignone Eq 22.5
   ssr      = csrsq / ( -gammasqr*csrsq + gammasqr ); // Mignone Eq 22.5

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   if ( ssl < (real)0.0  ||  ssr < (real)0.0 )  printf( "ssl = %14.7e, ssr = %14.7e\n", ssl, ssr );
#  endif

   const real lV2s       = lV2*lV2;
   const real rV2s       = rV2*rV2;
   const real lV3s       = lV3*lV3;
   const real rV3s       = rV3*rV3;
   const real __gammasql = (real)1.0 / gammasql;
   const real __gammasqr = (real)1.0 / gammasqr;
   const real deltal     = ssl*ssl + ssl*( __gammasql + lV2s + lV3s );
   const real deltar     = ssr*ssr + ssr*( __gammasqr + rV2s + rV3s );
   const real ssl__      = (real)1.0 + ssl;
   const real ssr__      = (real)1.0 + ssr;

   lmdapl = ( lV1 + SQRT(deltal) ) / ssl__;
   lmdaml = ( lV1 - SQRT(deltal) ) / ssl__;
   lmdapr = ( rV1 + SQRT(deltar) ) / ssr__;
   lmdamr = ( rV1 - SQRT(deltar) ) / ssr__;
   lmdal  = FMIN( lmdaml, lmdamr ); // Mignone Eq 21
   lmdar  = FMAX( lmdapl, lmdapr );


// 5. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18))
//    compute HLL conserved quantities using Mignone Eq 9
   Fl[0] = L[0]*lV1;
   Fl[1] = L[1]*lV1 + PL[4];
   Fl[2] = L[2]*lV1;
   Fl[3] = L[3]*lV1;
   Fl[4] = ( L[4] + PL[4] )*lV1;

   if ( lmdal >= (real)0.0 )
   {
//    Fl
//    intercell flux is left flux
      Flux_Out[0] = Fl[0];
      Flux_Out[1] = Fl[1];
      Flux_Out[2] = Fl[2];
      Flux_Out[3] = Fl[3];
      Flux_Out[4] = Fl[4];

//    evaluate the fluxes of passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
      {
         const real vx = Flux_Out[FLUX_DENS]/L[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
      }

      else
      {
         const real vx = Flux_Out[FLUX_DENS]/R[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
      }
#     endif

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );
      return;
   } // if ( lmdal >= (real)0.0 )

   Fr[0] = R[0]*rV1;
   Fr[1] = R[1]*rV1 + PR[4];
   Fr[2] = R[2]*rV1;
   Fr[3] = R[3]*rV1;
   Fr[4] = ( R[4] + PR[4] )*rV1;

   if ( lmdar <= (real)0.0 )
   {
//    Fr
//    intercell flux is right flux
      Flux_Out[0] = Fr[0];
      Flux_Out[1] = Fr[1];
      Flux_Out[2] = Fr[2];
      Flux_Out[3] = Fr[3];
      Flux_Out[4] = Fr[4];

//    evaluate the fluxes of passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
      {
         const real vx = Flux_Out[FLUX_DENS]/L[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
      }

      else
      {
         const real vx = Flux_Out[FLUX_DENS]/R[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
      }
#     endif

     Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );
     return;
   } // if ( lmdar <= (real)0.0 )


// 6. compute contact wave speed using larger root from Mignone Eq 18
//    --> physical root is the root with the minus sign
   lmdatlmda = lmdal*lmdar;

// quadratic formuLa calcuLation
   a = lmdar*L[1]
     - lmdal*R[1]
     + lmdatlmda*( R[4] + R[0] - L[4] - L[0] );

   b = lmdal*( L[4] + L[0] ) - L[1]
     - lmdar*( R[4] + R[0] ) + R[1]
     + lmdal*( R[1]*rV1 + PR[4] )
     - lmdar*( L[1]*lV1 + PL[4] )
     - lmdatlmda*( R[1] - L[1] );

   c = lmdar*R[1] - lmdal*L[1] - ( R[1]*rV1 + PR[4] ) + ( L[1]*lV1 + PL[4] );

   const real delta = b*b - (real)4.0*a*c;

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   if ( delta < (real)0.0 )   printf("delta=%14.7e\n", delta );
#  endif

   lmdas = - ( (real)2.0*c ) / ( b + SIGN(b)*SQRT(delta) );

   ps = lmdas*( ( R[4] + R[0] )*( rV1 - lmdar ) + PR[4]*rV1 ) - R[1]*( rV1 - lmdar ) - PR[4];
   ps /= ( lmdas*lmdar - (real)1.0 );

// ps = lmdas*( ( L[4] + L[0] )*( lV1 - lmdal ) + PL[4]*lV1 ) - L[1]*( lV1 - lmdal ) - PL[4];
// ps /= ( lmdas*lmdal - (real)1.0 );


// 7. determine intercell flux according to Mignone 13
   if ( lmdas >= (real)0.0 )
   {
//    Fls
//    now calculate Usl with Mignone Eq 16
      den = (real)1.0 / ( lmdal - lmdas );

      real factor0 = lmdal - lV1;
      real factor1 = lmdal*den - lV1*den;

      Usl[0] = L[0]*factor1;
      Usl[1] = ( L[1]*factor0 + ps - PL[4] )*den;
      Usl[2] = L[2]*factor1;
      Usl[3] = L[3]*factor1;
      Usl[4] = (  -PL[4]*lV1 + ( L[4]*factor0 + ps*lmdas )  )*den;

//    now calculate Fsr using Mignone Eq 14
      Flux_Out[0] = lmdal*( Usl[0] - L[0] ) + Fl[0];
      Flux_Out[1] = lmdal*( Usl[1] - L[1] ) + Fl[1];
      Flux_Out[2] = lmdal*( Usl[2] - L[2] ) + Fl[2];
      Flux_Out[3] = lmdal*( Usl[3] - L[3] ) + Fl[3];
      Flux_Out[4] = lmdal*( Usl[4] - L[4] ) + Fl[4];

//    evaluate the fluxes of passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
      {
         const real vx = Flux_Out[FLUX_DENS]/L[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
      }

      else
      {
         const real vx = Flux_Out[FLUX_DENS]/R[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
      }
#     endif

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );
      return;
   } // if ( lmdas >= (real)0.0 )

   else
   {
//    Frs
//    now calculate Usr with Mignone Eq 16
      den = (real)1.0 / ( lmdar - lmdas );

    const real factor0 = lmdar - rV1;
    const real factor1 = lmdar * den - rV1*den;

    Usr[0] = R[0]*factor1;
    Usr[1] = (  R[1]*factor0 + ( ps - PR[4] )  )*den;
    Usr[2] = R[2]*factor1;
    Usr[3] = R[3]*factor1;
    Usr[4] = (  -PR[4]*rV1 + ( R[4]*factor0 + ps*lmdas )  )*den;

//    now calculate Fsr using Mignone Eq 14
      Flux_Out[0] = lmdar*( Usr[0] - R[0] ) + Fr[0];
      Flux_Out[1] = lmdar*( Usr[1] - R[1] ) + Fr[1];
      Flux_Out[2] = lmdar*( Usr[2] - R[2] ) + Fr[2];
      Flux_Out[3] = lmdar*( Usr[3] - R[3] ) + Fr[3];
      Flux_Out[4] = lmdar*( Usr[4] - R[4] ) + Fr[4];

//    evaluate the fluxes of passive scalars
#     if ( NCOMP_PASSIVE > 0 )
      if ( Flux_Out[FLUX_DENS] >= (real)0.0 )
      {
         const real vx = Flux_Out[FLUX_DENS]/L[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = L[v]*vx;
      }

      else
      {
         const real vx = Flux_Out[FLUX_DENS]/R[0];

         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = R[v]*vx;
      }
#     endif

      Hydro_Rotate3D( Flux_Out, XYZ, false, MAG_OFFSET );
      return;
   } // if ( lmdas >= (real)0.0 ) ... else ...

} // FUNCTION : Hydro_RiemannSolver_HLLC
#endif // #ifdef SRHD



// ################
// ## HD version ##
// ################
#ifndef SRHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer extended to include the contact wave
//                (HD version)
//
// Note        :  See the SRHD routine
//
// Parameter   :  See the SRHD routine
//
// Return      :  See the SRHD routine
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const long PassiveFloor, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                               const real* const EoS_Table[EOS_NTABLE_MAX], const bool FreezeHydro )
{

   if ( FreezeHydro )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Out[v] = 0.0;
      return;
   }

// 1. reorder the input variables for different spatial directions
   real L[NCOMP_TOTAL], R[NCOMP_TOTAL];

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      L[v] = L_In[v];
      R[v] = R_In[v];
   }

   Hydro_Rotate3D( L, XYZ, true, MAG_OFFSET );
   Hydro_Rotate3D( R, XYZ, true, MAG_OFFSET );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical_Single( L[0], "density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( R[0], "density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif


// 2. estimate the maximum wave speeds
// 2-1. compute the left/right states
   const real ZERO             = (real)0.0;
   const real ONE              = (real)1.0;
#  if ( HLLC_WAVESPEED != HLL_WAVESPEED_DAVIS )
   const real _TWO             = (real)0.5;
#  endif
   const bool CheckMinPres_Yes = true;
   const real Emag             = NULL_REAL;

   real _RhoL, _RhoR, u_L, u_R, P_L, P_R, Cs_L, Cs_R, W_L, W_R;

   _RhoL = ONE / L[0];
   _RhoR = ONE / R[0];
   u_L   = _RhoL*L[1];
   u_R   = _RhoR*R[1];
   P_L   = Hydro_Con2Pres( L[0], L[1], L[2], L[3], L[4], L+NCOMP_FLUID, CheckMinPres_Yes, MinPres, PassiveFloor, Emag,
                           EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   P_R   = Hydro_Con2Pres( R[0], R[1], R[2], R[3], R[4], R+NCOMP_FLUID, CheckMinPres_Yes, MinPres, PassiveFloor, Emag,
                           EoS_DensEint2Pres, EoS_GuessHTilde, EoS_HTilde2Temp,
                           EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table, NULL );
   Cs_L  = SQRT(  EoS_DensPres2CSqr( L[0], P_L, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )  );
   Cs_R  = SQRT(  EoS_DensPres2CSqr( R[0], P_R, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )  );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical_Single( P_R, "pressure", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( P_L, "pressure", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif


// 2-2a. use the Roe average eigenvalues
#  if   ( HLLC_WAVESPEED == HLL_WAVESPEED_ROE )
   real H_L, H_R, RhoL_sqrt, RhoR_sqrt, _RhoL_sqrt, _RhoR_sqrt, _RhoLR_sqrt_sum;
   real u_Roe, v_Roe, w_Roe, H_Roe;

// Roe averages
   H_L             = ( L[4] + P_L )*_RhoL;
   H_R             = ( R[4] + P_R )*_RhoR;
   RhoL_sqrt       = SQRT( L[0] );
   RhoR_sqrt       = SQRT( R[0] );
   _RhoL_sqrt      = ONE / RhoL_sqrt;
   _RhoR_sqrt      = ONE / RhoR_sqrt;
   _RhoLR_sqrt_sum = ONE / ( RhoL_sqrt + RhoR_sqrt );

   u_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[1] + _RhoR_sqrt*R[1] );
   v_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[2] + _RhoR_sqrt*R[2] );
   w_Roe = _RhoLR_sqrt_sum*( _RhoL_sqrt*L[3] + _RhoR_sqrt*R[3] );
   H_Roe = _RhoLR_sqrt_sum*(  RhoL_sqrt*H_L  +  RhoR_sqrt*H_R  );

// sound speed
//###NOTE: we have assumed a constant-gamma EoS here
// --> otherwise, one needs to specify how to convert (H-0.5*V2, Rho) to Cs^2
// --> see Eq. [A4] in Coleman 2020
#  if ( EOS != EOS_GAMMA )
#     error : ERROR : HLL_WAVESPEED_ROE only works with EOS_GAMMA !!
#  endif

   const real Gamma    = (real)EoS_AuxArray_Flt[0];
   const real Gamma_m1 = (real)EoS_AuxArray_Flt[1];
   const real _Gamma   = (real)EoS_AuxArray_Flt[3];

   real V2_Roe, Cs2_Roe, Cs_Roe, TempRho, TempPres;

   V2_Roe   = SQR( u_Roe ) + SQR( v_Roe ) + SQR( w_Roe );
   Cs2_Roe  = Gamma_m1*( H_Roe - _TWO*V2_Roe );
   TempRho  = _TWO*( L[0] + R[0] );  // different from Roe aveage density and is only for applying pressure floor
   TempPres = Cs2_Roe*TempRho*_Gamma;
   TempPres = Hydro_CheckMinPres( TempPres, MinPres );
   Cs2_Roe  = Gamma*TempPres/TempRho;
   Cs_Roe   = SQRT( Cs2_Roe );

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical_Single( Cs2_Roe, "Cs2_Roe", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif

// maximum and minimum eigenvalues
   const real EVal_min = u_Roe - Cs_Roe;
   const real EVal_max = u_Roe + Cs_Roe;

// left/right maximum wave speeds
   W_L = FMIN( EVal_min, u_L-Cs_L );
   W_R = FMAX( EVal_max, u_R+Cs_R );


// 2-2b. use the primitive variable Riemann solver (PVRS)
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_PVRS )
   real Rho_PVRS, Cs_PVRS, RhoCs_PVRS, P_PVRS, Gamma_SL, Gamma_SR, q_L, q_R;

   Rho_PVRS    = _TWO*( L[0] + R[0] );
   Cs_PVRS     = _TWO*( Cs_L + Cs_R );
   RhoCs_PVRS  = Rho_PVRS * Cs_PVRS;
   P_PVRS      = _TWO*(  ( P_L + P_R ) + ( u_L - u_R )*RhoCs_PVRS  );
   P_PVRS      = Hydro_CheckMinPres( P_PVRS, MinPres );

// for EOS_GAMMA/EOS_ISOTHERMAL, the calculations of Gamma_SL/R can be greatly simplified
// --> results should be exactly the same except for round-off errors
#  if   ( EOS == EOS_GAMMA )
   Gamma_SL    = (real)EoS_AuxArray_Flt[0];
   Gamma_SR    = (real)EoS_AuxArray_Flt[0];
#  elif ( EOS == EOS_ISOTHERMAL )
   Gamma_SL    = ONE;
   Gamma_SR    = ONE;
#  else
   real u_PVRS, Rho_Cs_PVRS, Rho_SL, Rho_SR, _P;

   u_PVRS      = _TWO*(  ( u_L + u_R ) + ( P_L - P_R )/RhoCs_PVRS  );
   Rho_Cs_PVRS = Rho_PVRS / Cs_PVRS;
   Rho_SL      = L[0] + ( u_L - u_PVRS )*Rho_Cs_PVRS;
   Rho_SR      = R[0] + ( u_PVRS - u_R )*Rho_Cs_PVRS;
   Rho_SL      = FMAX( Rho_SL, MinDens );
   Rho_SR      = FMAX( Rho_SR, MinDens );
   _P          = ONE / P_PVRS;
// see Eq. [9.8] in Toro 1999 for passive scalars
   Gamma_SL    = EoS_DensPres2CSqr( Rho_SL, P_PVRS, L+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )*Rho_SL*_P;
   Gamma_SR    = EoS_DensPres2CSqr( Rho_SR, P_PVRS, R+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, EoS_Table )*Rho_SR*_P;
#  endif // EOS

   q_L = ( P_PVRS <= P_L ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SL + ONE )/Gamma_SL*( P_PVRS/P_L - ONE )  );
   q_R = ( P_PVRS <= P_R ) ? ONE : SQRT(  ONE + _TWO*( Gamma_SR + ONE )/Gamma_SR*( P_PVRS/P_R - ONE )  );
   W_L = u_L - Cs_L*q_L;
   W_R = u_R + Cs_R*q_R;

#  ifdef CHECK_UNPHYSICAL_IN_FLUID
   Hydro_IsUnphysical_Single( q_L, "q_L", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
   Hydro_IsUnphysical_Single( q_R, "q_R", (real)0.0, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#  endif


// 2-2c. use the min/max of the left and right eigenvalues
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_DAVIS )
   const real W_L1 = u_L - Cs_L;
   const real W_L2 = u_R - Cs_R;
   const real W_R1 = u_L + Cs_L;
   const real W_R2 = u_R + Cs_R;
   W_L = FMIN( W_L1, W_L2 );
   W_R = FMAX( W_R1, W_R2 );


#  else
#  error : ERROR : unsupported HLLC_WAVESPEED !!
#  endif // HLLC_WAVESPEED


// 3. evaluate the star-region velocity (V_S) and pressure (P_S)
   real temp1_L, temp1_R, temp2, V_S, P_S;

// do not use u_L-W_L and u_R-W_R to prevent from large round-off errors when Cs<<u~W
// ==> temp1_L ~ temp1_R ~ 0.0 ==> temp2 = inf
#  if   ( HLLC_WAVESPEED == HLL_WAVESPEED_ROE )
   temp1_L = L[0]*(  (EVal_min<u_L-Cs_L) ? (u_L-EVal_min) : (+Cs_L)  );
   temp1_R = R[0]*(  (EVal_max>u_R+Cs_R) ? (u_R-EVal_max) : (-Cs_R)  );
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_PVRS )
   temp1_L = +L[0]*( Cs_L*q_L );
   temp1_R = -R[0]*( Cs_R*q_R );
#  elif ( HLLC_WAVESPEED == HLL_WAVESPEED_DAVIS )
   temp1_L = +L[0]*(  ( W_L1 < W_L2 ) ? Cs_L : (u_L-u_R)+Cs_R  );
   temp1_R = -R[0]*(  ( W_R2 > W_R1 ) ? Cs_R : (u_L-u_R)+Cs_L  );
#  else
#  error : ERROR : unsupported HLLC_WAVESPEED !!
#  endif

   temp2 = ONE / ( temp1_L - temp1_R );
   V_S   = temp2*( P_L - P_R + temp1_L*u_L - temp1_R*u_R );
   P_S   = temp2*(  temp1_L*( P_R + temp1_R*u_R ) - temp1_R*( P_L + temp1_L*u_L )  );
   P_S   = Hydro_CheckMinPres( P_S, MinPres );


// 4. evaluate the weightings of the left/right fluxes and contact wave
   real Flux_LR[NCOMP_TOTAL], temp4, Coeff_LR, Coeff_S;  // use NCOMP_TOTAL for Flux_LR since it will be passed to Hydro_Con2Flux()

   if ( V_S >= ZERO )
   {
      const real MaxV_L = FMIN( W_L, ZERO );

      Hydro_Con2Flux( 0, Flux_LR, L, MinPres, PassiveFloor, NULL, NULL, NULL, NULL, &P_L );

      for (int v=0; v<NCOMP_FLUID; v++)   Flux_LR[v] -= MaxV_L*L[v];    // fluxes along the maximum wave speed

//    deal with the special case of V_S=MaxV_L=0
//###REVISE: should it return zero flux due to symmetry?
      if ( V_S == ZERO  &&  MaxV_L == ZERO )
      {
         Coeff_LR = ONE;
         Coeff_S  = ZERO;
      }

      else
      {
         temp4    = ONE / ( V_S - MaxV_L );
         Coeff_LR = temp4*V_S;
         Coeff_S  = -temp4*MaxV_L*P_S;
      }
   } // if ( V_S >= ZERO )

   else // V_S < 0.0
   {
      const real MaxV_R = FMAX( W_R, ZERO );

      Hydro_Con2Flux( 0, Flux_LR, R, MinPres, PassiveFloor, NULL, NULL, NULL, NULL, &P_R );

      for (int v=0; v<NCOMP_FLUID; v++)    Flux_LR[v] -= MaxV_R*R[v];   // fluxes along the maximum wave speed

      temp4    = ONE / ( V_S - MaxV_R );
      Coeff_LR = temp4*V_S;
      Coeff_S  = -temp4*MaxV_R*P_S;
   } // if ( V_S >= ZERO ) ... else ...


// 5. evaluate the HLLC fluxes
   for (int v=0; v<NCOMP_FLUID; v++)   Flux_Out[v] = Coeff_LR*Flux_LR[v];

   Flux_Out[1] += Coeff_S;
   Flux_Out[4] += Coeff_S*V_S;


// 6. evaluate the fluxes of passive scalars
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

} // FUNCTION : Hydro_RiemannSolver_HLLC
#endif // #ifndef SRHD



#endif // #if ( MODEL == HYDRO )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLC__
