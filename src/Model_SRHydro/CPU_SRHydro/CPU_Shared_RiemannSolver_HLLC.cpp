#ifndef __CUFLU_RIEMANNSOLVER_HLLC__
#define __CUFLU_RIEMANNSOLVER_HLLC__

#include "CUFLU.h"

// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
GPU_DEVICE
void QuadraticSolver (real A, real B, real C, real delta, real *x_plus, real *x_minus, const int line);

#else // #ifdef __CUDACC__

#include "../../../include/SRHydroPrototypes.h"
void QuadraticSolver (real A, real B, real C, real delta, real *x_plus, real *x_minus, const int line);

#endif // #ifdef __CUDACC__ ... else ...

#if ( MODEL == SR_HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_RiemannSolver_HLLC
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                The wave speed is estimated by the same formula in HLLE solver
//
// Note        :  1. The input data should be conserved variables
//                2. Ref : a. Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         b. Batten, P., Clarke, N., Lambert, C., & Causon, D. M. 1997, SIAM J. Sci. Comput.,
//                            18, 1553
//                3. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  [1] XYZ      : Target spatial direction : (0/1/2) --> (x/y/z)
//                [2] Flux_Out : Array to store the output flux
//                [3] L_In     : Input left  state (conserved variables)
//                [4] R_In     : Input right state (conserved variables)
//                [5] Gamma    : Ratio of specific heats
//                [6] MinTemp  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                 const real Gamma, const real MinTemp )
{
# if ( EOS ==  IDEAL_GAS)
  const real Gamma_m1 = Gamma - (real)1.0;
  real rhl, rhr;
# endif

# ifdef CHECK_NEGATIVE_IN_FLUID
  real lV2, lV3, rV2, rV3;
# endif
  real CL[NCOMP_TOTAL], CR[NCOMP_TOTAL]; /* conserved vars. */
  real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL]; /* primitive vars. */
  real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
  real Fhll[NCOMP_TOTAL], Uhll[NCOMP_TOTAL];
  real Usl[NCOMP_TOTAL], Usr[NCOMP_TOTAL];
  real cslsq, csrsq, gammasql, gammasqr;
  real ssl, ssr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  real lmdal,lmdar; /* Left and Right wave speeds */
  real lmdas; /* Contact wave speed */
  real ovlrmll;
  real a,b,c;
  real den,ps; /* Pressure in inner region */
  real lV1, rV1;
  real lFactor,rFactor; /* Lorentz factor */

/* 0. reorder the input conserved variables for different spatial directions */
   for(int v=0;v<NCOMP_TOTAL;v++){
       CL[v]=L_In[v];
       CR[v]=R_In[v];
   }

   SRHydro_Rotate3D( CL, XYZ, true );
   SRHydro_Rotate3D( CR, XYZ, true );

/* 1. compute primitive vars. from conserved vars. */
   SRHydro_Con2Pri (CL, PL, Gamma, MinTemp);
   SRHydro_Con2Pri (CR, PR, Gamma, MinTemp);

/* 2. Transform 4-velocity to 3-velocity */
   lFactor=1/SQRT(1+VectorDotProduct(PL, PL, 1, 3));
   rFactor=1/SQRT(1+VectorDotProduct(PR, PR, 1, 3));

   lV1=PL[1]*lFactor;

   rV1=PR[1]*rFactor;

#  ifdef CHECK_NEGATIVE_IN_FLUID
   lV2=PL[2]*lFactor;
   lV3=PL[3]*lFactor;

   rV2=PR[2]*rFactor;
   rV3=PR[3]*rFactor;

   real lV, rV;
   lV = SQRT(lV1*lV1 + lV2*lV2 + lV3*lV3);
   rV = SQRT(rV1*rV1 + rV2*rV2 + rV3*rV3);
  
   if ( lV >= 1.0 || rV >= 1.0 ) {
     printf( "function: %s: %d\n", __FUNCTION__, __LINE__);
     printf( "lV = %20.17e, rV = %20.17e\n", lV, rV);
     printf( "lUx = %20.17e, lUy = %20.17e, lUz = %20.17e\n", PL[1], PL[2], PL[3]);
     printf( "rUx = %20.17e, rUy = %20.17e, rUz = %20.17e\n", PR[1], PR[2], PR[3]);
   }
#  endif


/* 3. Compute the max and min wave speeds used in Mignone */
#  if ( EOS == RELATIVISTIC_IDEAL_GAS )
   //real nhl =  2.5*PL[4] + SQRT(2.25*SQR(PL[4]) + SQR(PL[0]));
   //real nhr =  2.5*PR[4] + SQRT(2.25*SQR(PR[4]) + SQR(PR[0]));
   real nhl =  FMA( 2.5, PL[4], SQRT( FMA( 2.25, SQR(PL[4]), SQR(PL[0]) ) ) );
   real nhr =  FMA( 2.5, PR[4], SQRT( FMA( 2.25, SQR(PR[4]), SQR(PR[0]) ) ) );

//   cslsq = ( PL[4]*( 5*nhl - 8*PL[4] ) ) / ((3*nhl)*( nhl - PL[4] ));
//   csrsq = ( PR[4]*( 5*nhr - 8*PR[4] ) ) / ((3*nhr)*( nhr - PR[4] ));

   cslsq = PL[4] * FMA( 4.5, PL[4], 5*SQRT( FMA( 2.25, SQR(PL[4]), SQR(PL[0]) ) ) ) 
        / ( 3*nhl* FMA( 1.5, PL[4],   SQRT( FMA( 2.25, SQR(PL[4]), SQR(PL[0]) ) ) ) );

   csrsq = PR[4] * FMA( 4.5, PR[4], 5*SQRT( FMA( 2.25, SQR(PR[4]), SQR(PR[0]) ) ) ) 
        / ( 3*nhr* FMA( 1.5, PR[4],   SQRT( FMA( 2.25, SQR(PR[4]), SQR(PR[0]) ) ) ) );

#  elif ( EOS ==  IDEAL_GAS)
   rhl = PL[0] + PL[4] * Gamma / Gamma_m1; /* Mignone Eq 3.5 */
   rhr = PR[0] + PR[4] * Gamma / Gamma_m1;

   cslsq = Gamma * PL[4] / rhl; /* Mignone Eq 4 */
   csrsq = Gamma * PR[4] / rhr;
#  endif

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( cslsq >= 1.0 || csrsq >= 1.0 || cslsq < 0.0 || csrsq < 0.0 )
     printf( "cslsq=%10.7e, cslrq=%10.7e\n", cslsq, csrsq);
#  endif


// square of Lorentz factor
   gammasql = 1.0 + VectorDotProduct(PL, PL, 1, 3);
   gammasqr = 1.0 + VectorDotProduct(PR, PR, 1, 3);

   ssl = cslsq / FMA( - gammasql, cslsq, gammasql ); /* Mignone Eq 22.5 */
   ssr = csrsq / FMA( - gammasqr, csrsq, gammasqr ); /* Mignone Eq 22.5 */

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( ( ssl < 0.0 ) || ( ssr < 0.0 ) ) printf("ssl = %14.7e, ssr = %14.7e\n", ssl, ssr);
#  endif

   real lV1s = lV1*lV1;
   real rV1s = rV1*rV1;
 
   real deltal = (1.0 - lV1s) + ssl;
   real deltar = (1.0 - rV1s) + ssr;

   QuadraticSolver(1.0 + ssl, -2*lV1, lV1s - ssl, deltal, &lmdapl, &lmdaml, __LINE__);
   QuadraticSolver(1.0 + ssr, -2*rV1, rV1s - ssr, deltar, &lmdapr, &lmdamr, __LINE__);

   lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
   lmdar = MAX(lmdapl, lmdapr);

    
/* 4. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18) 
 *    compute HLL conserved quantities using Mignone eq 9
 * */
   Fl[0] = CL[0] * lV1;
   Fl[1] = FMA( CL[1], lV1, PL[4] );
   Fl[2] = CL[2] * lV1;
   Fl[3] = CL[3] * lV1;
   Fl[4] = CL[1];

  if( lmdal >= 0.0){ /* Fl */
    /* intercell flux is left flux */
    Flux_Out[0] = Fl[0];
    Flux_Out[1] = Fl[1];
    Flux_Out[2] = Fl[2];
    Flux_Out[3] = Fl[3];
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = Fl[4];
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = Fl[4] - Fl[0];
#   endif

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

   Fr[0] = CR[0] * rV1;
   Fr[1] = FMA( CR[1], rV1, PR[4] );
   Fr[2] = CR[2] * rV1;
   Fr[3] = CR[3] * rV1;
   Fr[4] = CR[1];

   if( lmdar <= 0.0 ){ /* Fr */
    /* intercell flux is right flux */
    Flux_Out[0] = Fr[0];
    Flux_Out[1] = Fr[1];
    Flux_Out[2] = Fr[2];
    Flux_Out[3] = Fr[3];
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = Fr[4];
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = Fr[4] - Fr[0];
#   endif

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

/* 5. Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 *    Compute HLL conserved quantities using Mignone eq 9
 */
  ovlrmll = 1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;


  Fhll[0] = FMA( lmdatlmda, (CR[0] - CL[0]), FMA( lmdar, Fl[0], - lmdal*Fr[0] ) ) * ovlrmll;
  Fhll[1] = FMA( lmdatlmda, (CR[1] - CL[1]), FMA( lmdar, Fl[1], - lmdal*Fr[1] ) ) * ovlrmll;
  Fhll[2] = FMA( lmdatlmda, (CR[2] - CL[2]), FMA( lmdar, Fl[2], - lmdal*Fr[2] ) ) * ovlrmll;
  Fhll[3] = FMA( lmdatlmda, (CR[3] - CL[3]), FMA( lmdar, Fl[3], - lmdal*Fr[3] ) ) * ovlrmll;
# if ( CONSERVED_ENERGY == 1 )
  Fhll[4] = FMA( lmdatlmda, (CR[4] - CL[4]), FMA( lmdar, Fl[4], - lmdal*Fr[4] ) ) * ovlrmll;
# elif ( CONSERVED_ENERGY == 2 )
  Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] + CR[0] - CL[4] - CL[0])) * ovlrmll;
# endif

  Uhll[0] = ( FMA( lmdar, CR[0], Fl[0] ) - FMA( lmdal, CL[0], Fr[0] ) ) * ovlrmll;
  Uhll[1] = ( FMA( lmdar, CR[1], Fl[1] ) - FMA( lmdal, CL[1], Fr[1] ) ) * ovlrmll;
  Uhll[2] = ( FMA( lmdar, CR[2], Fl[2] ) - FMA( lmdal, CL[2], Fr[2] ) ) * ovlrmll;
  Uhll[3] = ( FMA( lmdar, CR[3], Fl[3] ) - FMA( lmdal, CL[3], Fr[3] ) ) * ovlrmll;
# if ( CONSERVED_ENERGY == 1 )
  Uhll[4] = ( FMA( lmdar, CR[4], Fl[4] ) - FMA( lmdal, CL[4], Fr[4] ) ) * ovlrmll;
# elif ( CONSERVED_ENERGY == 2 )
  Uhll[4] = (lmdar * ( CR[4] + CR[0] ) - lmdal * ( CL[4] + CL[0]) + Fl[4] - Fr[4]) * ovlrmll;
# endif

# ifdef CHECK_NEGATIVE_IN_FLUID
  SRHydro_CheckUnphysical(Uhll, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
# endif

/* 6. Compute contact wave speed using larger root from Mignone Eq 18
 *    Physical root is the root with the minus sign
 */
  /* quadratic formuLa calcuLation */
  a = Fhll[4];
  b = -(Uhll[4] + Fhll[1]);
  c = Uhll[1];

  real delta = b*b-4*a*c;
  real null;

  QuadraticSolver(a, b ,c, delta, &null, &lmdas, __LINE__);


 /* 7. Determine intercell flux according to Mignone 13
 */
    if( lmdas >= 0.0 ){ /* Fls */

    /* Mignone 2006 Eq 48 */
    ps = FMA( -Fhll[4], lmdas, Fhll[1]);

    /* now calculate Usl with Mignone Eq 16 */
    den = 1.0 / (lmdal - lmdas);

    real factor0 = lmdal - lV1;
    real factor1 = FMA( lmdal, den, -lV1*den );

    Usl[0] =  CL[0] * factor1;
    Usl[1] = FMA( CL[1], factor0, ps - PL[4] )* den;
    Usl[2] =  CL[2] * factor1;
    Usl[3] =  CL[3] * factor1;
#   if ( CONSERVED_ENERGY == 1 )
    Usl[4] = FMA( - PL[4], lV1, FMA( CL[4], factor0, ps * lmdas ) ) * den;
#   elif ( CONSERVED_ENERGY == 2 )
    Usl[4] = (( CL[4] + CL[0] ) * factor0 + ps * lmdas - PL[4] * lV1) * den;
#   endif

#   ifdef CHECK_NEGATIVE_IN_FLUID
    SRHydro_CheckUnphysical(Usl, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdal, Usl[0] - CL[0], Fl[0] );
    Flux_Out[1] = FMA( lmdal, Usl[1] - CL[1], Fl[1] );
    Flux_Out[2] = FMA( lmdal, Usl[2] - CL[2], Fl[2] );
    Flux_Out[3] = FMA( lmdal, Usl[3] - CL[3], Fl[3] );
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = FMA( lmdal, Usl[4] - CL[4], Fl[4] );
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = lmdal*(Usl[4] - CL[4] - CL[0]) + Fl[4] - Flux_Out[0];
#   endif
    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }
  else{ /* Frs */
    /* Mignone 2006 Eq 48 */
    ps = FMA( -Fhll[4], lmdas, Fhll[1] );
    /* now calculate Usr with Mignone Eq 16 */
    den = 1.0 / (lmdar - lmdas);

    real factor0 = lmdar - rV1;
    real factor1 = FMA( lmdar, den, -rV1*den );

    Usr[0] = CR[0] * factor1;
    Usr[1] = FMA( CR[1], factor0, ps - PR[4] ) * den;
    Usr[2] = CR[2] * factor1;
    Usr[3] = CR[3] * factor1;
#   if ( CONSERVED_ENERGY == 1 )
    Usr[4] = FMA( - PR[4], rV1, FMA( CR[4], factor0, ps * lmdas ) ) * den;
#   elif ( CONSERVED_ENERGY == 2 )
    Usr[4] = (( CR[4] + CR[0] ) * factor0 + ps * lmdas - PR[4] * rV1) * den;
#   endif
#   ifdef CHECK_NEGATIVE_IN_FLUID
    SRHydro_CheckUnphysical(Usr, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdar, Usr[0] - CR[0], + Fr[0] );
    Flux_Out[1] = FMA( lmdar, Usr[1] - CR[1], + Fr[1] );
    Flux_Out[2] = FMA( lmdar, Usr[2] - CR[2], + Fr[2] );
    Flux_Out[3] = FMA( lmdar, Usr[3] - CR[3], + Fr[3] );
    Flux_Out[4] = FMA( lmdar, Usr[4] - CR[4], + Fr[4] );
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = FMA( lmdar, Usr[4] - CR[4], + Fr[4] );
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = lmdar*(Usr[4] - CR[4] - CR[0]) + Fr[4] - Flux_Out[0];
#   endif

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

} // FUNCTION : SRHydro_RiemannSolver_HLLC



#endif // #if ( MODEL == SR_HYDRO )
#endif
