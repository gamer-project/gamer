#ifndef __CUFLU_RIEMANNSOLVER_HLLC__
#define __CUFLU_RIEMANNSOLVER_HLLC__

#include "CUFLU.h"

// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

#include "../../../include/SRHydroPrototypes.h"

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
  real nhl, nhr;
# if ( EOS ==  CONSTANT_GAMMA)
  const real Gamma_m1 = Gamma - (real)1.0;
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
  real lV1, rV1, lV2, rV2, lV3, rV3;
  real lFactor,rFactor; /* Lorentz factor */

/* 0. reorder the input conserved variables for different spatial directions */
   for(int v=0;v<NCOMP_TOTAL;v++){
       CL[v]=L_In[v];
       CR[v]=R_In[v];
   }

   SRHydro_Rotate3D( CL, XYZ, true );
   SRHydro_Rotate3D( CR, XYZ, true );

/* 1. compute primitive vars. from conserved vars. */
   lFactor = SRHydro_Con2Pri (CL, PL, Gamma, MinTemp);
   rFactor = SRHydro_Con2Pri (CR, PR, Gamma, MinTemp);

#  ifdef CHECK_FAILED_CELL_IN_FLUID
   SRHydro_CheckUnphysical(NULL, PL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
   SRHydro_CheckUnphysical(NULL, PR, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#  endif

/* 2. Transform 4-velocity to 3-velocity */
   lV1=PL[1]/lFactor;
   lV2=PL[2]/lFactor;
   lV3=PL[3]/lFactor;

   rV1=PR[1]/rFactor;
   rV2=PR[2]/rFactor;
   rV3=PR[3]/rFactor;



/* 3. Compute the max and min wave speeds used in Mignone */
   cslsq = SoundSpeedSquare( PL[4]/PL[0], Gamma);
   csrsq = SoundSpeedSquare( PR[4]/PR[0], Gamma);

#  ifdef CHECK_FAILED_CELL_IN_FLUID
   if ( cslsq >= 1.0 || csrsq >= 1.0 || cslsq < 0.0 || csrsq < 0.0 )
     printf( "cslsq=%10.7e, cslrq=%10.7e\n", cslsq, csrsq);
#  endif


// square of Lorentz factor
   gammasql = SQR(lFactor);
   gammasqr = SQR(rFactor);

   ssl = cslsq / FMA( - gammasql, cslsq, gammasql ); /* Mignone Eq 22.5 */
   ssr = csrsq / FMA( - gammasqr, csrsq, gammasqr ); /* Mignone Eq 22.5 */

#  ifdef CHECK_FAILED_CELL_IN_FLUID
   if ( ( ssl < (real)0.0 ) || ( ssr < (real)0.0 ) ) printf("ssl = %14.7e, ssr = %14.7e\n", ssl, ssr);
#  endif

   real lV2s = lV2*lV2;
   real rV2s = rV2*rV2;

   real lV3s = lV3*lV3;
   real rV3s = rV3*rV3;
 
   real __gammasql = (real)1.0 / gammasql;
   real __gammasqr = (real)1.0 / gammasqr;

   real deltal = ssl*ssl + ssl*( __gammasql + lV2s + lV3s );
   real deltar = ssr*ssr + ssr*( __gammasqr + rV2s + rV3s );

   real ssl__ = (real)1.0 + ssl;
   real ssr__ = (real)1.0 + ssr;


   lmdapl = ( lV1 + SQRT(deltal) ) / ssl__ ;
   lmdaml = ( lV1 - SQRT(deltal) ) / ssl__ ;

   lmdapr = ( rV1 + SQRT(deltar) ) / ssr__ ;
   lmdamr = ( rV1 - SQRT(deltar) ) / ssr__ ;

   lmdal = FMIN(lmdaml, lmdamr); /* Mignone Eq 21 */
   lmdar = FMAX(lmdapl, lmdapr);

    
/* 4. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18) 
 *    compute HLL conserved quantities using Mignone eq 9
 * */
   Fl[0] = CL[0] * lV1;
   Fl[1] = FMA( CL[1], lV1, PL[4] );
   Fl[2] = CL[2] * lV1;
   Fl[3] = CL[3] * lV1;
#  if ( CONSERVED_ENERGY == 1 )
   Fl[4] = CL[1];
#  elif ( CONSERVED_ENERGY == 2 )
   Fl[4] = ( CL[4] + PL[4] ) * lV1;
#  endif

  if( lmdal >= (real)0.0){ /* Fl */
    /* intercell flux is left flux */
    Flux_Out[0] = Fl[0];
    Flux_Out[1] = Fl[1];
    Flux_Out[2] = Fl[2];
    Flux_Out[3] = Fl[3];
    Flux_Out[4] = Fl[4];

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

   Fr[0] = CR[0] * rV1;
   Fr[1] = FMA( CR[1], rV1, PR[4] );
   Fr[2] = CR[2] * rV1;
   Fr[3] = CR[3] * rV1;
#  if ( CONSERVED_ENERGY == 1 )
   Fr[4] = CR[1];
#  elif ( CONSERVED_ENERGY == 2 )
   Fr[4] = ( CR[4] + PR[4] ) * rV1;
#  endif

   if( lmdar <= (real)0.0 ){ /* Fr */
    /* intercell flux is right flux */
    Flux_Out[0] = Fr[0];
    Flux_Out[1] = Fr[1];
    Flux_Out[2] = Fr[2];
    Flux_Out[3] = Fr[3];
    Flux_Out[4] = Fr[4];

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

/* 5. Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 *    Compute HLL conserved quantities using Mignone eq 9
 */
  ovlrmll = (real)1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;


  Fhll[0] = FMA( lmdatlmda, (CR[0] - CL[0]), FMA( lmdar, Fl[0], - lmdal*Fr[0] ) ) * ovlrmll;
  Fhll[1] = FMA( lmdatlmda, (CR[1] - CL[1]), FMA( lmdar, Fl[1], - lmdal*Fr[1] ) ) * ovlrmll;
  Fhll[2] = FMA( lmdatlmda, (CR[2] - CL[2]), FMA( lmdar, Fl[2], - lmdal*Fr[2] ) ) * ovlrmll;
  Fhll[3] = FMA( lmdatlmda, (CR[3] - CL[3]), FMA( lmdar, Fl[3], - lmdal*Fr[3] ) ) * ovlrmll;
  Fhll[4] = FMA( lmdatlmda, (CR[4] - CL[4]), FMA( lmdar, Fl[4], - lmdal*Fr[4] ) ) * ovlrmll;

  Uhll[0] = ( FMA( lmdar, CR[0], Fl[0] ) - FMA( lmdal, CL[0], Fr[0] ) ) * ovlrmll;
  Uhll[1] = ( FMA( lmdar, CR[1], Fl[1] ) - FMA( lmdal, CL[1], Fr[1] ) ) * ovlrmll;
  Uhll[2] = ( FMA( lmdar, CR[2], Fl[2] ) - FMA( lmdal, CL[2], Fr[2] ) ) * ovlrmll;
  Uhll[3] = ( FMA( lmdar, CR[3], Fl[3] ) - FMA( lmdal, CL[3], Fr[3] ) ) * ovlrmll;
  Uhll[4] = ( FMA( lmdar, CR[4], Fl[4] ) - FMA( lmdal, CL[4], Fr[4] ) ) * ovlrmll;

# ifdef CHECK_FAILED_CELL_IN_FLUID
  SRHydro_CheckUnphysical(Uhll, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
# endif

/* 6. Compute contact wave speed using larger root from Mignone Eq 18
 *    Physical root is the root with the minus sign
 */
  /* quadratic formuLa calcuLation */
  a = Fhll[4];
  b = -(Uhll[4] + Fhll[1]);
  c = Uhll[1];

  real delta = FMA( b, b, -(real)4*a*c );

# ifdef CHECK_FAILED_CELL_IN_FLUID
  if (delta < (real) 0.0) printf("delta=%f\n", delta);
# endif

    lmdas = - ((real)2.0 * c) / ( b + SIGN(b) * SQRT( delta ) );


 /* 7. Determine intercell flux according to Mignone 13
 */
    if( lmdas >= (real)0.0 ){ /* Fls */

    /* Mignone 2006 Eq 48 */
    ps = FMA( -Fhll[4], lmdas, Fhll[1]);

    /* now calculate Usl with Mignone Eq 16 */
    den = (real)1.0 / (lmdal - lmdas);

    real factor0 = lmdal - lV1;
    real factor1 = FMA( lmdal, den, -lV1*den );

    Usl[0] =  CL[0] * factor1;
    Usl[1] = FMA( CL[1], factor0, ps - PL[4] )* den;
    Usl[2] =  CL[2] * factor1;
    Usl[3] =  CL[3] * factor1;
    Usl[4] = FMA( - PL[4], lV1, FMA( CL[4], factor0, ps * lmdas ) ) * den;

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    SRHydro_CheckUnphysical(Usl, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdal, Usl[0] - CL[0], Fl[0] );
    Flux_Out[1] = FMA( lmdal, Usl[1] - CL[1], Fl[1] );
    Flux_Out[2] = FMA( lmdal, Usl[2] - CL[2], Fl[2] );
    Flux_Out[3] = FMA( lmdal, Usl[3] - CL[3], Fl[3] );
    Flux_Out[4] = FMA( lmdal, Usl[4] - CL[4], Fl[4] );

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }
  else{ /* Frs */
    /* Mignone 2006 Eq 48 */
    ps = FMA( -Fhll[4], lmdas, Fhll[1] );
    /* now calculate Usr with Mignone Eq 16 */
    den = (real)1.0 / (lmdar - lmdas);

    real factor0 = lmdar - rV1;
    real factor1 = FMA( lmdar, den, -rV1*den );

    Usr[0] = CR[0] * factor1;
    Usr[1] = FMA( CR[1], factor0, ps - PR[4] ) * den;
    Usr[2] = CR[2] * factor1;
    Usr[3] = CR[3] * factor1;
    Usr[4] = FMA( - PR[4], rV1, FMA( CR[4], factor0, ps * lmdas ) ) * den;

#   ifdef CHECK_FAILED_CELL_IN_FLUID
    SRHydro_CheckUnphysical(Usr, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#   endif

    /* now calculate Fsr using Mignone Eq 14 */
    Flux_Out[0] = FMA( lmdar, Usr[0] - CR[0], + Fr[0] );
    Flux_Out[1] = FMA( lmdar, Usr[1] - CR[1], + Fr[1] );
    Flux_Out[2] = FMA( lmdar, Usr[2] - CR[2], + Fr[2] );
    Flux_Out[3] = FMA( lmdar, Usr[3] - CR[3], + Fr[3] );
    Flux_Out[4] = FMA( lmdar, Usr[4] - CR[4], + Fr[4] );

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

} // FUNCTION : SRHydro_RiemannSolver_HLLC



#endif // #if ( MODEL == SR_HYDRO )
#endif
