#ifndef __CUFLU_RIEMANNSOLVER_HLLE__
#define __CUFLU_RIEMANNSOLVER_HLLE__

#include "CUFLU.h"

// external function
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else

#include "../../../include/SRHydroPrototypes.h"

#endif

#if ( MODEL == SR_HYDRO )

//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_RiemannSolver_HLLE
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
//                [6] MinTemp  : Minimum allowed temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                                 const real Gamma, const real MinTemp )
{
  real CL[NCOMP_TOTAL], CR[NCOMP_TOTAL]; /* conserved vars. */
  real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL]; /* primitive vars. */
  real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
  real Fhll[NCOMP_TOTAL];
  
  const real Gamma_m1 = Gamma - (real)1.0;
  real rhl, rhr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  real ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  real lmdal,lmdar; /* Left and Right wave speeds */
  real ovlrmll;
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
   if ( ( ssl < 0.0 ) || ( ssr < 0.0 ) ) printf("ssl = %14.7e, ssr = %14.7e\n", ssl, ssr);
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
     Fl[1] = FMA(CL[1], lV1, PL[4]);
     Fl[2] = CL[2] * lV1;
     Fl[3] = CL[3] * lV1;
     Fl[4] = CL[1];

     Fr[0] = CR[0] * rV1;
     Fr[1] = FMA(CR[1], rV1, PR[4]);
     Fr[2] = CR[2] * rV1;
     Fr[3] = CR[3] * rV1;
     Fr[4] = CR[1];
 /* 7. Determine intercell flux according to Mignone 13
 */
  if( lmdal >= (real)0.0 ){ /* Fl */
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
  else if( lmdal < (real)0.0 && lmdar > (real)0.0 ){ /* Fs */
/* 5. Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 *    Compute HLL conserved quantities using Mignone eq 9
 */
    ovlrmll = (real)1.0 / ( lmdar - lmdal );
    lmdatlmda = lmdal*lmdar;

    Fhll[0] = FMA( lmdatlmda, (CR[0] - CL[0]), FMA( lmdar, Fl[0], - lmdal*Fr[0] ) ) * ovlrmll;
    Fhll[1] = FMA( lmdatlmda, (CR[1] - CL[1]), FMA( lmdar, Fl[1], - lmdal*Fr[1] ) ) * ovlrmll;
    Fhll[2] = FMA( lmdatlmda, (CR[2] - CL[2]), FMA( lmdar, Fl[2], - lmdal*Fr[2] ) ) * ovlrmll;
    Fhll[3] = FMA( lmdatlmda, (CR[3] - CL[3]), FMA( lmdar, Fl[3], - lmdal*Fr[3] ) ) * ovlrmll;

#   if ( CONSERVED_ENERGY == 1 )
    Fhll[4] = FMA( lmdatlmda, (CR[4] - CL[4]), FMA( lmdar, Fl[4], - lmdal*Fr[4] ) ) * ovlrmll;
#   elif ( CONSERVED_ENERGY == 2 )
    Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] + CR[0] - CL[4] - CL[0])) * ovlrmll;
#   endif

    /* calculate Fs */
    Flux_Out[0] = Fhll[0];
    Flux_Out[1] = Fhll[1];
    Flux_Out[2] = Fhll[2];
    Flux_Out[3] = Fhll[3];
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = Fhll[4];
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = Fhll[4] - Fhll[0];
#   endif

    SRHydro_Rotate3D( Flux_Out, XYZ, false );
    return;
  }
  else{ /* Fr */
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

} // FUNCTION : SRHydro_RiemannSolver_HLLE


#endif // #if ( MODEL == SR_HYDRO )
#endif
