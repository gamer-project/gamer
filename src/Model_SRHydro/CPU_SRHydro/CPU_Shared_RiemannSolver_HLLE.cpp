#include "GAMER.h"
#include "CUFLU.h"
#include "../../../include/CPU_prototypes.h"

#if ( MODEL == SR_HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannSolver_HLLE
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
//                [6] MinPres  : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], 
                             const real Gamma, const real MinPres )
{
  real CL[NCOMP_TOTAL], CR[NCOMP_TOTAL]; /* conserved vars. */
  real PL[NCOMP_TOTAL], PR[NCOMP_TOTAL]; /* primitive vars. */
  real Fl[NCOMP_TOTAL], Fr[NCOMP_TOTAL];
  real Fhll[NCOMP_TOTAL];
  
  const real Gamma_m1 = Gamma - (real)1.0;
  double rhl, rhr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  double ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  double lmdal,lmdar; /* Left and Right wave speeds */
  double ovlrmll;
  double lV1, lV2, lV3, rV1, rV2, rV3;
  double lFactor,rFactor; /* Lorentz factor */

/* 0. reorder the input conserved variables for different spatial directions */
   for(int v=0;v<NCOMP_TOTAL;v++){
       CL[v]=L_In[v];
       CR[v]=R_In[v];
   }


   CPU_Rotate3D( CL, XYZ, true );
   CPU_Rotate3D( CR, XYZ, true );

/* 0.5 check negative conserved quanties */


/* 1. compute primitive vars. from conserved vars. */
   CPU_Con2Pri (CL, PL, Gamma);
   CPU_Con2Pri (CR, PR, Gamma);

/*  1.4 check negative primitive quanties*/



/* 2. Transform 4-velocity to 3-velocity */
   lFactor=1/SQRT(1+SQR(PL[1])+SQR(PL[2])+SQR(PL[3]));
   rFactor=1/SQRT(1+SQR(PR[1])+SQR(PR[2])+SQR(PR[3]));

   lV1=PL[1]*lFactor;
   lV2=PL[2]*lFactor;
   lV3=PL[3]*lFactor;

   rV1=PR[1]*rFactor;
   rV2=PR[2]*rFactor;
   rV3=PR[3]*rFactor;

/* 3. Compute the max and min wave speeds used in Mignone */
   rhl = PL[0] + PL[4] * Gamma / Gamma_m1; /* Mignone Eq 3.5 */
   rhr = PR[0] + PR[4] * Gamma / Gamma_m1;

   cslsq = Gamma * PL[4] / rhl; /* Mignone Eq 4 */
   csrsq = Gamma * PR[4] / rhr;

   vsql = SQR(lV1) + SQR(lV2) + SQR(lV3);
   vsqr = SQR(rV1) + SQR(rV2) + SQR(rV3);

   gammasql = 1.0 / (1.0 - vsql);
   gammasqr = 1.0 / (1.0 - vsqr);

   ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
   ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ((ssl*(1.0-SQR(lV1)+ssl) < 0))
      {
	Aux_Message (stderr, "\n\nerror:%s: %d\n",__FUNCTION__,__LINE__);
	Aux_Message (stderr, "ssl=%e\n",ssl);
	Aux_Message (stderr, "ssl*(1.0-SQR(lV1)+ssl = %e < 0.\n",ssl*(1.0-SQR(lV1)+ssl));
      }
   if ((ssr*(1.0-SQR(rV1)+ssr) < 0))
      {
	Aux_Message (stderr, "\n\nerror:%s: %d\n",__FUNCTION__,__LINE__);
	Aux_Message (stderr, "ssr=%e\n",ssr);
	Aux_Message (stderr, "ssr*(1.0-SQR(rV1)+ssr = %e < 0.\n",ssr*(1.0-SQR(rV1)+ssr));
      }
#  endif

   radl = SQRT( ssl*(1.0-SQR(lV1)+ssl) ); /* Mignone Eq 23 (radical part) */
   radr = SQRT( ssr*(1.0-SQR(rV1)+ssr) );

   lmdapl = (lV1 + radl) / (1.0 + ssl); /* Mignone Eq 23 */
   lmdapr = (rV1 + radr) / (1.0 + ssr);
   lmdaml = (lV1 - radl) / (1.0 + ssl);
   lmdamr = (rV1 - radr) / (1.0 + ssr);

   lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
   lmdar = MAX(lmdapl, lmdapr);
    
/* 4. compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18) 
 *    compute HLL conserved quantities using Mignone eq 9
 * */
     Fl[0] = CL[0] * lV1;
     Fl[1] = CL[1] * lV1 + PL[4];
     Fl[2] = CL[2] * lV1;
     Fl[3] = CL[3] * lV1;
     Fl[4] = CL[1];

     Fr[0] = CR[0] * rV1;
     Fr[1] = CR[1] * rV1 + PR[4];
     Fr[2] = CR[2] * rV1;
     Fr[3] = CR[3] * rV1;
     Fr[4] = CR[1];
 /* 7. Determine intercell flux according to Mignone 13
 */
  if( lmdal >= 0.0){ /* Fl */
    /* intercell flux is left flux */
    Flux_Out[0] = Fl[0];
    Flux_Out[1] = Fl[1];
    Flux_Out[2] = Fl[2];
    Flux_Out[3] = Fl[3];
    Flux_Out[4] = Fl[4];

   CPU_Rotate3D( Flux_Out, XYZ, false );
   return;
  }
  else if( lmdal <= 0.0 && lmdar >= 0.0 ){ /* Fs */
/* 5. Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 *    Compute HLL conserved quantities using Mignone eq 9
 */
    ovlrmll = 1.0 / ( lmdar - lmdal );
    lmdatlmda = lmdal*lmdar;

    Fhll[0] = (lmdar*Fl[0] - lmdal*Fr[0] + lmdatlmda * (CR[0] - CL[0])) * ovlrmll;
    Fhll[1] = (lmdar*Fl[1] - lmdal*Fr[1] + lmdatlmda * (CR[1] - CL[1])) * ovlrmll;
    Fhll[2] = (lmdar*Fl[2] - lmdal*Fr[2] + lmdatlmda * (CR[2] - CL[2])) * ovlrmll;
    Fhll[3] = (lmdar*Fl[3] - lmdal*Fr[3] + lmdatlmda * (CR[3] - CL[3])) * ovlrmll;
    Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] - CL[4])) * ovlrmll;

    /* calcULate Fs  */
    Flux_Out[0] = Fhll[0];
    Flux_Out[1] = Fhll[1];
    Flux_Out[2] = Fhll[2];
    Flux_Out[3] = Fhll[3];
    Flux_Out[4] = Fhll[4];

   CPU_Rotate3D( Flux_Out, XYZ, false );
    return;
  }
  else{ /* Fr */
    /* intercell flux is right flux */

    Flux_Out[0] = Fr[0];
    Flux_Out[1] = Fr[1];
    Flux_Out[2] = Fr[2];
    Flux_Out[3] = Fr[3];
    Flux_Out[4] = Fr[4];

   CPU_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

} // FUNCTION : CPU_RiemannSolver_HLLE



#endif // #if ( MODEL == SR_HYDRO )
