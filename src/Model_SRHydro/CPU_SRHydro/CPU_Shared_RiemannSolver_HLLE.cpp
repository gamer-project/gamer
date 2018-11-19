#include "GAMER.h"
#include "CUFLU.h"
#include "../../../include/CPU_prototypes.h"

#if ( MODEL == SR_HYDRO )

static void QuadraticSolver (real A, real B, real C, real *x_plus, real *x_minus);

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

/* 1. compute primitive vars. from conserved vars. */
   CPU_Con2Pri (CL, PL, Gamma);
   CPU_Con2Pri (CR, PR, Gamma);

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
#  if ( EOS == RELATIVISTIC_IDEAL_GAS )
   real nhl =  2.5*PL[4] + SQRT(2.25*SQR(PL[4]) + SQR(PL[0]));
   real nhr =  2.5*PR[4] + SQRT(2.25*SQR(PR[4]) + SQR(PR[0]));

   cslsq = ( PL[4]*( 5*nhl - 8*PL[4] ) ) / ((3*nhl)*( nhl - PL[4] ));
   csrsq = ( PR[4]*( 5*nhr - 8*PR[4] ) ) / ((3*nhr)*( nhr - PR[4] ));
#  elif ( EOS ==  IDEAL_GAS)
   rhl = PL[0] + PL[4] * Gamma / Gamma_m1; /* Mignone Eq 3.5 */
   rhr = PR[0] + PR[4] * Gamma / Gamma_m1;

   cslsq = Gamma * PL[4] / rhl; /* Mignone Eq 4 */
   csrsq = Gamma * PR[4] / rhr;
#  endif

// square of Lorentz factor
   gammasql = 1.0 + SQR(PL[1]) + SQR(PL[2]) + SQR(PL[3]);
   gammasqr = 1.0 + SQR(PR[1]) + SQR(PR[2]) + SQR(PR[3]);

   ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
   ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

   QuadraticSolver(1.0 + ssl, -2*lV1, lV1*lV1 - ssl, &lmdapl, &lmdaml);
   QuadraticSolver(1.0 + ssr, -2*rV1, rV1*rV1 - ssr, &lmdapr, &lmdamr);

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
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = Fl[4];
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = Fl[4] - Fl[0];
#   endif

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
#   if ( CONSERVED_ENERGY == 1 )
    Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] - CL[4])) * ovlrmll;
#   elif ( CONSERVED_ENERGY == 2 )
    Fhll[4] = (lmdar*Fl[4] - lmdal*Fr[4] + lmdatlmda * (CR[4] + CR[0] - CL[4] - CL[0])) * ovlrmll;
#   endif

    /* calcULate Fs  */
    Flux_Out[0] = Fhll[0];
    Flux_Out[1] = Fhll[1];
    Flux_Out[2] = Fhll[2];
    Flux_Out[3] = Fhll[3];
#   if ( CONSERVED_ENERGY == 1 )
    Flux_Out[4] = Fhll[4];
#   elif ( CONSERVED_ENERGY == 2 )
    Flux_Out[4] = Fhll[4] - Fhll[0];
#   endif

   CPU_Rotate3D( Flux_Out, XYZ, false );
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

   CPU_Rotate3D( Flux_Out, XYZ, false );
    return;
  }

} // FUNCTION : CPU_RiemannSolver_HLLE


//=====================================================
// Solve A*X^2 + B*x + C = 0
// delta = sqrt(B*B-4*A*C)
// x_plus  = ( -B + delta ) / (2*A)
// x_minus = ( -B - delta ) / (2*A)
//=====================================================
static void QuadraticSolver (real A, real B, real C, real *x_plus, real *x_minus)
{
  real delta = B*B-4*A*C;

  if (A != 0.0){
           if ( delta > 0.0 ) {

           real factor = -0.5*( B + SIGN(B) *  SQRT(delta) );
     
           if  ( B > 0.0 && C != 0.0 ){
             *x_plus   = C/factor;
     	     *x_minus  = factor/A;             return;
          }else if  ( B < 0.0 && C != 0.0 ){
     	     *x_plus   = factor/A;
     	     *x_minus  = C/factor;             return;
          }else if ( B == 0.0 && C < 0.0 ){
             *x_plus = SQRT(-C/A);
             *x_minus = -SQRT(-C/A);           return;
          }else if ( B > 0.0 && C == 0.0 ){
             *x_plus = 0.0;
             *x_minus = -B/A;                  return;
          }else if ( B < 0.0 && C == 0.0 ){
             *x_plus = -B/A;
             *x_minus = 0.0;                   return;
          }else if ( B == 0.0 && C == 0.0 ){
             *x_plus  = 0.0;
             *x_minus = 0.0;                   return;
          }else                                goto NO_REAL_SOLUTIONS;

     }else if ( delta == 0.0 ){
             *x_plus  = -0.5*B/A;
             *x_minus = -0.5*B/A;               return;
     }else                                      goto NO_REAL_SOLUTIONS;
  }else{ // if ( A == 0.0 )
        if ( B != 0.0 ){
	   *x_plus  = -C/B;
	   *x_minus = -C/B;                     return;
         }else                                  goto NO_REAL_SOLUTIONS;
  }

     NO_REAL_SOLUTIONS:
     {
        Aux_Message(stderr, "No real solution in Quadratic Solver!\n");
        Aux_Message(stderr, "A=%14.7e, B=%14.7e, C=%14.7e\n", A, B, C);
        Aux_Message(stderr, "B*B-4*A*C=%14.7e\n", B*B-4*A*C);  return;
     }
}

#endif // #if ( MODEL == SR_HYDRO )
