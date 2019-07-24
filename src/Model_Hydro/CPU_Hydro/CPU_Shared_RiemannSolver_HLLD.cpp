#ifndef __CUFLU_RIEMANNSOLVER_HLLD__
#define __CUFLU_RIEMANNSOLVER_HLLD__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined MHD )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_Con2Flux( const int XYZ, real Flux[], const real In[], const real Gamma_m1, const real MinPres );
real Hydro_CheckMinPres( const real InPres, const real MinPres );
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );

#endif // #ifdef __CUDACC__ ... else ...



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannSolver_HLLD
// Description :  Approximate Riemann solver of Harten, Lax, and van Leer.
//                The wave speed is estimated by the same formula in the HLLE solver
//
// Note        :  1. Input data should be conserved variables
//                2. Ref : (a) Riemann Solvers and Numerical Methods for Fluid Dynamics - A Practical Introduction
//                             ~ by Eleuterio F. Toro
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                         (c) Batten et al., SIAM J. Sci. Comput., 18, 1553 (1997)
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
void Hydro_RiemannSolver_HLLD( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres )
{

   const real MaxErr2         = SQR(MAX_ERROR);
   const real ZERO            = (real)0.0;
   const real ONE             = (real)1.0;
   const real _TWO            = (real)0.5;
   const real Gamma_m1        = Gamma - ONE;
   const bool NormPassive_No  = false;
   const bool JeansMinPres_No = false;
   const int  IdxBx           = MAG_OFFSET + 0;
   const int  IdxBy           = MAG_OFFSET + 1;
   const int  IdxBz           = MAG_OFFSET + 2;

   real Con_L[NCOMP_TOTAL_PLUS_MAG], Con_R[NCOMP_TOTAL_PLUS_MAG], Pri_L[NCOMP_TOTAL_PLUS_MAG], Pri_R[NCOMP_TOTAL_PLUS_MAG];

   for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
   {
      Con_L[v] = L_In[v];
      Con_R[v] = R_In[v];
   }

   Hydro_Rotate3D( Con_L, XYZ, true, IdxBx );
   Hydro_Rotate3D( Con_R, XYZ, true, IdxBx );

   const real BxL = Con_L[IdxBx];
   const real ByL = Con_L[IdxBy];
   const real BzL = Con_L[IdxBz];
   const real BxR = Con_R[IdxBx];
   const real ByR = Con_R[IdxBy];
   const real BzR = Con_R[IdxBz];

#  ifdef GAMER_DEBUG
   if ( BxL != BxR )
      printf( "ERROR : BxL (%24.17e) != BxR (%24.17e) for XYZ %d at file <%s>, line <%d>, function <%s>!!\n",
              BxL, BxR, XYZ, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   Hydro_Con2Pri( Con_L, Pri_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
   Hydro_Con2Pri( Con_R, Pri_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

   const real Vx_min = FMIN( Pri_L[1] , Pri_R[1] );
   const real Vx_max = FMAX( Pri_L[1] , Pri_R[1] );

   real tmp_1, tmp_2, crit, crit_Bx;
   real _RhoL, _RhoR;
   real sqrt_RhoLst, sqrt_RhoRst;
   real PT_L, PT_R, PT_st;
   real Bx, _Bx, Bx2, _Bx2, BtL2, BtR2, B2L_d2, B2R_d2;
   real a2 , Cf2 ,Cax2, Cat2, Ca2_plus_a2, Ca2_min_a2, Cf2_min_Cs2;
   real Cf_L, Cf_R, Cf_max;
   real Sd_L, Sd_R, Sdm_L, Sdm_R, SdL_SdmL, SdR_SdmR;
   real VBdot_Lst, VBdot_Rst;

   real Speed[5] = { ZERO, ZERO, ZERO, ZERO, ZERO };
   real Con_Lst[NCOMP_TOTAL_PLUS_MAG], Con_Ldst[NCOMP_TOTAL_PLUS_MAG];
   real Con_Rdst[NCOMP_TOTAL_PLUS_MAG], Con_Rst[NCOMP_TOTAL_PLUS_MAG];
   real Pri_Lst[NCOMP_TOTAL_PLUS_MAG], Pri_Rst[NCOMP_TOTAL_PLUS_MAG];
   real Flux_L[NCOMP_TOTAL_PLUS_MAG], Flux_R[NCOMP_TOTAL_PLUS_MAG];

   _RhoL       = ONE/Con_L[0];
   _RhoR       = ONE/Con_R[0];
   Bx          = BxL;
   _Bx         = ONE/Bx;
   Bx2         = SQR( Bx );
   _Bx2        = ONE/Bx2;
   BtL2        = SQR( ByL ) + SQR( BzL );
   BtR2        = SQR( ByR ) + SQR( BzR );
   B2L_d2      = _TWO*( Bx2 + BtL2 );
   B2R_d2      = _TWO*( Bx2 + BtR2 );
   PT_L        = Pri_L[4] + B2L_d2;
   PT_R        = Pri_R[4] + B2R_d2;

   a2          = Gamma*Pri_L[4]*_RhoL;
   Cax2        = Bx2*_RhoL;
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
   }

   Cf_L = SQRT( Cf2 );

   a2          = Gamma*Pri_R[4]*_RhoR;
   Cax2        = Bx2*_RhoR;
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
   }

   Cf_R = SQRT( Cf2 );

   Cf_max   = FMAX( Cf_L , Cf_R );
   Speed[0] = Vx_min - Cf_max;
   Speed[4] = Vx_max + Cf_max;

   Hydro_Con2Flux( 0, Flux_L, Con_L, Gamma_m1, MinPres );
   Hydro_Con2Flux( 0, Flux_R, Con_R, Gamma_m1, MinPres );


// return the upwind fluxes if flow is supersonic
   if ( Speed[0] >= ZERO )
   {
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   Flux_Out[v] = Flux_L[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false, IdxBx );

      return;
   }

   if ( Speed[4] <= ZERO )
   {
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   Flux_Out[v] = Flux_R[v];

      Hydro_Rotate3D( Flux_Out, XYZ, false, IdxBx );

      return;
   }

   Sd_L           = Speed[0] - Pri_L[1];
   Sd_R           = Speed[4] - Pri_R[1];
   tmp_1          = Sd_L*Pri_L[0];
   tmp_2          = Sd_R*Pri_R[0];
   Speed[2]       = ( tmp_2*Pri_R[1] - tmp_1*Pri_L[1] - PT_R + PT_L ) / ( tmp_2 - tmp_1);

   Sdm_L          = Speed[0] - Speed[2];
   Sdm_R          = Speed[4] - Speed[2];
   SdL_SdmL       = Sd_L / Sdm_L;
   SdR_SdmR       = Sd_R / Sdm_R;
   Con_Lst[0]     = Con_L[0]*SdL_SdmL;
   Con_Rst[0]     = Con_R[0]*SdR_SdmR;
   sqrt_RhoLst    = SQRT( Con_Lst[0] );
   sqrt_RhoRst    = SQRT( Con_Rst[0] );

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( Hydro_CheckNegative(Con_Lst[0]) )
      printf( "ERROR : invalid Con_Lst[0] (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Con_Lst[0], __FILE__, __LINE__, __FUNCTION__ );

   if ( Hydro_CheckNegative(Con_Rst[0]) )
      printf( "ERROR : invalid Con_Rst[0] (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Con_Rst[0], __FILE__, __LINE__, __FUNCTION__ );
#  endif

   tmp_1          = FABS( Bx );
   Speed[1]       = Speed[2] - tmp_1/sqrt_RhoLst;
   Speed[3]       = Speed[2] + tmp_1/sqrt_RhoRst;
   PT_st          = PT_L + Pri_L[0]*Sd_L*( Sd_L - Sdm_L );

   Con_Lst[    1] = Con_Lst[0]*Speed[2];
   Con_Lst[IdxBx] = Bx;

   tmp_1          = FMIN( B2L_d2 , B2R_d2 );

   if ( tmp_1 == ZERO )       crit_Bx = ZERO;
   else                       crit_Bx = FABS( _TWO*Bx2/tmp_1 );

   if ( crit_Bx < MaxErr2 )   crit = Con_L[0]*Sd_L*Sdm_L;
   else                       crit = Con_L[0]*Sd_L*Sdm_L*_Bx2 - ONE;

   if ( FABS( crit ) < MAX_ERROR )
   {
      Con_Lst[2] = Con_Lst[0]*Pri_L[2];
      Con_Lst[3] = Con_Lst[0]*Pri_L[3];

      Con_Lst[IdxBy] = ByL;
      Con_Lst[IdxBz] = BzL;
   }

   else
   {
      if ( crit_Bx < MaxErr2 )
      {
        Con_Lst[    2] = Con_Lst[0]*Pri_L[2];
        Con_Lst[    3] = Con_Lst[0]*Pri_L[3];
        Con_Lst[IdxBy] = ByL*SdL_SdmL;
        Con_Lst[IdxBz] = BzL*SdL_SdmL;
      }

      else
      {
        tmp_1          = ONE/crit;
        tmp_2          = ( Sd_L - Sdm_L )*_Bx*tmp_1;
        Con_Lst[    2] = Con_Lst[0]*( Pri_L[2] - Pri_L[IdxBy]*tmp_2 );
        Con_Lst[    3] = Con_Lst[0]*( Pri_L[3] - Pri_L[IdxBz]*tmp_2 );
        tmp_2          = ( Con_L[0]*Sd_L*Sd_L - Bx2 )*_Bx2*tmp_1;
        Con_Lst[IdxBy] = ByL*tmp_2;
        Con_Lst[IdxBz] = BzL*tmp_2;
      }
   } // if if ( FABS( crit ) < MAX_ERROR ) ... else ...

   VBdot_Lst      = ( Con_Lst[1]*Con_Lst[IdxBx] + Con_Lst[2]*Con_Lst[IdxBy] + Con_Lst[3]*Con_Lst[IdxBz] ) / Con_Lst[0];
   Con_Lst[    4] = (  Sd_L*Con_L[4] - PT_L*Pri_L[1] + PT_st*Speed[2] +
                       Bx*( Pri_L[1]*Pri_L[IdxBx] + Pri_L[2]*Pri_L[IdxBy] + Pri_L[3]*Pri_L[IdxBz] - VBdot_Lst )  ) / Sdm_L;

   Hydro_Con2Pri( Con_Lst, Pri_Lst, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

   Con_Rst[    1] = Con_Rst[0]*Speed[2];
   Con_Rst[IdxBx] = Bx;

   if ( crit_Bx < MaxErr2 )   crit = Con_R[0]*Sd_R*Sdm_R;
   else                       crit = Con_R[0]*Sd_R*Sdm_R*_Bx2 - ONE;

   if ( FABS( crit ) <  MAX_ERROR )
   {
      Con_Rst[    2] = Con_Rst[0]*Pri_R[2];
      Con_Rst[    3] = Con_Rst[0]*Pri_R[3];
      Con_Rst[IdxBy] = ByR;
      Con_Rst[IdxBz] = BzR;
   }

   else
   {
      if ( crit_Bx < MaxErr2 )
      {
        Con_Rst[    2] = Con_Rst[0]*Pri_R[2];
        Con_Rst[    3] = Con_Rst[0]*Pri_R[3];
        Con_Rst[IdxBy] = ByR*SdR_SdmR;
        Con_Rst[IdxBz] = BzR*SdR_SdmR;
      }

      else
      {
        tmp_1          = ONE/crit;
        tmp_2          = ( Sd_R - Sdm_R )*_Bx*tmp_1;
        Con_Rst[    2] = Con_Rst[0]*( Pri_R[2] - Pri_R[IdxBy]*tmp_2 );
        Con_Rst[    3] = Con_Rst[0]*( Pri_R[3] - Pri_R[IdxBz]*tmp_2 );
        tmp_2          = ( Con_R[0]*Sd_R*Sd_R - Bx2 )*_Bx2*tmp_1;
        Con_Rst[IdxBy] = ByR*tmp_2;
        Con_Rst[IdxBz] = BzR*tmp_2;
      }
   } // if ( FABS( crit ) <  MAX_ERROR ) ... else ...

   VBdot_Rst  = ( Con_Rst[1]*Con_Rst[IdxBx] + Con_Rst[2]*Con_Rst[IdxBy] + Con_Rst[3]*Con_Rst[IdxBz] ) / Con_Rst[0];
   Con_Rst[4] = (  Sd_R*Con_R[4] - PT_R*Pri_R[1] + PT_st*Speed[2] +
                   Bx*( Pri_R[1]*Pri_R[IdxBx] + Pri_R[2]*Pri_R[IdxBy] + Pri_R[3]*Pri_R[IdxBz] - VBdot_Rst )  ) / Sdm_R;

   Hydro_Con2Pri( Con_Rst, Pri_Rst, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

   if ( crit_Bx < MaxErr2 )
   {
     for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
     {
       Con_Ldst[v] = Con_Lst[v];
       Con_Rdst[v] = Con_Rst[v];
     }
   }

   else
   {
     const real invsumd = ONE/( sqrt_RhoLst + sqrt_RhoRst );
     const real Bxsig   = SIGN( Bx );

     Con_Ldst[0] = Con_Lst[0];
     Con_Rdst[0] = Con_Rst[0];

     Con_Ldst[1] = Con_Lst[1];
     Con_Rdst[1] = Con_Rst[1];

     tmp_1 = invsumd*(  sqrt_RhoLst*Pri_Lst[2] + sqrt_RhoRst*Pri_Rst[2] + Bxsig*( Pri_Rst[IdxBy] - Pri_Lst[IdxBy] )  );
     Con_Ldst[2] = Con_Ldst[0]*tmp_1;
     Con_Rdst[2] = Con_Rdst[0]*tmp_1;

     tmp_1 = invsumd*(  sqrt_RhoLst*Pri_Lst[3] + sqrt_RhoRst*Pri_Rst[3] + Bxsig*( Pri_Rst[IdxBz] - Pri_Lst[IdxBz] )  );
     Con_Ldst[3] = Con_Ldst[0]*tmp_1;
     Con_Rdst[3] = Con_Rdst[0]*tmp_1;

     tmp_1 = invsumd*(  sqrt_RhoLst*Pri_Rst[IdxBy] + sqrt_RhoRst*Pri_Lst[IdxBy] +
                        Bxsig*sqrt_RhoLst*sqrt_RhoRst*( Pri_Rst[2] - Pri_Lst[2] )  );
     Con_Ldst[IdxBy] = tmp_1;
     Con_Rdst[IdxBy] = tmp_1;

     tmp_1 = invsumd*(  sqrt_RhoLst*Pri_Rst[IdxBz] + sqrt_RhoRst*Pri_Lst[IdxBz] +
                        Bxsig*sqrt_RhoLst*sqrt_RhoRst*( Pri_Rst[3] - Pri_Lst[3] )  );
     Con_Ldst[IdxBz] = tmp_1;
     Con_Rdst[IdxBz] = tmp_1;

     tmp_1 = Speed[2]*Bx + ( Con_Ldst[2]*Con_Ldst[IdxBy] + Con_Ldst[3]*Con_Ldst[IdxBz] ) / Con_Ldst[0];
     Con_Ldst[4] = Con_Lst[4] - sqrt_RhoLst*Bxsig*( VBdot_Lst - tmp_1 );
     Con_Rdst[4] = Con_Rst[4] + sqrt_RhoRst*Bxsig*( VBdot_Rst - tmp_1 );
   } // if ( crit_Bx < MaxErr2 ) ... else ...


// evaluate the HLLD fluxes
   if ( Speed[1] >= ZERO )
   {
      for (int v=0; v<NCOMP_FLUID; v++)
      Flux_Out[    v] = Flux_L[v] + Speed[0]*( Con_Lst[v] - Con_L[v] );

      Flux_Out[IdxBx] = ZERO;
      Flux_Out[IdxBy] = Flux_L[IdxBy] + Speed[0]*(Con_Lst[IdxBy] - ByL);
      Flux_Out[IdxBz] = Flux_L[IdxBz] + Speed[0]*(Con_Lst[IdxBz] - BzL);
   }

   else if ( Speed[2] >= ZERO )
   {
      tmp_1 = Speed[1] - Speed[0];

      for (int v=0; v<NCOMP_FLUID; v++)
      Flux_Out[    v] = Flux_L[v] - Speed[0]*Con_L[v] - tmp_1*Con_Lst[v] + Speed[1]*Con_Ldst[v];

      Flux_Out[IdxBx] = ZERO;
      Flux_Out[IdxBy] = Flux_L[IdxBy] - Speed[0]*ByL - tmp_1*Con_Lst[IdxBy] + Speed[1]*Con_Ldst[IdxBy];
      Flux_Out[IdxBz] = Flux_L[IdxBz] - Speed[0]*BzL - tmp_1*Con_Lst[IdxBz] + Speed[1]*Con_Ldst[IdxBz];
   }

   else if ( Speed[3] > ZERO )
   {
      tmp_1 = Speed[3] - Speed[4];

      for (int v=0; v<NCOMP_FLUID; v++)
      Flux_Out[    v] = Flux_R[v] - Speed[4]*Con_R[v] - tmp_1*Con_Rst[v] + Speed[3]*Con_Rdst[v];

      Flux_Out[IdxBx] = ZERO;
      Flux_Out[IdxBy] = Flux_R[IdxBy] - Speed[4]*ByR - tmp_1*Con_Rst[IdxBy] + Speed[3]*Con_Rdst[IdxBy];
      Flux_Out[IdxBz] = Flux_R[IdxBz] - Speed[4]*BzR - tmp_1*Con_Rst[IdxBz] + Speed[3]*Con_Rdst[IdxBz];
   }

   else
   {
      for (int v=0; v<NCOMP_FLUID; v++)
      Flux_Out[    v] = Flux_R[v] + Speed[4]*( Con_Rst[v] - Con_R[v] );

      Flux_Out[IdxBx] = ZERO;
      Flux_Out[IdxBy] = Flux_R[IdxBy] + Speed[4]*( Con_Rst[IdxBy] - ByR );
      Flux_Out[IdxBz] = Flux_R[IdxBz] + Speed[4]*( Con_Rst[IdxBz] - BzR );
   } // if ( Speed[x] > ZERO ) ... else ...


// evaluate the fluxes for passive scalars
#  if ( NCOMP_PASSIVE > 0 )
   if ( Flux_Out[FLUX_DENS] >= ZERO )
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoL;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = Con_L[v]*vx;
   }

   else
   {
      const real vx = Flux_Out[FLUX_DENS]*_RhoR;

      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  Flux_Out[v] = Con_R[v]*vx;
   }
#  endif


// restore the correct order
   Hydro_Rotate3D( Flux_Out, XYZ, false, IdxBx );

} // FUNCTION : Hydro_RiemannSolver_HLLD



#endif // #if ( MODEL == HYDRO  &&  defined MHD )



#endif // #ifndef __CUFLU_RIEMANNSOLVER_HLLD__
