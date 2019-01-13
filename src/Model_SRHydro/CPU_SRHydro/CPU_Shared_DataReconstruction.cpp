#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP) )

static void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                        const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                        const real Gamma, const int XYZ, real LimitedSlope[], int iteration );

bool CPU_CheckUnphysical( const real Con[], const real Pri[], const char s[], const int line, bool show);
static bool boolean;
#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input and output data should be primitive variables
//                3. The PLM and PPM data reconstruction functions share the same function name
//                4. The face-centered variables will be advanced by half time-step for the CTU scheme
//                5. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   : [ 1] PriVar         : Array storing the input primitive variables
//               [ 2] FC_Var         : Array to store the output face-centered primitive variables
//               [ 3] NIn            : Size of the input array "PriVar" in one direction
//               [ 4] NGhost         : Size of the ghost zone
//                                     --> "NIn-2*NGhost" cells will be computed along each direction
//                                     --> The size of the output array "FC_Var" is assumed to be "(NIn-2*NGhost)^3"
//                                     --> The reconstructed data at cell (i,j,k) will be stored in the
//                                         array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//               [ 5] Gamma          : Ratio of specific heats
//               [ 6] LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//               [ 7] MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//               [ 8] EP_Coeff       : Coefficient of the extrema-preserving limiter
//               [ 9] dt             : Time interval to advance solution (for the CTU scheme)
//               [10] dh             : Grid size (for the CTU scheme)
//            [11/12] MinDens/Pres   : Minimum allowed density and pressure
//               [13] iteration      :
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                             const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres, int iteration )
{
   const int dr1[3] = { 1, NIn, NIn*NIn };
   const int NOut   = NIn - 2*NGhost;                    // number of output grids
   int  ID1, ID2, ID1_L, ID1_R, ID1_LL, ID1_RR, dL, dR;
   real Min, Max;
   real LimitedSlope[NCOMP_TOTAL] = { (real)0.0 };


   for (int k1=NGhost, k2=0;  k1<NGhost+NOut;  k1++, k2++)
   for (int j1=NGhost, j2=0;  j1<NGhost+NOut;  j1++, j2++)
   for (int i1=NGhost, i2=0;  i1<NGhost+NOut;  i1++, i2++)
   {
      ID1 = (k1*NIn  + j1)*NIn  + i1;
      ID2 = (k2*NOut + j2)*NOut + i2;

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {

//       (2-1) evaluate the monotonic slope
         dL    = 2*d;
         dR    = dL+1;
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];

         if ( LR_Limiter == EXTPRE )
         {
            ID1_LL = ID1 - 2*dr1[d];
            ID1_RR = ID1 + 2*dr1[d];

            LimitSlope( PriVar[ID1_LL], PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], PriVar[ID1_RR], LR_Limiter,
                        MinMod_Coeff, EP_Coeff, Gamma, d, LimitedSlope, iteration );
         }

         else
         {
            LimitSlope( NULL, PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], NULL, LR_Limiter,
                        MinMod_Coeff, NULL_REAL, Gamma, d, LimitedSlope, iteration );
         }


//       (2-2) get the face-centered primitive variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            FC_Var[ID2][dL][v] = PriVar[ID1][v] - (real)0.5*LimitedSlope[v];
            FC_Var[ID2][dR][v] = PriVar[ID1][v] + (real)0.5*LimitedSlope[v];
         }


//       (2-3) ensure the face-centered variables lie between neighboring cell-centered values
         if ( LR_Limiter != EXTPRE )
         {
            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               Min = ( PriVar[ID1][v] < PriVar[ID1_L][v] ) ? PriVar[ID1][v] : PriVar[ID1_L][v];
               Max = ( PriVar[ID1][v] > PriVar[ID1_L][v] ) ? PriVar[ID1][v] : PriVar[ID1_L][v];
               FC_Var[ID2][dL][v] = ( FC_Var[ID2][dL][v] > Min  ) ? FC_Var[ID2][dL][v] : Min;
               FC_Var[ID2][dL][v] = ( FC_Var[ID2][dL][v] < Max  ) ? FC_Var[ID2][dL][v] : Max;
               FC_Var[ID2][dR][v] = (real)2.0*PriVar[ID1][v] - FC_Var[ID2][dL][v];

               Min = ( PriVar[ID1][v] < PriVar[ID1_R][v] ) ? PriVar[ID1][v] : PriVar[ID1_R][v];
               Max = ( PriVar[ID1][v] > PriVar[ID1_R][v] ) ? PriVar[ID1][v] : PriVar[ID1_R][v];
               FC_Var[ID2][dR][v] = ( FC_Var[ID2][dR][v] > Min  ) ? FC_Var[ID2][dR][v] : Min;
               FC_Var[ID2][dR][v] = ( FC_Var[ID2][dR][v] < Max  ) ? FC_Var[ID2][dR][v] : Max;
               FC_Var[ID2][dL][v] = (real)2.0*PriVar[ID1][v] - FC_Var[ID2][dR][v];
            }
#          if ( defined ( CHECK_NEGATIVE_IN_FLUID ) &&  EXTRAPOLATE != CONSERVED_QUANTITIES )
           boolean = CPU_CheckUnphysical(NULL,FC_Var[ID2][dR], __FUNCTION__, __LINE__, true);
           boolean = CPU_CheckUnphysical(NULL,FC_Var[ID2][dL], __FUNCTION__, __LINE__, true);
#          elif ( defined ( CHECK_NEGATIVE_IN_FLUID ) &&  EXTRAPOLATE == CONSERVED_QUANTITIES )
           boolean = CPU_CheckUnphysical(FC_Var[ID2][dR], NULL, __FUNCTION__, __LINE__, true);
           boolean = CPU_CheckUnphysical(FC_Var[ID2][dL], NULL, __FUNCTION__, __LINE__, true);
#          endif
         }
         else // for the extrema-preserving limiter --> ensure positive density and pressure
         {
/*
            FC_Var[ID2][dL][0] = CPU_CheckMinDens( FC_Var[ID2][dL][0], MinDens );
            FC_Var[ID2][dR][0] = CPU_CheckMinDens( FC_Var[ID2][dR][0], MinDens );

            FC_Var[ID2][dL][4] = CPU_CheckMinPres( FC_Var[ID2][dL][4], MinPres );
            FC_Var[ID2][dR][4] = CPU_CheckMinPres( FC_Var[ID2][dR][4], MinPres );
*/
#          if ( defined ( CHECK_NEGATIVE_IN_FLUID ) &&  EXTRAPOLATE != CONSERVED_QUANTITIES )
           boolean = CPU_CheckUnphysical(NULL,FC_Var[ID2][dR], __FUNCTION__, __LINE__, true);
           boolean = CPU_CheckUnphysical(NULL,FC_Var[ID2][dL], __FUNCTION__, __LINE__, true);
#          elif ( defined ( CHECK_NEGATIVE_IN_FLUID ) && ( EXTRAPOLATE == CONSERVED_QUANTITIES ) )
           boolean = CPU_CheckUnphysical(FC_Var[ID2][dR], NULL, __FUNCTION__, __LINE__, true);
           boolean = CPU_CheckUnphysical(FC_Var[ID2][dL], NULL, __FUNCTION__, __LINE__, true);
#          endif
         }


      } // for (int d=0; d<3; d++)
   } // k,j,i

} // FUNCTION : CPU_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input and output data should be primitive variables
//                3. The PLM and PPM data reconstruction functions share the same function name
//                4. The face-centered variables will be advanced by half time-step for the CTU scheme
//                5. Currently the extrema-preserving limiter is not supported in PPM
//                6. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                7. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  PriVar         : Array storing the input primitive variables
//                FC_Var         : Array to store the output face-centered primitive variables
//                NIn            : Size of the input array "PriVar" in one direction
//                NGhost         : Size of the ghost zone
//                                  --> "NIn-2*NGhost" cells will be computed along each direction
//                                  --> The size of the output array "FC_Var" is assumed to be "(NIn-2*NGhost)^3"
//                                  --> The reconstructed data at cell (i,j,k) will be stored in the
//                                      array "FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma          : Ratio of specific heats
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter (useless in PPM)
//                dt             : Time interval to advance solution (for the CTU scheme)
//                dh             : Grid size (for the CTU scheme)
//                MinDens/Pres   : Minimum allowed density and pressure
//------------------------------------------------------------------------------------------------------
void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                             const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                             const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres, int iteration  )
{

// check
#  ifdef GAMER_DEBUG
   if ( LR_Limiter == EXTPRE )
      Aux_Error( ERROR_INFO, "PPM reconstruction does NOT support the extrema-preserving limiter !!\n");
#  endif


   const int NOut   = NIn - 2*NGhost;                    // number of output grids
   const int NSlope = NOut + 2;                          // number of grids required to store the slope data
   const int dr1[3] = { 1, NIn, NIn*NIn };
   const int dr3[3] = { 1, NSlope, NSlope*NSlope };

   int ID1, ID2, ID3, ID1_L, ID1_R, ID3_L, ID3_R, dL, dR;
   real LimitedSlope[NCOMP_TOTAL] = { (real)0.0 };
   real CC_L, CC_R, CC_C, dCC_L, dCC_R, dCC_C, FC_L, FC_R, dFC[NCOMP_TOTAL], dFC6[NCOMP_TOTAL], Max, Min;

   real (*Slope_PPM)[3][NCOMP_TOTAL] = new real [ NSlope*NSlope*NSlope ][3][NCOMP_TOTAL];



// (2-1) evaluate the monotonic slope
   for (int k1=NGhost-1, k2=0;  k1<NGhost-1+NSlope;  k1++, k2++)
   for (int j1=NGhost-1, j2=0;  j1<NGhost-1+NSlope;  j1++, j2++)
   for (int i1=NGhost-1, i2=0;  i1<NGhost-1+NSlope;  i1++, i2++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NSlope + j2)*NSlope + i2;

//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];

         if ( LR_Limiter == EXTPRE )
         {
            Aux_Error( ERROR_INFO, "PPM reconstruction does NOT support the extrema-preserving limiter !!\n");

            /*
            ID1_LL = ID1 - 2*dr1[d];
            ID1_RR = ID1 + 2*dr1[d];

            LimitSlope( PriVar[ID1_LL], PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], PriVar[ID1_RR], LR_Limiter,
                        MinMod_Coeff, EP_Coeff, Gamma, d, LimitedSlope, iteration );
            */
         }

         else
         {
            LimitSlope( NULL, PriVar[ID1_L], PriVar[ID1], PriVar[ID1_R], NULL, LR_Limiter,
                        MinMod_Coeff, NULL_REAL, Gamma, d, LimitedSlope, iteration );
         }


//       store the slope to the array "Slope_PPM"
         for (int v=0; v<NCOMP_TOTAL; v++)   Slope_PPM[ID2][d][v] = LimitedSlope[v];

      } // for (int d=0; d<3; d++)
   } // k,j,i


   for (int k1=NGhost, k2=0, k3=1;  k1<NGhost+NOut;  k1++, k2++, k3++)
   for (int j1=NGhost, j2=0, j3=1;  j1<NGhost+NOut;  j1++, j2++, j3++)
   for (int i1=NGhost, i2=0, i3=1;  i1<NGhost+NOut;  i1++, i2++, i3++)
   {
      ID1 = (k1*NIn    + j1)*NIn    + i1;
      ID2 = (k2*NOut   + j2)*NOut   + i2;
      ID3 = (k3*NSlope + j3)*NSlope + i3;


//    (2-3) get the face-centered primitive variables
//    loop over different spatial directions
      for (int d=0; d<3; d++)
      {
         dL    = 2*d;
         dR    = dL+1;
         ID1_L = ID1 - dr1[d];
         ID1_R = ID1 + dr1[d];
         ID3_L = ID3 - dr3[d];
         ID3_R = ID3 + dr3[d];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
//          (2-3-1) parabolic interpolation
            CC_L  = PriVar[ID1_L][v];
            CC_R  = PriVar[ID1_R][v];
            CC_C  = PriVar[ID1  ][v];

            dCC_L = Slope_PPM[ID3_L][d][v];
            dCC_R = Slope_PPM[ID3_R][d][v];
            dCC_C = Slope_PPM[ID3  ][d][v];

            FC_L  = (real)0.5*( CC_C + CC_L ) - (real)1.0/(real)6.0*( dCC_C - dCC_L );
            FC_R  = (real)0.5*( CC_C + CC_R ) - (real)1.0/(real)6.0*( dCC_R - dCC_C );


//          (2-3-2) monotonicity constraint
            dFC [v] = FC_R - FC_L;
            dFC6[v] = (real)6.0*(  CC_C - (real)0.5*( FC_L + FC_R )  );

            if (  ( FC_R - CC_C )*( CC_C - FC_L ) <= (real)0.0  )
            {
               FC_L = CC_C;
               FC_R = CC_C;
            }
            else if ( dFC[v]*dFC6[v] > +dFC[v]*dFC[v] )
               FC_L = (real)3.0*CC_C - (real)2.0*FC_R;
            else if ( dFC[v]*dFC6[v] < -dFC[v]*dFC[v] )
               FC_R = (real)3.0*CC_C - (real)2.0*FC_L;


//          (2-3-3) ensure the face-centered variables lie between neighboring cell-centered values
            Min  = ( CC_C < CC_L ) ? CC_C : CC_L;
            Max  = ( CC_C > CC_L ) ? CC_C : CC_L;
            FC_L = ( FC_L > Min  ) ? FC_L : Min;
            FC_L = ( FC_L < Max  ) ? FC_L : Max;

            Min  = ( CC_C < CC_R ) ? CC_C : CC_R;
            Max  = ( CC_C > CC_R ) ? CC_C : CC_R;
            FC_R = ( FC_R > Min  ) ? FC_R : Min;
            FC_R = ( FC_R < Max  ) ? FC_R : Max;


            FC_Var[ID2][dL][v] = FC_L;
            FC_Var[ID2][dR][v] = FC_R;

         } // for (int v=0; v<NCOMP_TOTAL; v++)

      } // for (int d=0; d<3; d++)
   } // k,j,i

   delete [] Slope_PPM;

} // FUNCTION : CPU_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )




//-------------------------------------------------------------------------------------------------------
// Function    :  LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//                2. The L2 and R2 elements are useful only for the extrema-preserving limiter
//
// Parameter   :  L2             : Element x-2
//                L1             : Element x-1
//                C0             : Element x
//                R1             : Element x+1
//                R2             : Element x+2
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                Gamma          : Ratio of specific heats
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ            : Target spatial direction : (0/1/2) --> (x/y/z)
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                LimitedSlope  : Array to store the output monotonic slope
//-------------------------------------------------------------------------------------------------------
void LimitSlope( const real L2[], const real L1[], const real C0[], const real R1[], const real R2[],
                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                 const real Gamma, const int XYZ, real LimitedSlope[], int iteration )
{

// check
#  ifdef GAMER_DEBUG
   if ( LR_Limiter == EXTPRE  &&  ( L2 == NULL || R2 == NULL )  )
      Aux_Error( ERROR_INFO, "input element == NULL !!\n" );
#  endif


   real Slope_L[NCOMP_TOTAL], Slope_R[NCOMP_TOTAL], Slope_C[NCOMP_TOTAL], Slope_A[NCOMP_TOTAL];
   real Slope_LL[NCOMP_TOTAL], Slope_RR[NCOMP_TOTAL], Slope_LR;
   real D2_L, D2_R, D2_C, D2_Sign, D2_Limiter, Slope_Sign;  // variables for the extrema-preserving limiter
   real beta_L, beta_R, Xi_L, Xi_R, delta[NCOMP_TOTAL];
   int v_min, v_max;

//   if( iteration ) 
//    {
//      v_min = MOMX;
//      v_max = ENGY;
//    }
//     else
//    {
      v_min = 0;
      v_max = NCOMP_TOTAL;
//    }

   for (int v=v_min; v<v_max; v++)
   {
      Slope_L[v] = C0[v] - L1[v];
      Slope_R[v] = R1[v] - C0[v];

//    evaluate different slopes
      switch ( LR_Limiter )
      {
//       generalized MinMod
         case GMINMOD: 
         Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );

           if (  Slope_L[v]*Slope_R[v] > (real)0.0  )
           {
             Slope_L[v] *= MinMod_Coeff;
             Slope_R[v] *= MinMod_Coeff;
             LimitedSlope[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
             LimitedSlope[v]  = FMIN(  FABS( Slope_C[v] ), LimitedSlope[v]  );
             LimitedSlope[v] *= SIGN( Slope_C[v] );
           } else LimitedSlope[v] = (real)0.0;
         break;

//       van-Leer + generalized MinMod
         case VL_GMINMOD:           
         Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );

           if (  Slope_L[v]*Slope_R[v] > (real)0.0 )
              Slope_A[v] = (real)2.0*Slope_L[v]*Slope_R[v]/( Slope_L[v] + Slope_R[v] );
           else
              Slope_A[v] = (real)0.0;

           if (  Slope_L[v]*Slope_R[v] > (real)0.0 )
           {
             Slope_L[v] *= MinMod_Coeff;
             Slope_R[v] *= MinMod_Coeff;
             LimitedSlope[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
             LimitedSlope[v]  = FMIN(  FABS( Slope_C[v] ), LimitedSlope[v]  );
             LimitedSlope[v]  = FMIN(  FABS( Slope_A[v] ), LimitedSlope[v]  );
             LimitedSlope[v] *= SIGN( Slope_C[v] );
           } else LimitedSlope[v] = (real)0.0;
         break;

//       extrema-preserving
         case EXTPRE:
         Slope_C[v] = (real)0.5*( Slope_L[v] + Slope_R[v] );

           Slope_LL[v] = L1[v] - L2[v];
           Slope_RR[v] = R2[v] - R1[v];

           if (  Slope_L[v]*Slope_R[v] > (real)0.0  &&  Slope_LL[v]*Slope_RR[v] > (real)0.0  )
           {
             Slope_L[v] *= MinMod_Coeff;
             Slope_R[v] *= MinMod_Coeff;
             LimitedSlope[v]  = FMIN(  FABS( Slope_L[v] ), FABS( Slope_R[v] )  );
             LimitedSlope[v]  = FMIN(  FABS( Slope_C[v] ), LimitedSlope[v]  );
             LimitedSlope[v] *= SIGN( Slope_C[v] );
           }
           else
           {
             D2_L = Slope_L [v] - Slope_LL[v];
             D2_R = Slope_RR[v] - Slope_R [v];
             D2_C = Slope_R [v] - Slope_L [v];

             D2_Sign    = SIGN( D2_C );
             Slope_Sign = SIGN( Slope_C[v] );

             D2_Limiter = FMIN(  FABS(D2_C), FMIN( FMAX(D2_Sign*D2_L, (real)0.0), FMAX(D2_Sign*D2_R, (real)0.0) )  );

             if ( D2_Sign*Slope_Sign < (real)0.0 )
                   LimitedSlope[v] = FMIN( (real)1.5*EP_Coeff*D2_Limiter, MinMod_Coeff*FABS(Slope_L[v]) );
             else  LimitedSlope[v] = FMIN( (real)1.5*EP_Coeff*D2_Limiter, MinMod_Coeff*FABS(Slope_R[v]) );

             LimitedSlope[v] = Slope_Sign * FMIN( FABS(Slope_C[v]), LimitedSlope[v] );
           }
         break;

//       van-Leer, Ref: eq.(14.54) in Toro
         case VANLEER:              

           if (  Slope_L[v]*Slope_R[v] > (real)0.0 )
            {
              LimitedSlope[v] = (real)2.0*Slope_LR/( Slope_L[v] + Slope_R[v] );
              LimitedSlope[v] *= MinMod_Coeff;
            } else LimitedSlope[v] = 0.0;
            break;

//       van-Leer Albada 1, Ref: eq.(14.55) in Toro
         case ALBADA:               

           Slope_LR = Slope_L[v]*Slope_R[v];
           if (  Slope_LR > (real)0.0 )
            {
              LimitedSlope[v] = Slope_LR*( Slope_L[v] + Slope_R[v] ) /
              ( Slope_L[v]*Slope_L[v] + Slope_R[v]*Slope_R[v] );
             
              LimitedSlope[v] *= MinMod_Coeff;

            } else LimitedSlope[v] = 0.0;
            break;

//       minbee, Ref: eq.(14.44) in Toro
         case MINBEE:

           if ( Slope_R[v] > (real)0.0 )
           {
             LimitedSlope[v] = FMIN( Slope_L[v],      Slope_R[v] );
             LimitedSlope[v] = FMAX(          0, LimitedSlope[v] );
             LimitedSlope[v] *= MinMod_Coeff;
           }
           else
           {
             LimitedSlope[v] = FMAX( Slope_L[v],      Slope_R[v] );
             LimitedSlope[v] = FMIN(          0, LimitedSlope[v] );
             LimitedSlope[v] *= MinMod_Coeff;
           }
           break;
      
//       superbee, Ref: eq.(14.44) in Toro
         case SUPERBEE:

           real a, b;
           if ( Slope_R[v] > (real)0.0 )
           {
             a = FMIN(     Slope_L[v], 2.0*Slope_R[v] );
             b = FMIN( 2.0*Slope_L[v],     Slope_R[v] );
             LimitedSlope[v] = FMAX(   a,               b );
             LimitedSlope[v] = FMAX( 0.0, LimitedSlope[v] );
             LimitedSlope[v] *= MinMod_Coeff;
           }
           else
           {
             a = FMAX(     Slope_L[v], 2.0*Slope_R[v] );
             b = FMAX( 2.0*Slope_L[v],     Slope_R[v] );
             LimitedSlope[v] = FMIN(   a,               b );
             LimitedSlope[v] = FMIN( 0.0, LimitedSlope[v] );
             LimitedSlope[v] *= MinMod_Coeff;
           }
           break;

//       piece-wise constant
         case CONSTANT: 
           LimitedSlope[v] = 0.0;
           break;

         default :
            Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "LR_Limiter", LR_Limiter );
      }
   } // for (int v=0; v<NCOMP_TOTAL; v++)

} // FUNCTION : LimitSlope



#endif // #if ( !defined GPU  &&  MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP) )
