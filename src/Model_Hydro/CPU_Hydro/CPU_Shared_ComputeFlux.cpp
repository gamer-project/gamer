#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



#if   ( RSOLVER == EXACT )
extern void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                         const bool NormPassive, const int NNorm, const int NormIdx[],
                         const bool JeansMinPres, const real JeansMinPres_Coeff );
extern void CPU_RiemannSolver_Exact( const int XYZ, real eival_out[], real L_star_out[], real R_star_out[],
                                     real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( RSOLVER == ROE )
extern void CPU_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                   const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLE )
extern void CPU_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLC )
extern void CPU_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                    const real Gamma, const real MinPres );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver
//
// Note        :  1. Currently support the exact, HLLC, HLLE, and Roe solvers
//                2. The size of the input array "FC_Var" is assumed to be N_FC_VAR^3
//                   --> "N_FC_VAR-1" fluxes will be computed along each direction
//
// Parameter   :  FC_Var          : Array storing the input face-centered conserved variables
//                FC_Flux         : Array to store the output face-centered flux
//                NFlux           : Size of the array FC_Flux in each direction (must be >= N_FC_VAR-1)
//                                  --> The (i,j,k) flux will be stored in the array FC_Flux with
//                                      the index "(k*NFlux+j)*NFlux+i"
//                                  --> The (i,j,k) FC_Flux_x/y/z are defined at the +x/+y/+z surfaces of the
//                                      cell (i,j,k)
//                Gap             : Number of grids to be skipped in the transverse direction
//                                  --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed on each surface
//                Gamma           : Ratio of specific heats
//                CorrHalfVel     : true --> correcting the half-step velocity by gravity (for UNSPLIT_GRAVITY only)
//                Pot_USG         : Array storing the input potential for CorrHalfVel     (for UNSPLIT_GRAVITY only)
//                                  --> must have the same size as FC_Var ( (PS2+2)^3 )
//                Corner          : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                dt              : Time interval to advance the full-step solution       (for UNSPLIT_GRAVITY only)
//                dh              : Grid size                                             (for UNSPLIT_GRAVITY only)
//                Time            : Current physical time                                 (for UNSPLIT_GRAVITY only)
//                GravityType     : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//                MinPres         : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void CPU_ComputeFlux( const real FC_Var[][6][NCOMP_TOTAL], real FC_Flux[][3][NCOMP_TOTAL], const int NFlux, const int Gap,
                      const real Gamma, const bool CorrHalfVel, const real Pot_USG[], const double Corner[],
                      const real dt, const real dh, const double Time, const OptGravityType_t GravityType,
                      const double ExtAcc_AuxArray[], const real MinPres )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )
   {
      if (  ( GravityType == GRAVITY_SELF || GravityType == GRAVITY_BOTH )  &&  ( Pot_USG == NULL )  )
         Aux_Error( ERROR_INFO, "Pot_USG == NULL !!\n" );

      if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  &&  ( Corner == NULL )  )
         Aux_Error( ERROR_INFO, "Corner == NULL !!\n" );

      if ( N_FC_VAR != PS2+2 )
         Aux_Error( ERROR_INFO, "N_FC_VAR (%d) != PS2+2 (%d) !!\n", N_FC_VAR, PS2+2 );
   }
#  else
   if ( CorrHalfVel )   Aux_Error( ERROR_INFO, "CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG


   const int dID2[3] = { 1, N_FC_VAR, N_FC_VAR*N_FC_VAR };

   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL];
   int  ID1, ID2, dL, dR, start2[3]={0}, end1[3]={0};

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif

#  ifdef UNSPLIT_GRAVITY
   const real   GraConst = -(real)0.5*dt/dh;
   const int    dID3[3]  = { 1, USG_NXT_F, USG_NXT_F*USG_NXT_F };
   const int    didx     = USG_GHOST_SIZE - 1;    // assuming FC_Var has one ghost zone on each side
   const double dh_half  = 0.5*(double)dh;
   const real   dt_half  = (real)0.5*dt;

   int    d1, d2, d3, ID3;
   real   eL, eR, Acc[3];
   double xyz[3], CrShift[3];

// assuming N_FC_VAR = PS2+2 (ghostzone = 1 on each side) --> CrShift is the central coordinates of FC_Var[0]
   if (  CorrHalfVel  &&  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  )
      for (int d=0; d<3; d++)    CrShift[d] = Corner[d] - (double)dh;
#  endif


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      dL = 2*d;
      dR = dL+1;

#     ifdef UNSPLIT_GRAVITY
      if ( CorrHalfVel )
      {
         d1 =  d;
         d2 = (d+1)%3;
         d3 = (d+2)%3;
      }
#     endif

      switch ( d )
      {
         case 0 : start2[0] = 0;                start2[1] = Gap;              start2[2] = Gap;
                  end1  [0] = N_FC_VAR-1;       end1  [1] = N_FC_VAR-2*Gap;   end1  [2] = N_FC_VAR-2*Gap;   break;

         case 1 : start2[0] = Gap;              start2[1] = 0;                start2[2] = Gap;
                  end1  [0] = N_FC_VAR-2*Gap;   end1  [1] = N_FC_VAR-1;       end1  [2] = N_FC_VAR-2*Gap;   break;

         case 2 : start2[0] = Gap;              start2[1] = Gap;              start2[2] = 0;
                  end1  [0] = N_FC_VAR-2*Gap;   end1  [1] = N_FC_VAR-2*Gap;   end1  [2] = N_FC_VAR-1;       break;
      }

      for (int k1=0, k2=start2[2];  k1<end1[2];  k1++, k2++)
      for (int j1=0, j2=start2[1];  j1<end1[1];  j1++, j2++)
      for (int i1=0, i2=start2[0];  i1<end1[0];  i1++, i2++)
      {
         ID1 = (k1*NFlux    + j1)*NFlux    + i1;
         ID2 = (k2*N_FC_VAR + j2)*N_FC_VAR + i2;

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = FC_Var[ ID2         ][dR][v];
            ConVar_R[v] = FC_Var[ ID2+dID2[d] ][dL][v];
         }


//       correct the half-step velocity by gravity
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            Acc[0] = (real)0.0;
            Acc[1] = (real)0.0;
            Acc[2] = (real)0.0;

//          external gravity
            if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
            {
               xyz[0]  = CrShift[0] + (double)(i2*dh);
               xyz[1]  = CrShift[1] + (double)(j2*dh);
               xyz[2]  = CrShift[2] + (double)(k2*dh);
               xyz[d] += dh_half;

               CPU_ExternalAcc( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );

               for (int d=0; d<3; d++)    Acc[d] *= dt_half;
            }

//          self-gravity
            if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
            {
               ID3      = ( (k2+didx)*USG_NXT_F + (j2+didx) )*USG_NXT_F + (i2+didx);

               Acc[d1] +=            GraConst*( Pot_USG[ ID3+dID3[d1] ] - Pot_USG[ ID3                   ] );
               Acc[d2] += (real)0.25*GraConst*( Pot_USG[ ID3+dID3[d2] ] + Pot_USG[ ID3+dID3[d2]+dID3[d1] ]
                                               -Pot_USG[ ID3-dID3[d2] ] - Pot_USG[ ID3-dID3[d2]+dID3[d1] ] );
               Acc[d3] += (real)0.25*GraConst*( Pot_USG[ ID3+dID3[d3] ] + Pot_USG[ ID3+dID3[d3]+dID3[d1] ]
                                               -Pot_USG[ ID3-dID3[d3] ] - Pot_USG[ ID3-dID3[d3]+dID3[d1] ] );
            }

//          store the internal energy density
            eL = ConVar_L[4] - (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            eR = ConVar_R[4] - (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];

//          advance velocity by gravity
            for (int t=0; t<3; t++)
            {
               ConVar_L[t+1] += ConVar_L[0]*Acc[t];
               ConVar_R[t+1] += ConVar_R[0]*Acc[t];
            }

//          update total energy density with the internal energy density fixed
            ConVar_L[4] = eL + (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            ConVar_R[4] = eR + (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY


#        if   ( RSOLVER == EXACT )
         const bool NormPassive_No  = false; // do NOT convert any passive variable to mass fraction for the Riemann solvers
         const bool JeansMinPres_No = false;

         CPU_Con2Pri( ConVar_L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
         CPU_Con2Pri( ConVar_R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

         CPU_RiemannSolver_Exact( d, NULL, NULL, NULL, FC_Flux[ID1][d], PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         CPU_RiemannSolver_Roe ( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         CPU_RiemannSolver_HLLE( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLC )
         CPU_RiemannSolver_HLLC( d, FC_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif
      } // i,j,k
   } // for (int d=0; d<3; d++)

} // FUNCTION : CPU_ComputeFlux



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  (FLU_SCHEME == MHM || MHM_RP || CTU) )
