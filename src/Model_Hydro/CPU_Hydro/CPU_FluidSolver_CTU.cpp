#include "GAMER.h"
#include "CUFLU.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )



extern void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                                    const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                    const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres );
extern void CPU_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                         const bool NormPassive, const int NNorm, const int NormIdx[],
                         const bool JeansMinPres, const real JeansMinPres_Coeff );
extern void CPU_Pri2Con( const real In[], real Out[], const real _Gamma_m1,
                         const bool NormPassive, const int NNorm, const int NormIdx[] );
extern void CPU_ComputeFlux( const real FC_Var[][6][NCOMP_TOTAL], real FC_Flux[][3][NCOMP_TOTAL], const int NFlux, const int Gap,
                             const real Gamma, const bool CorrHalfVel, const real Pot_USG[], const double Corner[],
                             const real dt, const real dh, const double Time, const OptGravityType_t GravityType,
                             const double ExtAcc_AuxArray[], const real MinPres );
extern void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                                const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                                const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                                const bool NormPassive, const int NNorm, const int NormIdx[] );
extern void CPU_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL] );
extern real CPU_CheckMinPres( const real InPres, const real MinPres );

static void TGradient_Correction( real FC_Var[][6][NCOMP_TOTAL], const real FC_Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                                  const real Gamma_m1, const real _Gamma_m1, const real MinDens, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver_CTU
// Description :  CPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  Ref : Stone et al., ApJS, 178, 137 (2008)
//
// Parameter   :  Flu_Array_In       : Array storing the input fluid variables
//                Flu_Array_Out      : Array to store the output fluid variables
//                DE_Array_Out       : Array to store the dual-energy status
//                Flux_Array         : Array to store the output fluxes
//                Corner_Array       : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                Pot_Array_USG      : Array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup        : Number of patch groups to be evaluated
//                dt                 : Time interval to advance solution
//                dh                 : Grid size
//                Gamma              : Ratio of specific heats
//                StoreFlux          : true --> store the coarse-fine fluxes
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                EP_Coeff           : Coefficient of the extrema-preserving limiter
//                Time               : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType        : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray    : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//                MinDens/Pres       : Minimum allowed density and pressure
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                              is equal to the gas mass density
//                NNorm              : Number of passive scalars to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver_CTU( const real Flu_Array_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out[][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          char DE_Array_Out[][ PS2*PS2*PS2 ],
                          real Flux_Array[][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup, const real dt, const real dh, const real Gamma,
                          const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                          const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
                          const double ExtAcc_AuxArray[], const real MinDens, const real MinPres,
                          const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                          const bool JeansMinPres, const real JeansMinPres_Coeff )
{

// check
#  ifdef GAMER_DEBUG
   if ( LR_Limiter != VANLEER  &&  LR_Limiter != GMINMOD  &&  LR_Limiter != ALBADA  &&  LR_Limiter != EXTPRE  &&
        LR_Limiter != VL_GMINMOD )
      Aux_Error( ERROR_INFO, "unsupported reconstruction limiter (%d) !!\n", LR_Limiter );
#  endif


#  pragma omp parallel
   {
      const real  Gamma_m1       = Gamma - (real)1.0;
      const real _Gamma_m1       = (real)1.0 / Gamma_m1;
#     ifdef UNSPLIT_GRAVITY
      const bool CorrHalfVel_Yes = true;
#     endif
      const bool CorrHalfVel_No  = false;

      real Input[NCOMP_TOTAL];
      int ID1;

//    FC: Face-Centered variables/fluxes
      real (*FC_Var )[6][NCOMP_TOTAL] = new real [ N_FC_VAR*N_FC_VAR*N_FC_VAR    ][6][NCOMP_TOTAL];
      real (*FC_Flux)[3][NCOMP_TOTAL] = new real [ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ][3][NCOMP_TOTAL];
      real (*PriVar)    [NCOMP_TOTAL] = new real [ FLU_NXT*FLU_NXT*FLU_NXT       ]   [NCOMP_TOTAL];


//    loop over all patch groups
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {

//       1. conserved variables --> primitive variables
         for (int k=0; k<FLU_NXT; k++)
         for (int j=0; j<FLU_NXT; j++)
         for (int i=0; i<FLU_NXT; i++)
         {
            ID1 = (k*FLU_NXT + j)*FLU_NXT + i;

            for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Flu_Array_In[P][v][ID1];

            CPU_Con2Pri( Input, PriVar[ID1], Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                         JeansMinPres, JeansMinPres_Coeff );
         }


//       2. evaluate the face-centered values at the half time-step
         CPU_DataReconstruction( PriVar, FC_Var, FLU_NXT, FLU_GHOST_SIZE-1, Gamma, LR_Limiter,
                                 MinMod_Coeff, EP_Coeff, dt, dh, MinDens, MinPres );


//       3. primitive face-centered variables --> conserved face-centered variables
         for (int k=0; k<N_FC_VAR; k++)
         for (int j=0; j<N_FC_VAR; j++)
         for (int i=0; i<N_FC_VAR; i++)
         {
            ID1 = (k*N_FC_VAR + j)*N_FC_VAR + i;

            for (int f=0; f<6; f++)
            {
               for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = FC_Var[ID1][f][v];

               CPU_Pri2Con( Input, FC_Var[ID1][f], _Gamma_m1, NormPassive, NNorm, NormIdx );
            }
         }


//       4. evaluate the face-centered half-step fluxes by solving the Riemann problem
         CPU_ComputeFlux( FC_Var, FC_Flux, N_HF_FLUX, 0, Gamma, CorrHalfVel_No, NULL, NULL,
                          NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );


//       5. correct the face-centered variables by the transverse flux gradients
         TGradient_Correction( FC_Var, FC_Flux, dt, dh, Gamma_m1, _Gamma_m1, MinDens, MinPres );


//       6. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
#        ifdef UNSPLIT_GRAVITY
         CPU_ComputeFlux( FC_Var, FC_Flux, N_FL_FLUX, 1, Gamma, CorrHalfVel_Yes, Pot_Array_USG[P][0][0], Corner_Array[P],
                          dt, dh, Time, GravityType, ExtAcc_AuxArray, MinPres );
#        else
         CPU_ComputeFlux( FC_Var, FC_Flux, N_FL_FLUX, 1, Gamma, CorrHalfVel_No,  NULL, NULL,
                          NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );
#        endif


//       7. full-step evolution
         CPU_FullStepUpdate( Flu_Array_In[P], Flu_Array_Out[P], DE_Array_Out[P],
                             FC_Flux, dt, dh, Gamma, MinDens, MinPres, DualEnergySwitch,
                             NormPassive, NNorm, NormIdx );


//       8. store the inter-patch fluxes
         if ( StoreFlux )
         CPU_StoreFlux( Flux_Array[P], FC_Flux );

      } // for (int P=0; P<NPatchGroup; P++)

      delete [] FC_Var;
      delete [] FC_Flux;
      delete [] PriVar;

   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  TGradient_Correction
// Description :  1. Correct the face-centered variables by the transverse flux gradients
//                2. This function assumes that "N_FC_VAR == N_FC_FLUX == NGrid"
//
// Parameter   :  FC_Var       : Array to store the input and output face-centered conserved variables
//                               --> Size is assumed to be N_FC_VAR
//                FC_Flux      : Array storing the input face-centered fluxes
//                               --> Size is assumed to be N_FC_FLUX
//                dt           : Time interval to advance solution
//                dh           : Grid size
//                Gamma_m1     : Gamma - 1
//                _Gamma_m1    : 1/(Gamma - 1)
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void TGradient_Correction( real FC_Var[][6][NCOMP_TOTAL], const real FC_Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                           const real Gamma_m1, const real _Gamma_m1, const real MinDens, const real MinPres )
{

   const int  NGrid  = N_FC_VAR;    // size of the arrays FC_Var and FC_Flux in each direction
   const int  dID[3] = { 1, NGrid, NGrid*NGrid };
   const real dt_dh2 = (real)0.5*dt/dh;

   real Correct, TGrad1, TGrad2;
   int dL, dR, ID, ID_L1, ID_L2, ID_R, TDir1, TDir2, Gap[3]={0};


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      dL    = 2*d;
      dR    = dL+1;
      TDir1 = (d+1)%3;  // transverse direction ONE
      TDir2 = (d+2)%3;  // transverse direction TWO

      switch ( d )
      {
         case 0 : Gap[0] = 0;   Gap[1] = 1;   Gap[2] = 1;   break;
         case 1 : Gap[0] = 1;   Gap[1] = 0;   Gap[2] = 1;   break;
         case 2 : Gap[0] = 1;   Gap[1] = 1;   Gap[2] = 0;   break;
      }

      for (int k=Gap[2]; k<NGrid-Gap[2]; k++)
      for (int j=Gap[1]; j<NGrid-Gap[1]; j++)
      for (int i=Gap[0]; i<NGrid-Gap[0]; i++)
      {
         ID    = (k*NGrid + j)*NGrid + i;
         ID_R  = ID;
         ID_L1 = ID_R - dID[TDir1];
         ID_L2 = ID_R - dID[TDir2];

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            TGrad1  = FC_Flux[ID_R][TDir1][v] - FC_Flux[ID_L1][TDir1][v];
            TGrad2  = FC_Flux[ID_R][TDir2][v] - FC_Flux[ID_L2][TDir2][v];
            Correct = -dt_dh2*( TGrad1 + TGrad2 );

            FC_Var[ID][dL][v] += Correct;
            FC_Var[ID][dR][v] += Correct;
         }

//       ensure positive density and pressure
         FC_Var[ID][dL][0] = FMAX( FC_Var[ID][dL][0], MinDens );
         FC_Var[ID][dR][0] = FMAX( FC_Var[ID][dR][0], MinDens );

         FC_Var[ID][dL][4] = CPU_CheckMinPresInEngy( FC_Var[ID][dL][0], FC_Var[ID][dL][1], FC_Var[ID][dL][2],
                                                     FC_Var[ID][dL][3], FC_Var[ID][dL][4], Gamma_m1, _Gamma_m1, MinPres );
         FC_Var[ID][dR][4] = CPU_CheckMinPresInEngy( FC_Var[ID][dR][0], FC_Var[ID][dR][1], FC_Var[ID][dR][2],
                                                     FC_Var[ID][dR][3], FC_Var[ID][dR][4], Gamma_m1, _Gamma_m1, MinPres );

#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++) {
         FC_Var[ID][dL][v] = FMAX( FC_Var[ID][dL][v], TINY_NUMBER );
         FC_Var[ID][dR][v] = FMAX( FC_Var[ID][dR][v], TINY_NUMBER ); }
#        endif
      } // i,j,k
   } // for (int d=0; d<3; d++)

} // FUNCTION : TGradient_Correction



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )
