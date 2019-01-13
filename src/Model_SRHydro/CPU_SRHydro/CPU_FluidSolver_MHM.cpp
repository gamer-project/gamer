#include "GAMER.h"
#include "CUFLU.h"
#include "../../../include/CPU_prototypes.h"

#if (  !defined GPU  &&  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )

#if   ( FLU_SCHEME == MHM_RP )
static void CPU_RiemannPredict( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                const real Half_Flux[][3][NCOMP_TOTAL], real Half_Var[][NCOMP_TOTAL], const real dt,
                                const real dh, const real Gamma, const real MinDens, const real MinPres );
static void CPU_RiemannPredict_Flux( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Half_Flux[][3][NCOMP_TOTAL],
                                     const real Gamma, const real MinPres );
#elif ( FLU_SCHEME == MHM )
static void CPU_HancockPredict( real FC_Var[][6][NCOMP_TOTAL], const real dt, const real dh, const real Gamma,
                                const real C_Var[][ FLU_NXT*FLU_NXT*FLU_NXT ], const real MinDens, const real MinPres );
#endif


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FluidSolver_MHM
// Description :  CPU fluid solver based on the MUSCL-Hancock scheme
//
// Note        :  1. The three-dimensional evolution is achieved by using the unsplit method
//                2. Two half-step prediction schemes are supported, including "MHM" and "MHM_RP"
//                   MHM    : use interpolated face-centered values to calculate the half-step fluxes
//                   MHM_RP : use Riemann solver to calculate the half-step fluxes
//                3. Ref :
//                   MHM    : "Riemann Solvers and Numerical Methods for Fluid Dynamics
//                             - A Practical Introduction ~ by Eleuterio F. Toro"
//                   MHM_RP : Stone & Gardiner, NewA, 14, 139 (2009)
//
// Parameter   : [ 1] Flu_Array_In       : Array storing the input fluid conserved variables
//               [ 2] Flu_Array_Out      : Array to store the output fluid conserved variables
//               [ 3] DE_Array_Out       : Array to store the dual-energy status
//               [ 4] Flux_Array         : Array to store the output fluxes
//               [ 5] Corner_Array       : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//               [ 6] Pot_Array_USG      : Array storing the input potential for UNSPLIT_GRAVITY
//               [ 7] NPatchGroup        : Number of patch groups to be evaluated
//               [ 8] dt                 : Time interval to advance solution
//               [ 9] dh                 : Grid size
//               [10] Gamma              : Ratio of specific heats
//               [11] StoreFlux          : true --> store the coarse-fine fluxes
//               [12] LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                         (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                        vanLeer + generalized MinMod/extrema-preserving) limiter
//               [13] MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//               [14] EP_Coeff           : Coefficient of the extrema-preserving limiter
//               [15] Time               : Current physical time                                     (for UNSPLIT_GRAVITY only)
//               [16] GravityType        : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//               [17] ExtAcc_AuxArray    : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//            [18/19] MinDens/Pres       : Minimum allowed density and pressure
//               [20] DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//               [21] NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                                  is equal to the gas mass density
//               [22] NNorm              : Number of passive scalars to be normalized
//                                         --> Should be set to the global variable "PassiveNorm_NVar"
//               [23] NormIdx            : Target variable indices to be normalized
//                                         --> Should be set to the global variable "PassiveNorm_VarIdx"
//               [24] JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//               [25] JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
void CPU_FluidSolver_MHM( const real Flu_Array_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                          real Flu_Array_Out[][NCOMP_TOTAL][ PS2*PS2*PS2 ],
                          char DE_Array_Out[][ PS2*PS2*PS2 ],
                          real Flux_Array[][9][NCOMP_TOTAL][ PS2*PS2 ],
                          const double Corner_Array[][3],
                          const real Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                          const int NPatchGroup,
                          const real dt,
                          const real dh,
                          const real Gamma,
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
      const bool CorrHalfVel_No  = false;
      const int Max = 2;
      int iteration;
      bool state;
      real MinMod_Coeff_temp;

      real Input[NCOMP_TOTAL];
      int ID1;

//    FC: Face-Centered variables/fluxes
      real (*FC_Var )[6][NCOMP_TOTAL] = new real [ N_FC_VAR*N_FC_VAR*N_FC_VAR    ][6][NCOMP_TOTAL];
      real (*FC_Flux)[3][NCOMP_TOTAL] = new real [ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ][3][NCOMP_TOTAL];   // also used by "Half_Flux"
      real (*PriVar)    [NCOMP_TOTAL] = new real [ FLU_NXT*FLU_NXT*FLU_NXT       ]   [NCOMP_TOTAL];   // also used by "Half_Var"

#     if ( FLU_SCHEME == MHM_RP )
      real (*const Half_Flux)[3][NCOMP_TOTAL] = FC_Flux;
      real (*const Half_Var)    [NCOMP_TOTAL] = PriVar;
#     endif


//    loop over all patch groups
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         iteration = 0;
//       1. half-step prediction
#        if ( FLU_SCHEME == MHM_RP ) // a. use Riemann solver to calculate the half-step fluxes

//       (1.a-1) evaluate the half-step first-order fluxes by Riemann solver
//       check unphysical cell before computing flux
         CPU_RiemannPredict_Flux( Flu_Array_In[P], Half_Flux, Gamma, MinPres );


//       (1.a-2) evaluate the half-step solutions
         CPU_RiemannPredict( Flu_Array_In[P], Half_Flux, Half_Var, dt, dh, Gamma, MinDens, MinPres );
//       check unphysical cell after prediction


//       (1.a-3) conserved variables --> primitive variables
#        if  ( EXTRAPOLATE != CONSERVED_QUANTITIES )
         for (int k=0; k<N_HF_VAR; k++)
         for (int j=0; j<N_HF_VAR; j++)
         for (int i=0; i<N_HF_VAR; i++)
         {
            ID1 = (k*N_HF_VAR + j)*N_HF_VAR + i;

#           if  ( EXTRAPOLATE == FOUR_VELOCITY )
            for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Half_Var[ID1][v];
            CPU_Con2Pri(Input, Half_Var[ID1], Gamma);

#           elif ( EXTRAPOLATE == THREE_VELOCITY )
            CPU_Con2Pri(Half_Var[ID1], Input, Gamma);
            CPU_4Velto3Vel(Input, Half_Var[ID1]);
#           endif
         }
#        endif

         do {
              MinMod_Coeff_temp = ( Max - iteration ) * ( MinMod_Coeff / (real) Max );
//           (1.a-4) evaluate the face-centered values by data reconstruction
             CPU_DataReconstruction( Half_Var, FC_Var, N_HF_VAR, FLU_GHOST_SIZE-2, Gamma, LR_Limiter,
                                     MinMod_Coeff_temp, EP_Coeff, NULL_REAL, NULL_INT, MinDens, MinPres, iteration );
//           check unphysical cell after data reconstruction


//           (1.a-5) primitive face-centered variables --> conserved face-centered variables
#            if ( EXTRAPOLATE != CONSERVED_QUANTITIES )
             for (int k=0; k<N_FC_VAR; k++)
             for (int j=0; j<N_FC_VAR; j++)
             for (int i=0; i<N_FC_VAR; i++)
             {
                ID1 = (k*N_FC_VAR + j)*N_FC_VAR + i;

                for (int f=0; f<6; f++)
                {
#                  if  ( EXTRAPOLATE == FOUR_VELOCITY )
                   for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = FC_Var[ID1][f][v];
                   CPU_Pri2Con( Input, FC_Var[ID1][f], Gamma);

#                  elif ( EXTRAPOLATE == THREE_VELOCITY )
                   CPU_3Velto4Vel(FC_Var[ID1][f], Input);
                   CPU_Pri2Con( Input, FC_Var[ID1][f], Gamma);
#                  endif
                }
             }
#            endif

#            elif ( FLU_SCHEME == MHM ) // b. use interpolated face-centered values to calculate the half-step fluxes

//           (1.b-1) conserved variables --> primitive variables
#            if  ( EXTRAPOLATE != CONSERVED_QUANTITIES ||  defined(CHECK_NEGATIVE_IN_FLUID) )
             for (int k=0; k<FLU_NXT; k++)
             for (int j=0; j<FLU_NXT; j++)
             for (int i=0; i<FLU_NXT; i++)
             {
                ID1 = (k*FLU_NXT + j)*FLU_NXT + i;

#               if  ( EXTRAPOLATE == FOUR_VELOCITY )
                for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Flu_Array_In[P][v][ID1];
#               ifdef CHECK_NEGATIVE_IN_FLUID
                CPU_CheckUnphysical(Input, NULL, __FUNCTION__, __LINE__, true);
#               endif
                CPU_Con2Pri(Input, PriVar[ID1], Gamma);

#               elif ( EXTRAPOLATE == THREE_VELOCITY )
                for (int v=0; v<NCOMP_TOTAL; v++)   PriVar[ID1][v] = Flu_Array_In[P][v][ID1];
#               ifdef CHECK_NEGATIVE_IN_FLUID
                CPU_CheckUnphysical(PriVar[ID1], NULL, __FUNCTION__, __LINE__, true);
#               endif
                CPU_Con2Pri(PriVar[ID1], Input, Gamma);
                CPU_4Velto3Vel(Input, PriVar[ID1]);

#               elif ( EXTRAPOLATE == CONSERVED_QUANTITIES )
                for (int v=0; v<NCOMP_TOTAL; v++)   PriVar[ID1][v] = Flu_Array_In[P][v][ID1];
                CPU_CheckUnphysical(PriVar[ID1], NULL, __FUNCTION__, __LINE__, true);
#               endif
             }
#           endif

         do {
              MinMod_Coeff_temp = ( Max - iteration ) * ( MinMod_Coeff / (real) Max );

//            (1.b-2) evaluate the face-centered values by data reconstruction
              CPU_DataReconstruction( PriVar, FC_Var, FLU_NXT, FLU_GHOST_SIZE-1, Gamma, LR_Limiter,
                                      MinMod_Coeff_temp, EP_Coeff, NULL_REAL, NULL_INT, MinDens, MinPres, iteration );
//            check unphysical cell after data reconstruction

//            (1.b-3) primitive face-centered variables --> conserved face-centered variables
#             if ( EXTRAPOLATE != CONSERVED_QUANTITIES )
              for (int k=0; k<N_FC_VAR; k++)
              for (int j=0; j<N_FC_VAR; j++)
              for (int i=0; i<N_FC_VAR; i++)
              {
                 ID1 = (k*N_FC_VAR + j)*N_FC_VAR + i;

                 for (int f=0; f<6; f++)
                 {
#                   if  ( EXTRAPOLATE == FOUR_VELOCITY )
                    for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = FC_Var[ID1][f][v];
                    CPU_Pri2Con( Input, FC_Var[ID1][f], Gamma);

#                   elif ( EXTRAPOLATE == THREE_VELOCITY )
                    CPU_3Velto4Vel(FC_Var[ID1][f], Input);
                    CPU_Pri2Con( Input, FC_Var[ID1][f], Gamma);
#                   endif
                 }
              }
#             endif


//            (1.b-4) evaluate the half-step solutions
              CPU_HancockPredict( FC_Var, dt, dh, Gamma, Flu_Array_In[P], MinDens, MinPres );

#             endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...

              CPU_ComputeFlux( FC_Var, FC_Flux, N_FL_FLUX, 1, Gamma, CorrHalfVel_No,  NULL, NULL,
                               NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );

//            3. full-step evolution
              state = CPU_FullStepUpdate( Flu_Array_In[P], Flu_Array_Out[P], DE_Array_Out[P],
					  FC_Flux, dt, dh, Gamma, MinDens, MinPres, DualEnergySwitch, 
					  NormPassive, NNorm, NormIdx);

              iteration++;
//              if (iteration >= 2) printf("%f\n",MinMod_Coeff_temp);

//       perform CPU_FullStepUpdate again if cell is unphysical and iteration < MAX
         }while( state && iteration <= Max );

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if(state) printf("Adaptive MinMod_Coeff is fail! %s: %d\n",__FUNCTION__, __LINE__) ;
#        endif


//       4. store the inter-patch fluxes
         if ( StoreFlux )
         CPU_StoreFlux( Flux_Array[P], FC_Flux );

      } // for (int P=0; P<NPatchGroup; P++)

      delete [] FC_Var;
      delete [] FC_Flux;
      delete [] PriVar;

   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//
// Parameter   :  Flu_Array_In : Array storing the input conserved variables
//                Half_Flux    : Array to store the output face-centered fluxes
//                               --> The size is assumed to be N_HF_FLUX^3
//                Gamma        : Ratio of specific heats
//                MinPres      : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannPredict_Flux( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                              real Half_Flux[][3][NCOMP_TOTAL],
                              const real Gamma,
                              const real MinPres )
{
   const int dr[3] = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   int ID1, ID2, dN[3]={ 0 };
   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL];

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      switch ( d )
      {
         case 0 : dN[0] = 0;  dN[1] = 1;  dN[2] = 1;  break;
         case 1 : dN[0] = 1;  dN[1] = 0;  dN[2] = 1;  break;
         case 2 : dN[0] = 1;  dN[1] = 1;  dN[2] = 0;  break;
      }

      for (int k1=0, k2=dN[2];  k1<N_HF_FLUX-dN[2];  k1++, k2++)
      for (int j1=0, j2=dN[1];  j1<N_HF_FLUX-dN[1];  j1++, j2++)
      for (int i1=0, i2=dN[0];  i1<N_HF_FLUX-dN[0];  i1++, i2++)
      {
         ID1 = (k1*N_HF_FLUX + j1)*N_HF_FLUX + i1;
         ID2 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

//       get the left and right states
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = Flu_Array_In[v][ ID2       ];
            ConVar_R[v] = Flu_Array_In[v][ ID2+dr[d] ];
         }
#        ifdef CHECK_MIN_TEMP
         ConVar_L[4] = CPU_CheckMinTempInEngy(ConVar_L);
         ConVar_R[4] = CPU_CheckMinTempInEngy(ConVar_R);
#        endif

//       check unphysical cells
#        ifdef CHECK_NEGATIVE_IN_FLUID
         CPU_CheckUnphysical(ConVar_L, NULL, __FUNCTION__, __LINE__, true);
         CPU_CheckUnphysical(ConVar_R, NULL, __FUNCTION__, __LINE__, true);
#        endif

//       invoke the Riemann solver
#        if ( RSOLVER == HLLC )
         CPU_RiemannSolver_HLLC ( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         CPU_RiemannSolver_HLLE ( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver !!
#        endif
      }
   } // for (int d=0; d<3; d++)

} // FUNCTION : CPU_RiemannPredict_Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_RiemannPredict
// Description :  Evolve the cell-centered variables by half time-step by using the Riemann solvers
//
// Note        :  Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//
// Parameter   :  Flu_Array_In : Array storing the input conserved variables
//                Half_Flux    : Array storing the input face-centered fluxes
//                               --> The size is assumed to be N_HF_FLUX^3
//                Half_Var     : Array to store the output conserved variables
//                               --> The size is assumed to be N_HF_VAR^3
//                dt           : Time interval to advance solution
//                dh           : Grid size
//                Gamma        : Ratio of specific heats
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_RiemannPredict( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                         const real Half_Flux[][3][NCOMP_TOTAL],
                         real Half_Var[][NCOMP_TOTAL], 
                         const real dt, 
                         const real dh, 
                         const real Gamma,
                         const real MinDens, 
                         const real MinPres )
{

   const int  dID3[3]   = { 1, N_HF_FLUX, N_HF_FLUX*N_HF_FLUX };
   const real dt_dh2    = (real)0.5*dt/dh;
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;

   real dF[3][NCOMP_TOTAL];
   int ID1, ID2, ID3;

   for (int k1=0, k2=1;  k1<N_HF_VAR;  k1++, k2++)
   for (int j1=0, j2=1;  j1<N_HF_VAR;  j1++, j2++)
   for (int i1=0, i2=1;  i1<N_HF_VAR;  i1++, i2++)
   {
      ID1 = (k1*N_HF_VAR  + j1)*N_HF_VAR  + i1;
      ID2 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;
      ID3 = (k1*N_HF_FLUX + j1)*N_HF_FLUX + i1;

      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)    dF[d][v] = Half_Flux[ ID3+dID3[d] ][d][v] - Half_Flux[ID3][d][v];

      for (int v=0; v<NCOMP_TOTAL; v++)
         Half_Var[ID1][v] = Flu_Array_In[v][ID2] - dt_dh2*( dF[0][v] + dF[1][v] + dF[2][v] );

//    ensure positive density and pressure
      Half_Var[ID1][0] = FMAX( Half_Var[ID1][0], MinDens );
#     ifdef CHECK_MIN_TEMP
      Half_Var[ID1][4] = CPU_CheckMinTempInEngy( Half_Var[ID1]);
#     endif

#     ifdef CHECK_NEGATIVE_IN_FLUID
      CPU_CheckUnphysical(Half_Var[ID1], NULL, __FUNCTION__, __LINE__, true);
#     endif
   } // i,j,k

} // FUNCTION : CPU_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_HancockPredict
// Description :  Evolve the face-centered variables by half time-step by calculating the face-centered fluxes
//                (no Riemann solver is required)
//
// Note        :  1. Work for the MHM scheme
//                2. Do NOT require data in the neighboring cells
//
// Parameter   :  FC_Var       : Face-centered conserved variables
//                               --> The size is assumed to be N_FC_VAR^3
//                dt           : Time interval to advance solution
//                dh           : Grid size
//                Gamma        : Ratio of specific heats
//                C_Var        : Array storing the conservative variables
//                               --> For checking negative density and pressure
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void CPU_HancockPredict( real FC_Var[][6][NCOMP_TOTAL], const real dt, const real dh, const real Gamma,
                         const real C_Var[][ FLU_NXT*FLU_NXT*FLU_NXT ], const real MinDens, const real MinPres )
{

   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
   const real dt_dh2    = (real)0.5*dt/dh;
   const int  NGhost    = FLU_GHOST_SIZE - 1;

   real Flux[6][NCOMP_TOTAL], dFlux[NCOMP_TOTAL];
   int ID1, ID2;


   for (int k1=0, k2=NGhost;  k1<N_FC_VAR;  k1++, k2++)
   for (int j1=0, j2=NGhost;  j1<N_FC_VAR;  j1++, j2++)
   for (int i1=0, i2=NGhost;  i1<N_FC_VAR;  i1++, i2++)
   {
      ID1 = (k1*N_FC_VAR + j1)*N_FC_VAR + i1;
      ID2 = (k2*FLU_NXT  + j2)*FLU_NXT  + i2;

      for (int f=0; f<6; f++)    CPU_Con2Flux( f/2, Flux[f], FC_Var[ID1][f], Gamma_m1, MinPres );

      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         dFlux[v] = dt_dh2 * ( Flux[1][v] - Flux[0][v] + Flux[3][v] - Flux[2][v] + Flux[5][v] - Flux[4][v] );

         for (int f=0; f<6; f++)    FC_Var[ID1][f][v] -= dFlux[v];
      }

//    check unphyscial results
      for (int f=0; f<6; f++)
      {
         if ( CPU_CheckUnphysical(FC_Var[ID1][f], NULL, __FUNCTION__, __LINE__, false) )
         {
//          set to the values before update
            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               FC_Var[ID1][0][v] = FC_Var[ID1][1][v] = FC_Var[ID1][2][v] = FC_Var[ID1][3][v] =
               FC_Var[ID1][4][v] = FC_Var[ID1][5][v] = C_Var[v][ID2];
            }

            break;
         }
      }

   } // i,j,k

} // FUNCTION : CPU_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )
#endif // #if (  !defined GPU  &&  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
