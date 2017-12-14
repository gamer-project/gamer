#include "GAMER.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



extern void CPU_DataReconstruction( const real PriVar[][NCOMP_TOTAL], real FC_Var[][6][NCOMP_TOTAL], const int NIn, const int NGhost,
                                    const real Gamma, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                    const real EP_Coeff, const real dt, const real dh, const real MinDens, const real MinPres );
extern void CPU_Con2Flux( const int XYZ, real Flux[], const real Input[], const real Gamma_m1, const real MinPres );
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
#if   ( RSOLVER == EXACT )
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
extern real CPU_CheckMinPres( const real InPres, const real MinPres );

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
void CPU_FluidSolver_MHM( const real Flu_Array_In[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
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
#     else
      const bool CorrHalfVel_No  = false;
#     endif

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

//       1. half-step prediction
#        if ( FLU_SCHEME == MHM_RP ) // a. use Riemann solver to calculate the half-step fluxes

//       (1.a-1) evaluate the half-step first-order fluxes by Riemann solver
         CPU_RiemannPredict_Flux( Flu_Array_In[P], Half_Flux, Gamma, MinPres );


//       (1.a-2) evaluate the half-step solutions
         CPU_RiemannPredict( Flu_Array_In[P], Half_Flux, Half_Var, dt, dh, Gamma, MinDens, MinPres );


//       (1.a-3) conserved variables --> primitive variables
         for (int k=0; k<N_HF_VAR; k++)
         for (int j=0; j<N_HF_VAR; j++)
         for (int i=0; i<N_HF_VAR; i++)
         {
            ID1 = (k*N_HF_VAR + j)*N_HF_VAR + i;

            for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Half_Var[ID1][v];

            CPU_Con2Pri( Input, Half_Var[ID1], Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                         JeansMinPres, JeansMinPres_Coeff );
         }


//       (1.a-4) evaluate the face-centered values by data reconstruction
         CPU_DataReconstruction( Half_Var, FC_Var, N_HF_VAR, FLU_GHOST_SIZE-2, Gamma, LR_Limiter,
                                 MinMod_Coeff, EP_Coeff, NULL_REAL, NULL_INT, MinDens, MinPres );


//       (1.a-5) primitive face-centered variables --> conserved face-centered variables
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

#        elif ( FLU_SCHEME == MHM ) // b. use interpolated face-centered values to calculate the half-step fluxes

//       (1.b-1) conserved variables --> primitive variables
         for (int k=0; k<FLU_NXT; k++)
         for (int j=0; j<FLU_NXT; j++)
         for (int i=0; i<FLU_NXT; i++)
         {
            ID1 = (k*FLU_NXT + j)*FLU_NXT + i;

            for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Flu_Array_In[P][v][ID1];

            CPU_Con2Pri( Input, PriVar[ID1], Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                         JeansMinPres, JeansMinPres_Coeff );
         }


//       (1.b-2) evaluate the face-centered values by data reconstruction
         CPU_DataReconstruction( PriVar, FC_Var, FLU_NXT, FLU_GHOST_SIZE-1, Gamma, LR_Limiter,
                                 MinMod_Coeff, EP_Coeff, NULL_REAL, NULL_INT, MinDens, MinPres );


//       (1.b-3) primitive face-centered variables --> conserved face-centered variables
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


//       (1.b-4) evaluate the half-step solutions
         CPU_HancockPredict( FC_Var, dt, dh, Gamma, Flu_Array_In[P], MinDens, MinPres );

#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


//       2. evaluate the full-step fluxes
#        ifdef UNSPLIT_GRAVITY
         CPU_ComputeFlux( FC_Var, FC_Flux, N_FL_FLUX, 1, Gamma, CorrHalfVel_Yes, Pot_Array_USG[P][0][0], Corner_Array[P],
                          dt, dh, Time, GravityType, ExtAcc_AuxArray, MinPres );
#        else
         CPU_ComputeFlux( FC_Var, FC_Flux, N_FL_FLUX, 1, Gamma, CorrHalfVel_No,  NULL, NULL,
                          NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );
#        endif


//       3. full-step evolution
         CPU_FullStepUpdate( Flu_Array_In[P], Flu_Array_Out[P], DE_Array_Out[P],
                             FC_Flux, dt, dh, Gamma, MinDens, MinPres, DualEnergySwitch,
                             NormPassive, NNorm, NormIdx );


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
void CPU_RiemannPredict_Flux( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Half_Flux[][3][NCOMP_TOTAL],
                              const real Gamma, const real MinPres )
{

   const int dr[3] = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   int ID1, ID2, dN[3]={ 0 };
   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL];

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif


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

//       invoke the Riemann solver
#        if   ( RSOLVER == EXACT )
         const bool NormPassive_No  = false;  // do NOT convert any passive variable to mass fraction for the Riemann solvers
         const bool JeansMinPres_No = false;

         CPU_Con2Pri( ConVar_L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
         CPU_Con2Pri( ConVar_R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

         CPU_RiemannSolver_Exact( d, NULL, NULL, NULL, Half_Flux[ID1][d], PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         CPU_RiemannSolver_Roe ( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         CPU_RiemannSolver_HLLE( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLC )
         CPU_RiemannSolver_HLLC( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
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
void CPU_RiemannPredict( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], const real Half_Flux[][3][NCOMP_TOTAL],
                         real Half_Var[][NCOMP_TOTAL], const real dt, const real dh, const real Gamma,
                         const real MinDens, const real MinPres )
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
      Half_Var[ID1][4] = CPU_CheckMinPresInEngy( Half_Var[ID1][0], Half_Var[ID1][1], Half_Var[ID1][2],
                                                 Half_Var[ID1][3], Half_Var[ID1][4], Gamma_m1, _Gamma_m1, MinPres );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      Half_Var[ID1][v] = FMAX( Half_Var[ID1][v], TINY_NUMBER );
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

//    check the negative density and energy
      for (int f=0; f<6; f++)
      {
         if ( FC_Var[ID1][f][0] <= (real)0.0  ||  FC_Var[ID1][f][4] <= (real)0.0 )
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

//    ensure positive density and pressure
      for (int f=0; f<6; f++)
      {
         FC_Var[ID1][f][0] = FMAX( FC_Var[ID1][f][0], MinDens );
         FC_Var[ID1][f][4] = CPU_CheckMinPresInEngy( FC_Var[ID1][f][0], FC_Var[ID1][f][1], FC_Var[ID1][f][2],
                                                     FC_Var[ID1][f][3], FC_Var[ID1][f][4], Gamma_m1, _Gamma_m1, MinPres );
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         FC_Var[ID1][f][v] = FMAX( FC_Var[ID1][f][v], TINY_NUMBER );
#        endif
      }
   } // i,j,k

} // FUNCTION : CPU_HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



#endif // #if (  !defined GPU  &&  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
