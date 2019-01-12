#include "GAMER.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == SR_HYDRO  && ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )

bool CPU_CheckUnphysical( const real Con[], const real Pri[], const char s[], const int line, bool show);
real CPU_CheckMinTempInEngy (const real Con[]);
static bool boolean;
//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Parameter   :  [ 1] Input            : Array storing the input initial data
//                [ 2] Output           : Array to store the ouptut updated data
//                [ 3] DE_Status        : Array to store the dual-energy status
//                [ 4] Flux             : Array storing the input face-centered flux
//                                        --> Size is assumed to be N_FL_FLUX^3
//                [ 5] dt               : Time interval to advance solution
//                [ 6] dh               : Grid size
//                [ 7] Gamma            : Ratio of specific heats
//                [ 8] MinDens          : Minimum allowed density
//                [ 9] MinPres          : Minimum allowed pressure
//                [10] DualEnergySwitch : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                [11] NormPassive      : true --> normalize passive scalars so that the sum of their mass density
//                                                 is equal to the gas mass density
//                [12] NNorm            : Number of passive scalars to be normalized
//                                        --> Should be set to the global variable "PassiveNorm_NVar"
//                [13] NormIdx          : Target variable indices to be normalized
//                                        --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                         const real Flux[][3][NCOMP_TOTAL], const real dt, const real dh,
                         const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                         const bool NormPassive, const int NNorm, const int NormIdx[], bool state )
{
   const int  dID1[3]   = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh     = dt/dh;

   int  ID1, ID2, ID3;
   real dF[3][NCOMP_TOTAL];

   for (int k1=0, k2=FLU_GHOST_SIZE;  k1<PS2;  k1++, k2++)
   for (int j1=0, j2=FLU_GHOST_SIZE;  j1<PS2;  j1++, j2++)
   for (int i1=0, i2=FLU_GHOST_SIZE;  i1<PS2;  i1++, i2++)
   {

      ID1 = (k1*N_FL_FLUX + j1)*N_FL_FLUX + i1;
      ID2 = (k1*PS2       + j1)*PS2       + i1;
      ID3 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

#     ifdef CHECK_NEGATIVE_IN_FLUID
      real Con[NCOMP_FLUID];
      for (int v = 0;v<NCOMP_FLUID;v++) Con[v] = Input[v][ID3];
      boolean = CPU_CheckUnphysical(Con, NULL, __FUNCTION__, __LINE__, true);
#     endif

      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)   
         dF[d][v] = Flux[ ID1+dID1[d] ][d][v] - Flux[ID1][d][v];

      for (int v=0; v<NCOMP_TOTAL; v++)
         Output[v][ID2] = Input[v][ID3] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );

#     if ( defined CHECK_MIN_TEMP ) || (defined CHECK_NEGATIVE_IN_FLUID )
      real Cons[NCOMP_FLUID];
      for (int v = 0;v<NCOMP_FLUID;v++) Cons[v] = Output[v][ID2];
#     ifdef CHECK_MIN_TEMP
      Output[ENGY][ID2] = CPU_CheckMinTempInEngy( Cons );
      Cons[ENGY] = Output[ENGY][ID2];
#     endif
#     ifdef CHECK_NEGATIVE_IN_FLUID
      state = CPU_CheckUnphysical(Cons, NULL, __FUNCTION__, __LINE__, true);
#     endif
#     else
      state = CPU_CheckUnphysical(Cons, NULL, __FUNCTION__, __LINE__, false);
#     endif
   } // i,j,k

} // FUNCTION : CPU_FullStepUpdate



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_StoreFlux
// Description :  Store the inter-patch fluxes for the AMR fix-up operation
//
// Parameter   :  Output   : Array to store the inter-patch fluxes
//                FC_Flux  : Array storing the face-centered fluxes
//                           --> Size is assumed to be N_FL_FLUX^3
//-------------------------------------------------------------------------------------------------------
void CPU_StoreFlux( real Flux_Array[][NCOMP_TOTAL][ PS2*PS2 ], const real FC_Flux[][3][NCOMP_TOTAL]  )
{

   int Face, ID1, ID2[9];

   for (int m=0; m<PS2; m++)
   for (int n=0; n<PS2; n++)
   {
      ID1    = m*PS2 + n;

      ID2[0] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + 0;
      ID2[1] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS1;
      ID2[2] = (  m*N_FL_FLUX +   n)*N_FL_FLUX + PS2;
      ID2[3] = (  m*N_FL_FLUX +   0)*N_FL_FLUX + n;
      ID2[4] = (  m*N_FL_FLUX + PS1)*N_FL_FLUX + n;
      ID2[5] = (  m*N_FL_FLUX + PS2)*N_FL_FLUX + n;
      ID2[6] = (  0*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[7] = (PS1*N_FL_FLUX +   m)*N_FL_FLUX + n;
      ID2[8] = (PS2*N_FL_FLUX +   m)*N_FL_FLUX + n;

      for (int t=0; t<9; t++)
      {
         Face = t/3;

         for (int v=0; v<NCOMP_TOTAL; v++)   Flux_Array[t][v][ID1] = FC_Flux[ ID2[t] ][Face][v];
      }
   }

} // FUNCTION : CPU_StoreFlux



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

