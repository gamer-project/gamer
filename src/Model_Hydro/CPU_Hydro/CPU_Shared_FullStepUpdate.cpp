#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

#if (  !defined GPU  &&  MODEL == HYDRO  &&  \
       ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU )  )




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Parameter   :  Input    : Array storing the input initial data
//                Output   : Array to store the ouptut updated data
//                Flux     : Array storing the input face-centered flux
//                           --> Size is assumed to be N_FL_FLUX^3
//                dt       : Time interval to advance solution
//                dh       : Grid size
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
void CPU_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ],
                         const real Flux[][3][5], const real dt, const real dh,
                         const real Gamma )
{

   const int  dID1[3] = { 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX };
   const real dt_dh   = dt/dh;

   int ID1, ID2, ID3;
   real dF[3][5];

   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;


   for (int k1=0, k2=FLU_GHOST_SIZE;  k1<PS2;  k1++, k2++)
   for (int j1=0, j2=FLU_GHOST_SIZE;  j1<PS2;  j1++, j2++)
   for (int i1=0, i2=FLU_GHOST_SIZE;  i1<PS2;  i1++, i2++)
   {

      ID1 = (k1*N_FL_FLUX + j1)*N_FL_FLUX + i1;
      ID2 = (k1*PS2       + j1)*PS2       + i1;
      ID3 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

      for (int d=0; d<3; d++)
      for (int v=0; v<5; v++)    dF[d][v] = Flux[ ID1+dID1[d] ][d][v] - Flux[ID1][d][v];

      for (int v=0; v<5; v++)
         Output[v][ID2] = Input[v][ID3] - dt_dh*( dF[0][v] + dF[1][v] + dF[2][v] );

//    we no longer check negative density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//    --> because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
      /*
//    ensure positive density and pressure
      Output[0][ID2] = FMAX( Output[0][ID2], MinDens );
      Output[4][ID2] = CPU_CheckMinPresInEngy( Output[0][ID2], Output[1][ID2], Output[2][ID2], Output[3][ID2], Output[4][ID2],
                                               Gamma_m1, _Gamma_m1, MinPres );
      */

//    check the negative density and energy
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CPU_CheckNegative(Output[0][ID2]) )
         Aux_Message( stderr, "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      Output[0][ID2], __FILE__, __LINE__, __FUNCTION__ );

      if ( CPU_CheckNegative(Output[4][ID2]) )
         Aux_Message( stderr, "ERROR : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                      Output[4][ID2], __FILE__, __LINE__, __FUNCTION__ );
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
void CPU_StoreFlux( real Flux_Array[][5][ PS2*PS2 ], const real FC_Flux[][3][5]  )
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

         for (int v=0; v<5; v++)    Flux_Array[t][v][ID1] = FC_Flux[ ID2[t] ][Face][v];
      }
   }

} // FUNCTION : CPU_StoreFlux



#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

