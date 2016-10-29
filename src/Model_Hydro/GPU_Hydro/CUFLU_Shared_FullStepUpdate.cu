#include "Copyright.h"
#ifndef __CUFLU_FULLSTEPUPDATE_CU__
#define __CUFLU_FULLSTEPUPDATE_CU__



#include "Macro.h"
#include "CUFLU.h"

static __device__ void CUFLU_FullStepUpdate( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                             real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                             const real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                             const real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                             const real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                             const real dt, const real _dh, const real Gamma );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FullStepUpdate
// Description :  Evaluate the full-step solution
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  g_Fluid_In  : Global memory array storing the input fluid variables
//                g_Fluid_Out : Global memory array to store the output fluid variables
//                g_FC_Flux_x : Global memory array storing the input face-centered fluxes in the x direction
//                g_FC_Flux_y : Global memory array storing the input face-centered fluxes in the y direction
//                g_FC_Flux_z : Global memory array storing the input face-centered fluxes in the z direction
//                dt          : Time interval to advance solution
//                _dh         : 1 / grid size
//                Gamma       : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_FullStepUpdate( const real g_Fluid_In [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                            real g_Fluid_Out[][5][ PS2*PS2*PS2 ],
                                      const real g_FC_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                      const real g_FC_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                      const real g_FC_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                      const real dt, const real _dh, const real Gamma )
{

   /*
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
   */
   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID_Out  = blockDim.x;
   const uint3 dID_F    = make_uint3( 1, N_FL_FLUX, N_FL_FLUX*N_FL_FLUX );
   const real  dt_dh    = dt*_dh;

   uint   ID_In, ID_F, ID_Out;
   uint3  ID3d;
   FluVar ConVar;
   real   FluxDiff;


   ID_Out = tx;

// loop over all cells
   while ( ID_Out < PS2*PS2*PS2 )
   {
      ID3d.x = ID_Out%PS2;
      ID3d.y = ID_Out%(PS2*PS2)/PS2;
      ID3d.z = ID_Out/(PS2*PS2);
      ID_In  = __umul24( __umul24( ID3d.z+FLU_GHOST_SIZE, FLU_NXT  ) + ID3d.y+FLU_GHOST_SIZE, FLU_NXT  )
               + ID3d.x+FLU_GHOST_SIZE;
      ID_F   = __umul24( __umul24( ID3d.z, N_FL_FLUX ) + ID3d.y, N_FL_FLUX ) + ID3d.x;


//    get the full-step solution
      ConVar.Rho = g_Fluid_In[bx][0][ID_In];
      ConVar.Px  = g_Fluid_In[bx][1][ID_In];
      ConVar.Py  = g_Fluid_In[bx][2][ID_In];
      ConVar.Pz  = g_Fluid_In[bx][3][ID_In];
      ConVar.Egy = g_Fluid_In[bx][4][ID_In];

#     define Update( comp, v )                                                                        \
      {                                                                                               \
         FluxDiff = dt_dh * (  g_FC_Flux_x[bx][v][ID_F+dID_F.x] - g_FC_Flux_x[bx][v][ID_F] +          \
                               g_FC_Flux_y[bx][v][ID_F+dID_F.y] - g_FC_Flux_y[bx][v][ID_F] +          \
                               g_FC_Flux_z[bx][v][ID_F+dID_F.z] - g_FC_Flux_z[bx][v][ID_F]  );        \
         ConVar.comp -= FluxDiff;                                                                     \
      } // Update

      Update( Rho, 0 );
      Update( Px,  1 );
      Update( Py,  2 );
      Update( Pz,  3 );
      Update( Egy, 4 );

#     undef Update


//    we no longer check negative density and pressure here
//    --> these checks have been moved to Flu_Close()->CorrectUnphysical()
//    --> because we want to apply 1st-order-flux correction BEFORE setting a minimum density and pressure
      /*
//    enforce positive density and pressure
      ConVar.Rho = FMAX( ConVar.Rho, MinDens );
      ConVar.Egy = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );
      */


//    check negative density and energy
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( CUFLU_CheckNegative(ConVar.Rho) )
      printf( "WARNING : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              ConVar.Rho, __FILE__, __LINE__, __FUNCTION__ );

      if ( CUFLU_CheckNegative(ConVar.Egy) )
      printf( "WARNING : negative energy (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              ConVar.Egy, __FILE__, __LINE__, __FUNCTION__ );
#     endif


//    save the updated data back to the output global array
      g_Fluid_Out[bx][0][ID_Out] = ConVar.Rho;
      g_Fluid_Out[bx][1][ID_Out] = ConVar.Px;
      g_Fluid_Out[bx][2][ID_Out] = ConVar.Py;
      g_Fluid_Out[bx][3][ID_Out] = ConVar.Pz;
      g_Fluid_Out[bx][4][ID_Out] = ConVar.Egy;


      ID_Out += dID_Out;

   } // while ( ID_Out < PS2*PS2*PS2 )

} // FUNCTION : CUFLU_FullStepUpdate



#endif // #ifndef __CUFLU_FULLSTEPUPDATE_CU__
