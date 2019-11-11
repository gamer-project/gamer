#include "CUPOT.h"
#include "CUFLU.h"

#if ( MODEL == SR_HYDRO  &&  defined GRAVITY )

// external functions
#ifdef __CUDACC__

#include "../GPU_SRHydro/CUFLU_Shared_FluUtility.cu"

#else
# include "../../../include/SRHydroPrototypes.h"
#endif // #ifdef __CUDACC__ ... else ...


// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"


// variables reside in constant memory
__constant__ double c_ExtAcc_AuxArray[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_SetConstMem_HydroGravitySolver
// Description :  Set the constant memory used by CUPOT_HydroGravitySolver()
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUPOT_SetConstMem_HydroGravitySolver( double h_ExtAcc_AuxArray[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( c_ExtAcc_AuxArray, h_ExtAcc_AuxArray, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_SetConstMem_HydroGravitySolver

#endif // ifdef __CUDACC__




//-----------------------------------------------------------------------------------------
// Function    :  CPU/CUPOT_HydroGravitySolver
// Description :  Advances the momentum and energy density of a group of patches by gravitational acceleration
//                (including external gravity)
//
// Note        :  1. Currently this function does NOT ensure the consistency between Etot-Ekin and
//                   the dual-energy variable (either internal energy of entropy)
//                   --> This consistency breaks only for cells with the dual-energy status labelled
//                       as DE_UPDATED_BY_ETOT_GRA
//                   --> We restore this consistency in Gra_Close()
//                2. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_New   : Array to store the input and output fluid variables
//                g_Pot_Array_New   : Array storing the input potential (at the current step)
//                                    --> _New: to be distinguishable from g_Pot_Array_USG[], which is defined at the previous step
//                g_Corner_Array    : Array storing the physical corner coordinates of each patch
//                g_Pot_Array_USG   : Array storing the input potential          for UNSPLIT_GRAVITY (at the previous step)
//                g_Flu_Array_USG   : Array storing the input density + momentum for UNSPLIT_GRAVITY (at the previous step)
//                g_DE_Array        : Array storing the dual-energy status (for both input and output)
//                NPatchGroup       : Number of input patch groups (for CPU only)
//                dt                : Time interval to advance solution
//                dh                : Cell size
//                P5_Gradient       : Use 5-points stencil to evaluate the potential gradient
//                GravityType       : Types of gravity --> self-gravity, external gravity, both
//                c_ExtAcc_AuxArray : Auxiliary array for adding external acceleration (for CPU only)
//                                    --> When using GPU, this array is stored in the constant memory and does
//                                        not need to be passed as a function argument
//                                        --> Declared on top of this file with the prefix "c_" to
//                                            highlight that this is a constant variable on GPU
//                TimeNew           : Physical time at the current  step (for the external gravity solver)
//                TimeOld           : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//                MinEint           : Minimum allowed internal energy (== MIN_PRES / (GAMMA-1))
//
// Return      :  g_Flu_Array_New, g_DE_Array
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUPOT_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN_USG][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const real dt, const real dh, const bool P5_Gradient,
   const OptGravityType_t GravityType,
   const double TimeNew, const double TimeOld, const real MinEint )
#else
void CPU_HydroGravitySolver(
         real   g_Flu_Array_New[][GRA_NIN][ CUBE(PS1) ],
   const real   g_Pot_Array_New[][ CUBE(GRA_NXT) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_G) ],
   const real   g_Flu_Array_USG[][GRA_NIN_USG][ CUBE(PS1) ],
         char   g_DE_Array     [][ CUBE(PS1) ],
   const int NPatchGroup,
   const real dt, const real dh, const bool P5_Gradient,
   const OptGravityType_t GravityType, const double c_ExtAcc_AuxArray[],
   const double TimeNew, const double TimeOld, const real MinEint )
#endif
{

// check
#  ifdef GAMER_DEBUG
   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeNew < 0.0 )
      printf( "ERROR : incorrect TimeNew (%14.7e) !!\n", TimeNew );

#  ifdef UNSPLIT_GRAVITY
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   if ( g_Pot_Array_USG == NULL  ||  g_Flu_Array_USG == NULL )
      printf( "ERROR : g_Pot_Array_USG == NULL  ||  g_Flu_Array_USG == NULL !!\n" );

   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeOld >= TimeNew  ||  TimeOld < 0.0 )
      printf( "ERROR : incorrect time (TimeOld %14.7e, TimeNew = %14.7e) !!\n", TimeOld, TimeNew );
#  endif

#  endif // #ifdef GAMER_DEBUG


   const real Gra_Const   = ( P5_Gradient ) ? -dt/(12.0*dh) : -dt/(2.0*dh);
   const int  PS1_sqr     = SQR(PS1);
   const int  didx_new[3] = { 1, GRA_NXT,   SQR(GRA_NXT)   };
#  ifdef UNSPLIT_GRAVITY
   const int  didx_old[3] = { 1, USG_NXT_G, SQR(USG_NXT_G) };
#  endif


// load potential from global to shared memory to improve the GPU performance
#  ifdef __CUDACC__
   __shared__ real s_pot_new[ CUBE(GRA_NXT) ];
#  ifdef UNSPLIT_GRAVITY
   __shared__ real s_pot_old[ CUBE(USG_NXT_G) ];
#  endif

   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      for (int t=threadIdx.x; t<CUBE(GRA_NXT); t+=GRA_BLOCK_SIZE)
         s_pot_new[t] = g_Pot_Array_New[blockIdx.x][t];

#     ifdef UNSPLIT_GRAVITY
      for (int t=threadIdx.x; t<CUBE(USG_NXT_G); t+=GRA_BLOCK_SIZE)
         s_pot_old[t] = g_Pot_Array_USG[blockIdx.x][t];
#     endif
   }

   __syncthreads();
#  endif // #ifdef __CUDACC__


// loop over all patches
// --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//     to work on different patches
#  ifdef __CUDACC__
   const int P = blockIdx.x;
#  else
#  pragma omp parallel for schedule( runtime )
   for (int P=0; P<NPatchGroup*8; P++)
#  endif
   {
//    point to the potential array of the target patch
#     ifdef __CUDACC__
      const real *const pot_new = s_pot_new;
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = s_pot_old;
#     endif
#     else // #ifdef __CUDACC__
      const real *const pot_new = g_Pot_Array_New[P];
#     ifdef UNSPLIT_GRAVITY
      const real *const pot_old = g_Pot_Array_USG[P];
#     endif
#     endif // #ifdef __CUDACC__ ... else ...


//    loop over all cells of the target patch
//    _g0: indices for the arrays without any ghost zone
      CGPU_LOOP( idx_g0, CUBE(PS1) )
      {
         real acc_new[3]={0.0, 0.0, 0.0};
         real Con_new[NCOMP_FLUID], Pri_new[NCOMP_FLUID];

		 real LorentzFactor_new, n_new , Ux_new, Uy_new, Uz_new, P_new;
		 real Uxx_new, Uyy_new, Uzz_new, Uxy_new, Uxz_new, Uyz_new;
		 real Const1_new, Const2_new, Const3_new;

#        ifdef UNSPLIT_GRAVITY
         real acc_old[3]={0.0, 0.0, 0.0}, Con_old[NCOMP_FLUID], Pri_old[NCOMP_FLUID];
		 real LorentzFactor_old, n_old , Ux_old, Uy_old, Uz_old, P_old;
		 real Uxx_old, Uyy_old, Uzz_old, Uxy_old, Uxz_old, Uyz_old;
		 real Const1_old, Const2_old, Const3_old;
#        endif

         const int i_g0    = idx_g0 % PS1;
         const int j_g0    = idx_g0 % PS1_sqr / PS1;
         const int k_g0    = idx_g0 / PS1_sqr;

         const int i_new   = i_g0 + GRA_GHOST_SIZE;
         const int j_new   = j_g0 + GRA_GHOST_SIZE;
         const int k_new   = k_g0 + GRA_GHOST_SIZE;
         const int idx_new = IDX321( i_new, j_new, k_new, GRA_NXT, GRA_NXT );

#        ifdef UNSPLIT_GRAVITY
         const int i_old   = i_g0 + USG_GHOST_SIZE;
         const int j_old   = j_g0 + USG_GHOST_SIZE;
         const int k_old   = k_g0 + USG_GHOST_SIZE;
         const int idx_old = IDX321( i_old, j_old, k_old, USG_NXT_G, USG_NXT_G );
#        endif


//       external gravity
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         {
            double x, y, z;

            x = g_Corner_Array[P][0] + (double)(i_g0*dh);
            y = g_Corner_Array[P][1] + (double)(j_g0*dh);
            z = g_Corner_Array[P][2] + (double)(k_g0*dh);

            ExternalAcc( acc_new, x, y, z, TimeNew, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_new[d] *= dt;

#           ifdef UNSPLIT_GRAVITY
            ExternalAcc( acc_old, x, y, z, TimeOld, c_ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    acc_old[d] *= dt;
#           endif
         }


//       self-gravity
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
         {
            const int ip1_new = idx_new + didx_new[0];
            const int jp1_new = idx_new + didx_new[1];
            const int kp1_new = idx_new + didx_new[2];
            const int im1_new = idx_new - didx_new[0];
            const int jm1_new = idx_new - didx_new[1];
            const int km1_new = idx_new - didx_new[2];

#           ifdef UNSPLIT_GRAVITY
            const int ip1_old = idx_old + didx_old[0];
            const int jp1_old = idx_old + didx_old[1];
            const int kp1_old = idx_old + didx_old[2];
            const int im1_old = idx_old - didx_old[0];
            const int jm1_old = idx_old - didx_old[1];
            const int km1_old = idx_old - didx_old[2];
#           endif

            if ( P5_Gradient )
            {
               const real Const_8 = (real)8.0;

               const int ip2_new = ip1_new + didx_new[0];
               const int jp2_new = jp1_new + didx_new[1];
               const int kp2_new = kp1_new + didx_new[2];
               const int im2_new = im1_new - didx_new[0];
               const int jm2_new = jm1_new - didx_new[1];
               const int km2_new = km1_new - didx_new[2];

#              ifdef UNSPLIT_GRAVITY
               const int ip2_old = ip1_old + didx_old[0];
               const int jp2_old = jp1_old + didx_old[1];
               const int kp2_old = kp1_old + didx_old[2];
               const int im2_old = im1_old - didx_old[0];
               const int jm2_old = jm1_old - didx_old[1];
               const int km2_old = km1_old - didx_old[2];
#              endif

               acc_new[0] += Gra_Const*( - pot_new[ip2_new] + Const_8*pot_new[ip1_new] - Const_8*pot_new[im1_new] + pot_new[im2_new] );
               acc_new[1] += Gra_Const*( - pot_new[jp2_new] + Const_8*pot_new[jp1_new] - Const_8*pot_new[jm1_new] + pot_new[jm2_new] );
               acc_new[2] += Gra_Const*( - pot_new[kp2_new] + Const_8*pot_new[kp1_new] - Const_8*pot_new[km1_new] + pot_new[km2_new] );

#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( - pot_old[ip2_old] + Const_8*pot_old[ip1_old] - Const_8*pot_old[im1_old] + pot_old[im2_old] );
               acc_old[1] += Gra_Const*( - pot_old[jp2_old] + Const_8*pot_old[jp1_old] - Const_8*pot_old[jm1_old] + pot_old[jm2_old] );
               acc_old[2] += Gra_Const*( - pot_old[kp2_old] + Const_8*pot_old[kp1_old] - Const_8*pot_old[km1_old] + pot_old[km2_old] );
#              endif
            } // if ( P5_Gradient )

            else
            {
               acc_new[0] += Gra_Const*( pot_new[ip1_new] - pot_new[im1_new] );
               acc_new[1] += Gra_Const*( pot_new[jp1_new] - pot_new[jm1_new] );
               acc_new[2] += Gra_Const*( pot_new[kp1_new] - pot_new[km1_new] );
               
#              ifdef UNSPLIT_GRAVITY
               acc_old[0] += Gra_Const*( pot_old[ip1_old] - pot_old[im1_old] );
               acc_old[1] += Gra_Const*( pot_old[jp1_old] - pot_old[jm1_old] );
               acc_old[2] += Gra_Const*( pot_old[kp1_old] - pot_old[km1_old] );
#              endif
            } // if ( P5_Gradient ) ... else ...
         } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )

         Con_new[DENS] = g_Flu_Array_New[P][DENS][idx_g0];
         Con_new[MOMX] = g_Flu_Array_New[P][MOMX][idx_g0];
         Con_new[MOMY] = g_Flu_Array_New[P][MOMY][idx_g0];
         Con_new[MOMZ] = g_Flu_Array_New[P][MOMZ][idx_g0];
         Con_new[ENGY] = g_Flu_Array_New[P][ENGY][idx_g0];

//       conserved vars --> primitive vars
         LorentzFactor_new = SRHydro_Con2Pri( Con_new, Pri_new, (real)1.333333333, (real)0.0);

         n_new    = Pri_new[0];
         Ux_new   = Pri_new[1];
         Uy_new   = Pri_new[2];
         Uz_new   = Pri_new[3];
         P_new    = Pri_new[4];

         Const1_new = n_new*( (real)2.0 * SQR(LorentzFactor_new)-(real)1.0 );
         Const2_new = ( SQR(Ux_new) + SQR(Uy_new) + SQR(Uz_new) ) * P_new;
         
         Uxx_new = Ux_new*Ux_new;
         Uyy_new = Uy_new*Uy_new;
         Uzz_new = Uz_new*Uz_new;
         Uxy_new = Ux_new*Uy_new;
         Uxz_new = Ux_new*Uz_new;
         Uyz_new = Uy_new*Uz_new;

//       backup temerature
         real Temperature = Pri_new[4]/Pri_new[0];

//       advance fluid
#        ifdef UNSPLIT_GRAVITY
         Con_old[DENS] = g_Flu_Array_USG[P][DENS][idx_g0];
         Con_old[MOMX] = g_Flu_Array_USG[P][MOMX][idx_g0];
         Con_old[MOMY] = g_Flu_Array_USG[P][MOMY][idx_g0];
         Con_old[MOMZ] = g_Flu_Array_USG[P][MOMZ][idx_g0];
         Con_old[ENGY] = g_Flu_Array_USG[P][ENGY][idx_g0];

//       conserved vars --> primitive vars
         LorentzFactor_old = SRHydro_Con2Pri( Con_old, Pri_old, (real)1.333333333, (real)0.0);

         n_old    = Pri_old[0];
         Ux_old   = Pri_old[1];
         Uy_old   = Pri_old[2];
         Uz_old   = Pri_old[3];
         P_old    = Pri_old[4];

         Const1_old = n_old*( (real)2.0 * SQR(LorentzFactor_old)-(real)1.0 );
         Const2_old = ( SQR(Ux_old) + SQR(Uy_old) + SQR(Uz_old) ) * P_old;
         
         Uxx_old = Ux_old*Ux_old;
         Uyy_old = Uy_old*Uy_old;
         Uzz_old = Uz_old*Uz_old;
         Uxy_old = Ux_old*Uy_old;
         Uxz_old = Ux_old*Uz_old;
         Uyz_old = Uy_old*Uz_old;


//       update the momentum density
         Con_new[MOMX] += - Uxx_old * acc_old[0] - Uxy_old * acc_old[1] - Uxz_old * acc_old[2] + Const1_old * acc_old[0] + Const2_old * acc_old[0];
         Con_new[MOMY] += - Uxy_old * acc_old[0] - Uyy_old * acc_old[1] - Uyz_old * acc_old[2] + Const1_old * acc_old[1] + Const2_old * acc_old[1];
         Con_new[MOMZ] += - Uxz_old * acc_old[0] - Uyz_old * acc_old[1] - Uzz_old * acc_old[2] + Const1_old * acc_old[2] + Const2_old * acc_old[2];

         Con_new[MOMX] += - Uxx_new * acc_new[0] - Uxy_new * acc_new[1] - Uxz_new * acc_new[2] + Const1_new * acc_new[0] + Const2_new * acc_new[0];
         Con_new[MOMY] += - Uxy_new * acc_new[0] - Uyy_new * acc_new[1] - Uyz_new * acc_new[2] + Const1_new * acc_new[1] + Const2_new * acc_new[1];
         Con_new[MOMZ] += - Uxz_new * acc_new[0] - Uyz_new * acc_new[1] - Uzz_new * acc_new[2] + Const1_new * acc_new[2] + Const2_new * acc_new[2];

         Con_new[MOMX] *= (real)0.5;
         Con_new[MOMY] *= (real)0.5;
         Con_new[MOMZ] *= (real)0.5;


//       update the total energy density
#        if (  CONSERVED_ENERGY == 1 )
         Const3_old = LorentzFactor_old * ( n_old + P_old );
         Const3_new = LorentzFactor_new * ( n_new + P_new );

         Con_new[ENGY] += Const3_old * ( Ux_old * acc_old[0] + Uy_old * acc_old[1] + Uz_old * acc_old[2] );
         Con_new[ENGY] += Const3_new * ( Ux_new * acc_new[0] + Uy_new * acc_new[1] + Uz_new * acc_new[2] );

         Con_new[ENGY] *= (real)0.5;
#        elif ( CONSERVED_ENERGY == 2 )
#        error: MODIFY __FUNCTION__!
#        else
#        error: CONSERVED_ENERGY must be 1 or 2!
#        endif

#        else  // #ifdef UNSPLIT_GRAVITY

////      1. update the momentum density
         Con_new[MOMX] += - n_new * Uxx_new * acc_new[0] - n_new * Uxy_new * acc_new[1] - n_new * Uxz_new * acc_new[2] + Const1_new * acc_new[0] + Const2_new * acc_new[0];
         Con_new[MOMY] += - n_new * Uxy_new * acc_new[0] - n_new * Uyy_new * acc_new[1] - n_new * Uyz_new * acc_new[2] + Const1_new * acc_new[1] + Const2_new * acc_new[1];
         Con_new[MOMZ] += - n_new * Uxz_new * acc_new[0] - n_new * Uyz_new * acc_new[1] - n_new * Uzz_new * acc_new[2] + Const1_new * acc_new[2] + Const2_new * acc_new[2];
 
//       2. update the momentum density ( remove velocity terms )
//         Con_new[MOMX] +=  n_new  * acc_new[0];
//         Con_new[MOMY] +=  n_new  * acc_new[1];
//         Con_new[MOMZ] +=  n_new  * acc_new[2];


////       3. update the total energy density
//         Const3_new = LorentzFactor_new * ( n_new + P_new );
//         Con_new[ENGY] += Const3_new * ( Ux_new * acc_new[0] + Uy_new * acc_new[1] + Uz_new * acc_new[2] );

//       4. update the total energy density ( assuming temperature is unchanged under gravity )
		 real Msqr = SQR(Con_new[MOMX]) + SQR(Con_new[MOMY]) + SQR(Con_new[MOMZ]);

#        if (  CONSERVED_ENERGY == 1 )
         real h = SpecificEnthalpy( NULL, Temperature, (real)1.333333333 );
		 real Dh = Con_new[DENS]*h;
         real factor = SQRT(Dh*Dh + Msqr);
		 Con_new[ENGY] = factor -  Con_new[DENS] * Dh * Temperature / factor;
#        elif ( CONSERVED_ENERGY == 2 )
		 real Dsqr = SQR(Con_new[DENS]);
         real HTilde = SRHydro_Temperature2HTilde( Temperature );
		 real h = HTilde + (real)1.0;
		 real factor = SQRT( Dsqr*h*h + Msqr );
		 Con_new[ENGY]  = Dsqr * SQR(HTilde) + (real)2.0*HTilde*Dsqr + Msqr;
	     Con_new[ENGY] /= ( factor + Con_new[DENS] );
		 Con_new[ENGY] -= Dsqr*h*Temperature / factor;
#        endif


#        endif // #ifdef UNSPLIT_GRAVITY ... else ...

//       store the updated fluid variables
         g_Flu_Array_New[P][MOMX][idx_g0] = Con_new[MOMX];
         g_Flu_Array_New[P][MOMY][idx_g0] = Con_new[MOMY];
         g_Flu_Array_New[P][MOMZ][idx_g0] = Con_new[MOMZ];
		 g_Flu_Array_New[P][ENGY][idx_g0] = Con_new[ENGY];

#        ifdef CHECK_FAILED_CELL_IN_FLUID
         SRHydro_CheckUnphysical(Con_new, NULL, (real)1.333333333, (real)0.0, __FUNCTION__, __LINE__, true);
#        endif


      } // CGPU_LOOP( idx_g0, CUBE(PS1) )
   } // for (int P=0; P<NPatchGroup*8; P++)

} // FUNCTION : CPU/CUPOT_HydroGravitySolver



#endif // #if ( MODEL == SR_HYDRO  &&  defined GRAVITY )
